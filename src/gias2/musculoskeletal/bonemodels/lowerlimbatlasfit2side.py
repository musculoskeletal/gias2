"""
FILE: lowerlimbatlasfit.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: Functions for 2-side lower limb atlas landmark registration using
PCA fitting

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""
import logging

import numpy as np
import sys
from scipy import optimize

from gias2.registration import alignment_fitting

log = logging.getLogger(__name__)


def _make_x0(ll, npcs, landmark_names, target_landmarks, source_landmarks, init_pc_weights=None):
    """ Generate initial parameters for lower limb atlas fitting.
    The pelvis landmarks are rigidly registered to get initial
    rigid transformation parameters
    """
    reg_landmarks = ('pelvis-LASIS', 'pelvis-RASIS', 'pelvis-Sacral')

    # check if we can do an auto registration of the pelvis
    do_rigid = True
    if 'pelvis-Sacral' not in landmark_names:
        # calculate sacral from LPSIS and RPSIS
        if ('pelvis-LPSIS' in landmark_names) and ('pelvis-RPSIS' in landmark_names):
            landmark_names.append('pelvis-Sacral')
            targ_lpsis = target_landmarks[landmark_names.index('pelvis-LPSIS')]
            targ_rpsis = target_landmarks[landmark_names.index('pelvis-RPSIS')]
            targ_sac = (np.array(targ_lpsis) + np.array(targ_rpsis)) * 0.5
            target_landmarks = np.vstack([target_landmarks, targ_sac])
            src_lpsis = target_landmarks[landmark_names.index('pelvis-LPSIS')]
            src_rpsis = target_landmarks[landmark_names.index('pelvis-RPSIS')]
            src_sac = (np.array(src_lpsis) + np.array(src_rpsis)) * 0.5
            source_landmarks = np.vstack([source_landmarks, src_sac])

    for name in reg_landmarks:
        if name not in landmark_names:
            log.debug('Could not find {} pelvis landmarks. Using zero initial rigid transformations'.format(name))
            do_rigid = False

    if do_rigid:
        # get pelvis landmark coords
        _target = np.array([target_landmarks[landmark_names.index(n), :] for n in reg_landmarks])
        _source = np.array([source_landmarks[landmark_names.index(n), :] for n in reg_landmarks])

        rx, ry, rz = [0, 0, 0]
        t0 = np.hstack([
            _target.mean(0) - _source.mean(0),
            [rx, ry, rz],
        ])
        init_rigid, fitted_landmarks, fit_errors = alignment_fitting.fitRigid(
            _source, _target,
            t0=t0,
            rotcentre=0.5 * (_source[0] + _source[1]),
            xtol=1e-9, maxfev=1e6,
            maxfun=1e6,
            epsfcn=1e-6,
            outputErrors=1,
        )
        log.debug('init rigid transform: {}'.format(init_rigid))
        log.debug('init rigid fit error {} -> {}'.format(*fit_errors))
        # fitted_errors = np.sqrt(((_target - fitted_landmarks)**2.0).sum(1))
        # log.debug('fitted errors: {}'.format(fitted_errors))
        # log.debug('fitted ldmks: {}'.format(fitted_landmarks))
    else:
        init_rigid = np.zeros(6, dtype=float)

    if init_pc_weights is None:
        init_pc_weights = np.zeros(npcs, dtype=float)

    init_x = np.hstack([init_pc_weights,
                        init_rigid,
                        [0, ] * ll.N_PARAMS_HIP_L,
                        [0, ] * ll.N_PARAMS_HIP_R,
                        [0, ] * ll.N_PARAMS_KNEE_L,
                        [0, ] * ll.N_PARAMS_KNEE_R,
                        ])

    return init_x


def _make_source_landmark_getter(landmark_names):
    """ Creates a function to return the coordinates of landmarks
    in a multibone atlas.

    inputs
    ------
    landmark_names : list of str
        List of landmark names, should be suffixed with '-l' or '-r'
        for side. e.g. ['femur-HC-l', 'femur-MEC-r']
    """
    landmark_labels = []
    for ln in landmark_names:
        terms = ln.split('-')
        if len(terms) == 2:
            body_name, landmark_name = terms
            landmark_labels.append((body_name, ln))
        elif len(terms) == 3:
            body_name, landmark_name, side = terms
            landmark_labels.append(('{}-{}'.format(body_name, side), '{}-{}'.format(body_name, landmark_name)))

    def get_source_landmarks(model, coords):
        for i, (bn, ln) in enumerate(landmark_labels):
            coords[i] = model.models[bn].landmarks[ln]

        return coords

    return get_source_landmarks


def _lower_limb_atlas_landmark_fit(ll_model, target_landmark_coords, landmark_names,
                                   pc_modes, mweight, initial_pc_weights=None, x0=None,
                                   callback=None, minimise_args={}):
    if initial_pc_weights is not None:
        if len(initial_pc_weights) != len(pc_modes):
            raise ValueError('Length of initial_pc_weights not equal to length of pc_modes')

    if x0 is not None:
        if len(x0) != (ll_model.N_PARAMS_RIGID + len(pc_modes)):
            raise ValueError('Incorrect number of elements in x0, need {}, given {}'.format(
                ll_model.N_PARAMS_RIGID + len(pc_modes), len(x0))
            )

    n_pc_modes = len(pc_modes)
    get_source_landmarks = _make_source_landmark_getter(landmark_names)
    source_landmark_coords = np.zeros((len(landmark_names), 3), dtype=float)
    source_landmark_coords = get_source_landmarks(ll_model, source_landmark_coords)
    hip_rigid_x_index = n_pc_modes
    hip_rot_l_x_index = n_pc_modes + ll_model.N_PARAMS_PELVIS
    hip_rot_r_x_index = n_pc_modes + ll_model.N_PARAMS_PELVIS + \
                        ll_model.N_PARAMS_HIP_L
    knee_rot_l_x_index = n_pc_modes + ll_model.N_PARAMS_PELVIS + \
                         ll_model.N_PARAMS_HIP_L + ll_model.N_PARAMS_HIP_R
    knee_rot_r_x_index = n_pc_modes + ll_model.N_PARAMS_PELVIS + \
                         ll_model.N_PARAMS_HIP_L + ll_model.N_PARAMS_HIP_R + \
                         ll_model.N_PARAMS_KNEE_L

    x_history = []

    def _x_splitter(x):
        return [x[:n_pc_modes],
                x[hip_rigid_x_index:hip_rot_l_x_index],
                x[hip_rot_l_x_index:hip_rot_r_x_index],
                x[hip_rot_r_x_index:knee_rot_l_x_index],
                x[knee_rot_l_x_index:knee_rot_r_x_index],
                x[knee_rot_r_x_index:],
                ]

    # =========================================================================#
    # define objective function
    def lower_limb_landmark_reg_obj(x):
        """Main objective function
        """
        # update model geometry
        x_split = _x_splitter(x)
        ll_model.update_all_models(x_split[0], pc_modes,
                                   x_split[1],
                                   x_split[2],
                                   x_split[3],
                                   x_split[4],
                                   x_split[5],
                                   )

        # get source landmark coords
        source_x = get_source_landmarks(ll_model, source_landmark_coords)

        # calc sum of squared distance between target and source landmarks
        ssdist = ((target_landmark_coords - source_x) ** 2.0).sum()

        # calc squared mahalanobis distance
        m2 = mweight * (x_split[0] ** 2.0).sum()
        sys.stdout.write('SSDist: {:12.6f}, mdistance: {:8.5f}\r'.format(ssdist, m2))
        sys.stdout.flush()
        # log.debug('SSDist: {:12.6f}, mdistance: {:8.5f}'.format(ssdist, m2))

        # total error
        total_error = ssdist + m2

        return total_error

    # def default_callback(x):
    #     # update model geometry
    #     x_split = _x_splitter(x)
    #     ll_model.update_all_models(x_split[0], pc_modes, 
    #                                x_split[1],
    #                                x_split[2],
    #                                x_split[3],
    #                                )

    #     # get source landmark coords
    #     source_x = get_source_landmarks(ll_model, source_landmark_coords)

    #     # calc sum of squared distance between target and source landmarks
    #     ssdist = ((target_landmark_coords - source_x)**2.0).sum()

    #     # calc squared mahalanobis distance
    #     m2 = mweight * (x_split[0]**2.0).sum()
    #     log.debug('SSDist: {:12.6f}, mdistance: {:8.5f}'.format(ssdist, m2))

    # =========================================================================#
    # make initial parameters
    if x0 is None:
        x0 = _make_x0(ll_model, n_pc_modes, landmark_names, target_landmark_coords,
                      source_landmark_coords, initial_pc_weights)
    else:
        x0 = np.array(x0)

    x_history.append(_x_splitter(x0))

    # if callback is None:
    #     callback = default_callback

    # run minimisation
    log.debug(' ')
    opt_results_1 = optimize.minimize(lower_limb_landmark_reg_obj,
                                      x0, callback=callback,
                                      **minimise_args
                                      )
    log.debug(' ')
    xopt1_split = _x_splitter(opt_results_1['x'])
    x_history.append(xopt1_split)
    ll_model.update_all_models(xopt1_split[0], pc_modes,
                               xopt1_split[1],
                               xopt1_split[2],
                               xopt1_split[3],
                               xopt1_split[4],
                               xopt1_split[5],
                               )

    # calc final landmark error
    opt_source_landmarks = get_source_landmarks(ll_model, source_landmark_coords)
    opt_landmark_dist = np.sqrt(((target_landmark_coords - opt_source_landmarks) ** 2.0).sum(1))
    opt_landmark_rmse = np.sqrt((opt_landmark_dist ** 2.0).mean())

    # prepare output
    output_info = {'source_landmark_getter': get_source_landmarks,
                   'obj': lower_limb_landmark_reg_obj,
                   'min_results': opt_results_1,
                   'opt_source_landmarks': opt_source_landmarks,
                   'mahalanobis_distance': np.sqrt((xopt1_split[0] ** 2.0).sum()),
                   }

    return x_history, opt_landmark_dist, opt_landmark_rmse, output_info


def fit(ll_model, target_landmark_coords, landmark_names,
        pc_modes, mweight, initial_pc_weights=None, x0=None,
        callback=None, minimise_args={}, verbose=False):
    """Fit a lower limb atlas model to landmarks.

    Inputs:
    ll_model: LowerLimbAtlas instance
    target_landmark_coords [m*3 array]: coordinates of m landmarks
    landmark_names [list of strings]: names of landmarks suffixed with 'l' or
        'r'. Refer to gias2.musculoskeletal.fw_model_landmarks
    pc_modes [list of ints]: the principal components to use in the fitting.
        If list is nested, fits will be successively carried out using
        the modes defined in each nested list, and the fitted parameters in
        the previous fit passed as initial parameters to the next fit.
    mweight [float or list]: mahalanobis distance penalty weight. If a list
        if given, each float should be the mweight corresponding to a fit
        defined in pc_modes.
    initial_pc_weights [list]: optional. Initial weights for each mode in the
        first fit.
    x0 [list]: initial fitting parameters for the 1st fit.
    minimise_args [dictionary or list of dictionaries]: keyword arguments for
        scipy.optimize.minimise. A set can be provide for each fit if a list
        of dicts is given.
    verbose [bool]: print fitting progress

    Returns:
    Each returned variable will be in list form if there are multiple fits.

    x_history [1-d array]: list of fitted parameters through each fit
    opt_landmark_dist [1-d array]: fitted distance to each target landmark
    opt_landmark_rmse [float]: fitted RMS distance to landmarks
    output_info [list]: list of dictionaries containing detailed information
        from each fit. Keywords:
            source_landmark_getter: function for getting model landmark
                coordinates
            obj: the objective function
            min_results: full output by scipy.optimize.minimize
            opt_source_landmarks: fitted model landmark coordinates

    """

    if len(target_landmark_coords) != len(landmark_names):
        raise ValueError('Number of target landmarks not equal to number of landmark names')

    if isinstance(pc_modes[0], (list, tuple, np.ndarray)):
        multi_fit = True
        n_iterations = len(pc_modes)

        # prepare mweights
        if isinstance(mweight, (list, tuple, np.ndarray)):
            if len(mweight) != n_iterations:
                raise ValueError('Length of mweight and pc_modes do not match')
        else:
            mweight = [mweight, ] * n_iterations

        # prepare minimise args
        if isinstance(minimise_args, (list, tuple)):
            if len(minimise_args) != n_iterations:
                raise ValueError('Length of minimise_args and pc_modes do not match')
        else:
            minimise_args = [minimise_args, ] * n_iterations
    else:
        multi_fit = False

    total_rigid_params = ll_model.N_PARAMS_RIGID

    # Run fitting
    if not multi_fit:
        # run single fit
        if verbose:
            log.debug('Running single lower limb fit')

        return _lower_limb_atlas_landmark_fit(
            ll_model, target_landmark_coords, landmark_names, pc_modes,
            mweight, initial_pc_weights=initial_pc_weights,
            x0=x0, callback=callback, minimise_args=minimise_args)
    else:
        # run multi-stage fit
        if verbose:
            log.debug(('Running {}-stage lower limb fit'.format(n_iterations)))

        x_history = []
        opt_landmark_dist = []
        opt_landmark_rmse = []
        output_info = []
        if x0 is not None:
            x_history.append(x0)

        for it in range(n_iterations):
            if it > 0:
                initial_pc_weights = None
                # check if new x0 is of right length
                n_modes_diff = len(pc_modes[it]) - len(pc_modes[it - 1])
                if n_modes_diff > 0:
                    # pad
                    x0 = np.hstack([
                        np.hstack([
                            x0[:-total_rigid_params],
                            np.zeros(n_modes_diff, dtype=float)
                        ]),
                        x0[-total_rigid_params:]
                    ])
                elif n_modes_diff < 0:
                    # truncate
                    x0 = np.hstack([
                        x0[:len(pc_modes)],
                        x0[-total_rigid_params:]
                    ])

            if verbose:
                log.debug(('it: {}'.format(it + 1)))
                log.debug(('modes: {}'.format(pc_modes[it])))
                log.debug(('mweight: {}'.format(mweight[it])))
                log.debug(('minargs: {}'.format(minimise_args[it])))
                log.debug(('x0: {}'.format(x0)))

            x_hist_it, opt_dist_it, opt_rmse_it, info_it = _lower_limb_atlas_landmark_fit(
                ll_model, target_landmark_coords, landmark_names,
                pc_modes[it], mweight[it],
                initial_pc_weights=initial_pc_weights,
                x0=x0, callback=callback,
                minimise_args=minimise_args[it])

            x0 = np.hstack(x_hist_it[-1])

            x_history.append(x_hist_it[-1])
            opt_landmark_dist.append(opt_dist_it)
            opt_landmark_rmse.append(opt_rmse_it)
            output_info.append(info_it)

            if verbose:
                log.debug(('it: {}, landmark rmse: {}'.format(it + 1, opt_rmse_it)))

    return x_history, opt_landmark_dist, opt_landmark_rmse, output_info
