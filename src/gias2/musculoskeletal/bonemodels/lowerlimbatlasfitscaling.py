"""
FILE: lowerlimbatlasfitscaling.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: Functions for lower limb atlas landmark registration using only
rigid and scaling transforms. Each model is scaled individually.

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

from gias2.musculoskeletal.bonemodels import modelcore
from gias2.registration import alignment_fitting

log = logging.getLogger(__name__)


# import pdb

def _make_x0(ll, target_landmarks, source_landmarks, init_scalings):
    """ Generate initial parameters for lower limb atlas fitting.
    The pelvis landmarks are rigidly registered to get initial
    rigid transformation parameters
    """
    init_rigid = alignment_fitting.fitRigid(
        source_landmarks[:3, :], target_landmarks[:3, :],
        rotcentre=0.5 * (source_landmarks[0] + source_landmarks[1]),
        xtol=1e-6, maxfev=999999,
        maxfun=9999999999
    )[0]

    init_x = np.hstack([init_scalings,
                        init_rigid,
                        [0, ] * ll.N_PARAMS_HIP,
                        [0, ] * ll.N_PARAMS_KNEE,
                        ])

    return init_x


def _lower_limb_atlas_landmark_fit_multi_scaling(ll_model, target_landmark_coords, landmark_names,
                                                 bones_to_scale=None, x0=None, minimise_args={}):
    if bones_to_scale is None:
        bones_to_scale = ll_model.bone_names
    else:
        for b in bones_to_scale:
            if b not in ll_model.bone_names:
                raise ValueError('{} not a valid bone name, must be one of {}'.format(b, ll_model.bone_names))

    n_models = len(bones_to_scale)  # this number of scaling factors

    if x0 is not None:
        if len(x0) != (ll_model.N_PARAMS_RIGID + n_models):
            raise ValueError('Incorrect number of elements in x0, need {}, given {}'.format(
                ll_model.N_PARAMS_RIGID + n_models, len(x0))
            )
        else:
            scaling_params = [bones_to_scale, x0[:n_models]]
    else:
        scaling_params = [bones_to_scale, [1.0, ] * n_models]

    get_source_landmarks = modelcore.make_source_landmark_getter(landmark_names)
    source_landmark_coords = np.zeros((len(landmark_names), 3), dtype=float)
    source_landmark_coords = get_source_landmarks(ll_model, source_landmark_coords)
    hip_rigid_x_index = n_models  # six elements
    hip_rot_x_index = n_models + ll_model.N_PARAMS_PELVIS  # 3 elements
    knee_rot_x_index = n_models + ll_model.N_PARAMS_PELVIS + ll_model.N_PARAMS_HIP  # 2 elements

    def _x_splitter(x):
        scaling_params[1] = x[:n_models]
        return [scaling_params,
                x[hip_rigid_x_index:hip_rot_x_index],
                x[hip_rot_x_index:knee_rot_x_index],
                x[knee_rot_x_index:],
                ]

    # =========================================================================#
    # define objective function
    def lower_limb_landmark_reg_obj(x):
        """Main objective function
        """
        # update model geometry
        x_split = _x_splitter(x)
        ll_model.update_all_models_multi_scaling(x_split[0],
                                                 x_split[1],
                                                 x_split[2],
                                                 x_split[3],
                                                 )

        # get source landmark coords
        source_x = get_source_landmarks(ll_model, source_landmark_coords)

        # calc sum of squared distance between target and source landmarks
        ssdist = ((target_landmark_coords - source_x) ** 2.0).sum()

        # log.debug(('SSDist: {:12.6f}'.format(ssdist)))
        sys.stdout.write('SSDist: {:12.6f}\r'.format(ssdist))
        sys.stdout.flush()

        return ssdist

    # =========================================================================#
    # make initial parameters
    if x0 is None:
        x0 = _make_x0(ll_model, target_landmark_coords,
                      source_landmark_coords, scaling_params[1])
    else:
        x0 = np.array(x0)

    x_history = [_x_splitter(x0), ]

    # run minimisation
    opt_results = optimize.minimize(lower_limb_landmark_reg_obj,
                                    x0, **minimise_args
                                    )
    xopt_split = _x_splitter(opt_results['x'])
    x_history.append(xopt_split)
    ll_model.update_all_models_multi_scaling(xopt_split[0],
                                             xopt_split[1],
                                             xopt_split[2],
                                             xopt_split[3],
                                             )

    # calc final landmark error
    opt_source_landmarks = get_source_landmarks(ll_model, source_landmark_coords)
    opt_landmark_dist = np.sqrt(((target_landmark_coords - opt_source_landmarks) ** 2.0).sum(1))
    opt_landmark_rmse = np.sqrt((opt_landmark_dist ** 2.0).mean())

    # prepare output
    output_info = {'source_landmark_getter': get_source_landmarks,
                   'obj': lower_limb_landmark_reg_obj,
                   'min_results': opt_results,
                   'opt_source_landmarks': opt_source_landmarks,
                   }

    return x_history, opt_landmark_dist, opt_landmark_rmse, output_info


def _lower_limb_atlas_landmark_fit_uniform_scaling(ll_model, target_landmark_coords, landmark_names,
                                                   x0=None, minimise_args={}):
    n_models = 1  # this number of scaling factors
    if x0 is not None:
        if len(x0) != (ll_model.N_PARAMS_RIGID + n_models):
            raise ValueError('Incorrect number of elements in x0, need {}, given {}'.format(
                ll_model.N_PARAMS_RIGID + n_models, len(x0))
            )
    get_source_landmarks = modelcore.make_source_landmark_getter(landmark_names)
    source_landmark_coords = np.zeros((len(landmark_names), 3), dtype=float)
    source_landmark_coords = get_source_landmarks(ll_model, source_landmark_coords)
    hip_rigid_x_index = n_models
    hip_rot_x_index = n_models + ll_model.N_PARAMS_PELVIS  # 3 elements
    knee_rot_x_index = n_models + ll_model.N_PARAMS_PELVIS + ll_model.N_PARAMS_HIP  # 2 elements

    def _x_splitter(x):
        return [x[0],
                x[hip_rigid_x_index:hip_rot_x_index],
                x[hip_rot_x_index:knee_rot_x_index],
                x[knee_rot_x_index:],
                ]

    # =========================================================================#
    # define objective function
    def lower_limb_landmark_reg_obj(x):
        """Main objective function
        """
        # update model geometry
        x_split = _x_splitter(x)
        ll_model.update_all_models_uniform_scaling(x_split[0],
                                                   x_split[1],
                                                   x_split[2],
                                                   x_split[3],
                                                   )

        # get source landmark coords
        source_x = get_source_landmarks(ll_model, source_landmark_coords)

        # calc sum of squared distance between target and source landmarks
        ssdist = ((target_landmark_coords - source_x) ** 2.0).sum()

        # log.debug(('SSDist: {:12.6f}'.format(ssdist)))
        sys.stdout.write('SSDist: {:12.6f}\r'.format(ssdist))
        sys.stdout.flush()

        return ssdist

    # =========================================================================#
    # make initial parameters
    if x0 is None:
        x0 = _make_x0(ll_model, target_landmark_coords,
                      source_landmark_coords, 1.0)
    else:
        x0 = np.array(x0)

    x_history = [_x_splitter(x0), ]

    # run minimisation
    opt_results = optimize.minimize(lower_limb_landmark_reg_obj,
                                    x0, **minimise_args
                                    )
    xopt_split = _x_splitter(opt_results['x'])
    x_history.append(xopt_split)
    ll_model.update_all_models_uniform_scaling(xopt_split[0],
                                               xopt_split[1],
                                               xopt_split[2],
                                               xopt_split[3],
                                               )

    # calc final landmark error
    opt_source_landmarks = get_source_landmarks(ll_model, source_landmark_coords)
    opt_landmark_dist = np.sqrt(((target_landmark_coords - opt_source_landmarks) ** 2.0).sum(1))
    opt_landmark_rmse = np.sqrt((opt_landmark_dist ** 2.0).mean())

    # prepare output
    output_info = {'source_landmark_getter': get_source_landmarks,
                   'obj': lower_limb_landmark_reg_obj,
                   'min_results': opt_results,
                   'opt_source_landmarks': opt_source_landmarks,
                   }

    return x_history, opt_landmark_dist, opt_landmark_rmse, output_info


def fit(ll_model, target_landmark_coords, landmark_names,
        bones_to_scale='uniform', x0=None,
        minimise_args={}, verbose=False):
    """Fit a lower limb atlas model to landmarks. Only rigid body
    transformations and isotropic scaling for each bone..

    Inputs:
    ll_model: LowerLimbLeftAtlas instance
    target_landmark_coords [m*3 array]: coordinates of m landmarks
    landmark_names [list of strings]: names of landmarks. Refer to 
        gias.musculoskeletal.fw_model_landmarks
    bones_to_scale: if 'uniform' - uniform scaling for all bones.
                    if list of bone names - non-uniform scaling for each named bone.
                    if lists of the 2 above - multi-stage reg 
    x0 [list]: Optional. initial fitting parameters for the 1st fit.
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

    # parse bones_to_scale input
    if isinstance(bones_to_scale, (tuple, list)):
        # if a list of lists or uniform then other elements
        if isinstance(bones_to_scale[0], (list, tuple)) or bones_to_scale[0] == 'uniform':
            multi_fit = True
            n_iterations = len(bones_to_scale)
            # prepare minimise args
            if isinstance(minimise_args, (list, tuple)):
                if len(minimise_args) != n_iterations:
                    raise ValueError('Length of minimise_args and initial_scaling do not match')
            else:
                minimise_args = [minimise_args, ] * n_iterations
        elif isinstance(bones_to_scale[0], str):
            # a list of bone names: nonuniform 1 iteration registration
            uniform_scale = False
            multi_fit = False
            n_iterations = 1
            if not isinstance(minimise_args, dict):
                raise ValueError('minimise_args must be a dictionary')
        else:
            raise ValueError('Unknow input for bones_to_scale: {}'.format(bones_to_scale))
    elif bones_to_scale == 'uniform':
        # uniform 1 interation registration
        uniform_scale = True
        multi_fit = False
        n_iterations = 1
        if not isinstance(minimise_args, dict):
            raise ValueError('minimise_args must be a dictionary')
    else:
        raise ValueError('Unknow input for bones_to_scale: {}'.format(bones_to_scale))

    total_rigid_params = ll_model.N_PARAMS_RIGID

    # Run fitting
    if not multi_fit:
        # run single fit
        if verbose:
            log.debug('Running single lower limb fit')

        if uniform_scale:
            return _lower_limb_atlas_landmark_fit_uniform_scaling(
                ll_model, target_landmark_coords, landmark_names,
                x0=x0, minimise_args=minimise_args)
        else:
            return _lower_limb_atlas_landmark_fit_multi_scaling(
                ll_model, target_landmark_coords, landmark_names,
                bones_to_scale=bones_to_scale,
                x0=x0, minimise_args=minimise_args)
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
            if bones_to_scale[it] == 'uniform':
                uniform_scale = True
            else:
                uniform_scale = False

            if it > 0:
                # prepare new x0 from previous results

                # if previous was a uniform fit and now is another uniform
                if uniform_scale and previous_uniform_scale:
                    x0 = np.hstack(x_hist_it[-1])
                # if previous was a uniform fit and now is a multi scale fit
                elif (not uniform_scale) and previous_uniform_scale:
                    x0 = np.hstack([
                        [x_hist_it[-1][0], ] * len(bones_to_scale[it]),
                        np.hstack(x_hist_it[-1][1:]),
                    ])
                    # pdb.set_trace()
                # if previous and current are both multi scale fits
                elif (not uniform_scale) and (not previous_uniform_scale):

                    new_s = []
                    for b in bones_to_scale[it]:
                        if b in bones_to_scale[it - 1]:
                            # pdb.set_trace()
                            new_s.append(x_hist_it[-1][0][1][bones_to_scale[it - 1].index(b)])
                        else:
                            new_s.append(1.0)

                    x0 = np.hstack([
                        new_s,
                        np.hstack(x_hist_it[-1][1:]),
                    ])
                # illegal case
                else:
                    # elif uniform_scale and (not previous_uniform_scale):
                    raise ValueError('Uniform fit cannot follow multiscale fit')

            if verbose:
                log.debug(('it: {}'.format(it + 1)))
                log.debug(('scaling: {}'.format(bones_to_scale[it])))
                log.debug(('minargs: {}'.format(minimise_args[it])))
                log.debug(('x0: {}'.format(x0)))

            if uniform_scale:
                previous_uniform_scale = True
                x_hist_it, opt_dist_it, \
                opt_rmse_it, info_it = _lower_limb_atlas_landmark_fit_uniform_scaling(
                    ll_model, target_landmark_coords, landmark_names,
                    x0=x0, minimise_args=minimise_args[it])
            else:
                previous_uniform_scale = False
                x_hist_it, opt_dist_it, \
                opt_rmse_it, info_it = _lower_limb_atlas_landmark_fit_multi_scaling(
                    ll_model, target_landmark_coords, landmark_names,
                    bones_to_scale=bones_to_scale[it],
                    x0=x0, minimise_args=minimise_args[it])

            x_history.append(x_hist_it[-1])
            opt_landmark_dist.append(opt_dist_it)
            opt_landmark_rmse.append(opt_rmse_it)
            output_info.append(info_it)

            if verbose:
                log.debug(('it: {}, landmark rmse: {}'.format(it + 1, opt_rmse_it)))

    return x_history, opt_landmark_dist, opt_landmark_rmse, output_info
