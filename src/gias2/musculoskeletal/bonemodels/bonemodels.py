"""
FILE: bonemodels.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: Classes and functions for specific bone models

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""
from os import path

import numpy as np
from scipy.spatial import cKDTree

from gias2.common import transform3D, math
from gias2.fieldwork.field import geometric_field
from gias2.musculoskeletal import model_alignment
from gias2.musculoskeletal.bonemodels import modelcore


# =======================================#
# Individual bone models                #
# =======================================#
class PelvisModel(modelcore.BoneModel):
    model_landmarks = [
        'pelvis-LASIS', 'pelvis-RASIS',
        'pelvis-Sacral',
        'pelvis-LHJC', 'pelvis-RHJC',
        'pelvis-LPSIS', 'pelvis-RPSIS',
        'pelvis-SacPlat',
    ]

    def __init__(self, name, gf, side=None):
        super(PelvisModel, self).__init__(name, gf)
        self.init_landmarks(self.model_landmarks)

    def update_acs(self):
        self.acs.update(*model_alignment.createPelvisACSISB(
            self.landmarks['pelvis-LASIS'],
            self.landmarks['pelvis-RASIS'],
            self.landmarks['pelvis-LPSIS'],
            self.landmarks['pelvis-RPSIS'],
        )
                        )


class FemurModel(modelcore.BoneModel):
    KNEE_ELEMS = [46, 43, 44, 47, 26, 27, 28, 29, 30, 52, 51, 48, 49]
    model_landmarks = [
        'femur-HC', 'femur-MEC',
        'femur-LEC', 'femur-GT',
    ]

    def __init__(self, name, gf, knee_surf_disc=4.0, side='left'):
        super(FemurModel, self).__init__(name, gf)
        self.side = side
        self.init_landmarks(self.model_landmarks, side=self.side)
        self.knee_surf_disc = knee_surf_disc
        self.knee_surf_evaluator = geometric_field.makeGeometricFieldElementsEvaluatorSparse(
            self.gf, self.KNEE_ELEMS,
            self.knee_surf_disc)

    def update_acs(self):
        self.acs.update(*model_alignment.createFemurACSISB(
            self.landmarks['femur-HC'],
            self.landmarks['femur-MEC'],
            self.landmarks['femur-LEC'],
            side=self.side
        )
                        )

    def evaluate_knee_surface(self):
        return self.knee_surf_evaluator(self.gf.get_field_parameters()).T


class TibiaFibulaModel(modelcore.BoneModel):
    KNEE_ELEMS = [37, 38, 39, 40, 41, 42, 43, 44, 45]
    model_landmarks = [
        'tibiafibula-LC', 'tibiafibula-MC',
        'tibiafibula-LM', 'tibiafibula-MM',
        'tibiafibula-TT',
    ]

    def __init__(self, name, gf, knee_surf_disc=4.0, side='left'):
        super(TibiaFibulaModel, self).__init__(name, gf)
        self.side = side
        self.init_landmarks(self.model_landmarks, side=self.side)
        self.knee_surf_disc = knee_surf_disc
        self.knee_surf_evaluator = geometric_field.makeGeometricFieldElementsEvaluatorSparse(
            self.gf, self.KNEE_ELEMS,
            self.knee_surf_disc)

    def update_acs(self):
        self.acs.update(*model_alignment.createTibiaFibulaACSISB(
            self.landmarks['tibiafibula-MM'],
            self.landmarks['tibiafibula-LM'],
            self.landmarks['tibiafibula-MC'],
            self.landmarks['tibiafibula-LC'],
            side=self.side
        )
                        )

    def evaluate_knee_surface(self):
        return self.knee_surf_evaluator(self.gf.get_field_parameters()).T


class PatellaModel(modelcore.BoneModel):
    model_landmarks = [
        'patella-sup', 'patella-inf', 'patella-lat'
    ]

    def __init__(self, name, gf, side='left'):
        super(PatellaModel, self).__init__(name, gf)
        self.side = side
        self.init_landmarks(
            self.model_landmarks,
            side=self.side
        )

    def update_acs(self):
        self.acs.update(*model_alignment.createPatellaACSTest(
            self.landmarks['patella-sup'],
            self.landmarks['patella-inf'],
            self.landmarks['patella-lat'],
            side=self.side
        )
                        )


# ==================================================#
# Multiple bone models                             #
# ==================================================#
MAP_DIR = path.expanduser('~/work_projects/MAP/')


class LowerLimbLeftAtlas(modelcore.MultiBoneAtlas):
    bone_names = ('pelvis', 'femur', 'patella', 'tibiafibula')
    combined_pcs_filename = path.join(MAP_DIR, 'models/lower_limb/shapemodels/LLP26_rigid.pc')
    bone_classes = {'pelvis': PelvisModel,
                    'femur': FemurModel,
                    'patella': PatellaModel,
                    'tibiafibula': TibiaFibulaModel,
                    }
    bone_files = {
        'pelvis': (path.join(MAP_DIR, 'models/pelvis/template_meshes/pelvis_combined_cubic_mean_rigid_LLP26.geof'),
                   path.join(MAP_DIR, 'models/pelvis/template_meshes/pelvis_combined_cubic_flat.ens'),
                   path.join(MAP_DIR, 'models/pelvis/template_meshes/pelvis_combined_cubic_flat.mesh'),
                   ),
        'femur': (path.join(MAP_DIR, 'models/femur/template_meshes/femur_left_mean_rigid_LLP26.geof'),
                  path.join(MAP_DIR, 'models/femur/template_meshes/femur_left_quartic_flat.ens'),
                  path.join(MAP_DIR, 'models/femur/template_meshes/femur_left_quartic_flat.mesh'),
                  ),
        'patella': (path.join(MAP_DIR, 'models/patella/template_meshes/patella_left_mean_rigid_LLP26.geof'),
                    path.join(MAP_DIR, 'models/patella/template_meshes/patella_11_left.ens'),
                    path.join(MAP_DIR, 'models/patella/template_meshes/patella_11_left.mesh'),
                    ),
        'tibiafibula': (
        path.join(MAP_DIR, 'models/tibiafibula/template_meshes/tibia_fibula_cubic_left_mean_rigid_LLP26.geof'),
        path.join(MAP_DIR, 'models/tibiafibula/template_meshes/tibia_fibula_left_cubic_flat.ens'),
        path.join(MAP_DIR, 'models/tibiafibula/template_meshes/tibia_fibula_left_cubic_flat.mesh'),
        ),
        }
    combined_model_field_basis = {'tri10': 'simplex_L3_L3',
                                  'tri15': 'simplex_L4_L4',
                                  'quad44': 'quad_L3_L3',
                                  'quad55': 'quad_L4_L4',
                                  'quad54': 'quad_L4_L3',
                                  }

    HJC = 'pelvis-LHJC'
    side = 'left'

    tp_lat_node = 265
    tp_med_node = 293
    default_knee_shift = -20.0
    knee_gap = 5.0
    patella_shift = np.array([50.0, 50.0, -10.0])
    _allow_knee_adduction_dof = False
    _allow_knee_adduction_correction = False
    _neutral_params = [[0, ], [0, ], [0, 0, 0, 0, 0, 0], [0, 0, 0], [0, ]]
    N_PARAMS_PELVIS = 6
    N_PARAMS_HIP = 3
    N_PARAMS_KNEE = 1
    N_PARAMS_RIGID = N_PARAMS_PELVIS + N_PARAMS_HIP + N_PARAMS_KNEE

    hip_flex_coeff = 1.0
    hip_abd_coeff = -1.0
    hip_rot_coeff = -1.0
    knee_flex_coeff = 1.0
    knee_abd_coeff = -1.0
    knee_rot_coeff = 1.0

    def load_bones(self):
        """ Load atlas bone models and set to neutral positions
        """
        self.load_models(self.bone_names, self.bone_classes, self.bone_files)
        self.load_combined_pcs(self.combined_pcs_filename)

        self.models['femur'].side = self.side
        self.models['femur'].init_landmarks(
            self.models['femur'].model_landmarks, side=self.side
        )
        self.models['femur'].update_acs()

        self.models['patella'].side = self.side
        self.models['patella'].init_landmarks(
            self.models['patella'].model_landmarks, side=self.side
        )
        self.models['patella'].update_acs()

        self.models['tibiafibula'].side = self.side
        self.models['tibiafibula'].init_landmarks(
            self.models['tibiafibula'].model_landmarks, side=self.side
        )
        self.models['tibiafibula'].update_acs()

        # zero shape model
        self.update_models_by_pcweights_sd([0, ], [0, ])
        # get pelvis origin from which rigid transformation will be made
        self.models['pelvis'].update_landmarks()
        self.models['pelvis'].update_acs()
        self.pelvis_origin = np.array(self.models['pelvis'].acs.o)
        # zero model again, this time with angles
        self.update_all_models(*self._neutral_params)

    def set_bone_gfield(self, name, gf):
        """
        Set the geometric field of a bone

        Input
        -----
        name : string
            Name of the bone. Must be one of "pelvis", "femur", "patella",
            "tibiafibula".
        gf : geometric_field instance
            geometric_field of the bone
        side : str
            one of 'left' or 'right'
        """

        if name not in self.bone_classes:
            raise ValueError('Invalid bone name. Must be one of {}'.format(str(list(self.bone_classes.keys()))))

        if not self.models.get(name):
            self.models[name] = self.bone_classes[name](name, gf)
        else:
            self.models[name].update_gf(gf.field_parameters.copy())

        self.models[name].side = self.side
        self.models[name].init_landmarks(
            self.models[name].model_landmarks, side=self.side
        )
        self.models[name].update_acs()

    def _set_kneegap_function(self, f):
        self._reset_tibia_kneegap = f

    def enable_knee_adduction_correction(self):
        self._allow_knee_adduction_correction = True
        self._set_kneegap_function(self._reset_tibia_kneegap_2)

    def disable_knee_adduction_correction(self):
        self._allow_knee_adduction_correction = False
        self._set_kneegap_function(self._reset_tibia_kneegap_1)

    def enable_knee_adduction_dof(self):
        self._allow_knee_adduction_dof = True
        self._neutral_params = [[0, ], [0, ], [0, 0, 0, 0, 0, 0], [0, 0, 0], [0, 0]]
        self.N_PARAMS_KNEE = 2
        self.N_PARAMS_RIGID = self.N_PARAMS_PELVIS + self.N_PARAMS_HIP + self.N_PARAMS_KNEE

    def disable_knee_adduction_dof(self):
        self._allow_knee_adduction_dof = False
        self._neutral_params = [[0, ], [0, ], [0, 0, 0, 0, 0, 0], [0, 0, 0], [0, ]]
        self.N_PARAMS_KNEE = 1
        self.N_PARAMS_RIGID = self.N_PARAMS_PELVIS + self.N_PARAMS_HIP + self.N_PARAMS_KNEE

    def _reset_femur_hip(self):
        """Reset femur position to have 0 rotations at the hip.
        i.e. align femur ACS with pelvis ACS at the LHJC
        """
        # translate axes to femur system
        op = self.models['pelvis'].landmarks[self.HJC]
        cs_targ = np.array([op,
                            op + self.models['pelvis'].acs.x,
                            op + self.models['pelvis'].acs.y,
                            op + self.models['pelvis'].acs.z,
                            ])

        of = self.models['femur'].landmarks['femur-HC']
        cs_source = np.array([of,
                              of + self.models['femur'].acs.x,
                              of + self.models['femur'].acs.y,
                              of + self.models['femur'].acs.z,
                              ])

        T = transform3D.directAffine(cs_source, cs_targ)
        self.models['femur'].gf.transformAffine(T)
        self.models['femur'].update_landmarks()
        self.models['femur'].update_acs()

    def _reset_tibia_knee(self):
        """Reset tibia position to have 0 rotations at the knee.
        i.e. align tibia ACS with femur ACS at the femur origin
        """
        # translate axes to femur system
        of = self.models['femur'].acs.o
        cs_targ = np.array([of,
                            of + self.models['femur'].acs.x,
                            of + self.models['femur'].acs.y,
                            of + self.models['femur'].acs.z,
                            ])

        ot = 0.5 * (self.models['tibiafibula'].landmarks['tibiafibula-MC'] +
                    self.models['tibiafibula'].landmarks['tibiafibula-LC'])
        cs_source = np.array([ot,
                              ot + self.models['tibiafibula'].acs.x,
                              ot + self.models['tibiafibula'].acs.y,
                              ot + self.models['tibiafibula'].acs.z,
                              ])

        T = transform3D.directAffine(cs_source, cs_targ)
        self.models['tibiafibula'].gf.transformAffine(T)
        self.models['tibiafibula'].update_landmarks()
        self.models['tibiafibula'].update_acs()

        self._reset_tibia_kneegap_1()

    def _calc_knee_gap(self):
        """Calculate the closest distance between points on the femoral
        knee articular surface and the tibial knee articular surface.
        """

        # evaluate points
        femur_points = self.models['femur'].evaluate_knee_surface()
        tibia_points = self.models['tibiafibula'].evaluate_knee_surface()

        # find closest neighbours
        femur_tree = cKDTree(femur_points)
        dist, femur_points_i = femur_tree.query(tibia_points)

        # calc distance in the tibial Y for each pair
        tibia_y_dist = np.dot(femur_points[femur_points_i] - tibia_points,
                              self.models['tibiafibula'].acs.y
                              )
        return tibia_y_dist.min()

    def _calc_knee_gap_2(self, tp_plane, femur_points):
        D = tp_plane.calcDistanceToPlane(femur_points)
        return D.min()

    def _reset_tibia_kneegap_1(self):
        """shift tibia along tibia y to maintain knee joint gap
        """
        current_knee_gap = self._calc_knee_gap()
        knee_shift = current_knee_gap - self.knee_gap  # if current_knee_gap > self.knee_gap, move in +ve tibia y to close the gap
        shift_t = knee_shift * self.models['tibiafibula'].acs.y
        self.models['tibiafibula'].gf.transformTranslate(shift_t)
        self.models['tibiafibula'].update_landmarks()
        self.models['tibiafibula'].update_acs()

    def _reset_tibia_kneegap_2(self):
        """ shift and reorient tibia to fit to femoral condyles. Allow rotation
        in its x axis
        """

        # place tibia as per normal
        self._reset_tibia_kneegap_1()

        # build kdtree of femoral condyle points
        femur_points = self.models['femur'].evaluate_knee_surface()
        femur_tree = cKDTree(femur_points)

        # calculate and apply varus-valgus angle (about floating-x)
        varus_angle = self._calc_varus_angle(femur_tree)[0]
        floating_x = self._get_knee_floating_x()
        self.models['tibiafibula'].gf.transformRotateAboutAxis(
            -varus_angle,
            self.models['femur'].acs.o,
            self.models['femur'].acs.o + floating_x,
        )
        self.models['tibiafibula'].update_landmarks()
        self.models['tibiafibula'].update_acs()

        # shift tibia to maintain defined knee joint gap
        # shift in the direction normal tp_v and floating x
        # tp_lat = self.models['tibiafibula'].gf.get_point_position(self.tp_lat_node)
        # tp_med = self.models['tibiafibula'].gf.get_point_position(self.tp_med_node)
        # tp_v = tp_lat - tp_med
        # knee_gap_v = math.norm(np.cross(tp_v, floating_x)) # points superiorly
        # lat_d, med_d = self._calc_varus_angle(femur_tree)[1]
        # knee_gap_o = 0.5*(tp_lat+tp_med)
        # tp_plane = geoprimitives.Plane(knee_gap_o, knee_gap_v, tp_v, floating_x)

        # self._reset_tibia_kneegap_1()

        # current_knee_gap = self._calc_knee_gap_2(tp_plane, femur_points)
        # knee_shift = -(self.knee_gap - current_knee_gap) # in the negative y of the tibia acs
        # print('knee_shift {}'.format(knee_shift))
        # shift_t = knee_shift * knee_gap_v
        # self.models['tibiafibula'].gf.transformTranslate(shift_t)
        # self.models['tibiafibula'].update_landmarks()
        # self.models['tibiafibula'].update_acs()

    def _get_knee_floating_x(self):
        return math.norm(
            np.cross(
                self.models['tibiafibula'].acs.y,
                self.models['femur'].acs.z
            )
        )

    def _get_hip_floating_x(self):
        return math.norm(
            np.cross(
                self.models['femur'].acs.y,
                self.models['pelvis'].acs.z,
            )
        )

    def _get_knee_cs(self):
        """Returns knee joint coordinate system:
        [origin, flexion axis, rotation axis, abduction axis]
        """
        o = self.models['femur'].acs.o  # origin
        abd = self._get_knee_floating_x()  # abduction
        rot = self.models['tibiafibula'].acs.y  # rotation
        flex = self.models['femur'].acs.z  # flexion
        return np.array([o, flex, rot, abd])

    def _get_hip_cs(self):
        """Returns hip joint coordinate system:
        [origin, flexion axis, rotation axis, abduction axis]
        """
        o = self.models['pelvis'].landmarks[self.HJC]  # origin
        abd = self._get_hip_floating_x()  # abduction
        rot = self.models['femur'].acs.y  # rotation
        flex = self.models['pelvis'].acs.z  # flexion
        return np.array([o, flex, rot, abd])

    def _calc_varus_angle(self, femur_points_tree):

        # evaluate tibial plateau points

        tp_lat = self.models['tibiafibula'].gf.get_point_position(self.tp_lat_node)
        tp_med = self.models['tibiafibula'].gf.get_point_position(self.tp_med_node)
        # tibial plateau vector
        tp_v = tp_med - tp_lat

        # find closest femur condyle points to tp_lat and tp_med
        (med_d, lat_d), (fc_med_i, fc_lat_i) = femur_points_tree.query([tp_med, tp_lat])

        # femoral condyle vector
        fc_v = femur_points_tree.data[fc_med_i] - femur_points_tree.data[fc_lat_i]

        # project vectors to Z-Y plane
        Z = self.models['femur'].acs.z
        Y = self.models['tibiafibula'].acs.y
        tp_v_zy = np.array([np.dot(tp_v, Z), np.dot(tp_v, Y)])
        fc_v_zy = np.array([np.dot(fc_v, Z), np.dot(fc_v, Y)])

        # calc angle in Z-Y plane
        tp_fc_angle = abs(math.angle(tp_v_zy, fc_v_zy))

        # if lat condyle is higher than med condyle, negative rotation is needed
        if lat_d > med_d:
            tp_fc_angle *= -1.0

        # anticlock wise is positive
        # if (2.0*tp_fc_angle)>np.pi:
        #     tp_fc_angle = tp_fc_angle - np.pi

        # print('lat_d {}, med_d {}'.format(lat_d, med_d))
        # print('varus angle {}'.format(tp_fc_angle))

        return tp_fc_angle, (lat_d, med_d)

    def update_pelvis(self, pelvis_rigid):
        """Update position and orientation of the pelvis.
        Inputs:
        pelvis_rigid [list]: list of 6 floats describing a rigid
                             transformation in the global coordinate system 
                             of the pelvis - [tx, ty, tz, rx, ry, rz].
                             Rotation is about the origin of the pelvis
                             anatomic coordinate system.
        """
        # rigid transform pelvis
        self.models['pelvis'].gf.transformRigidRotateAboutP(
            pelvis_rigid,
            self.pelvis_origin,
        )
        self.models['pelvis'].update_landmarks()
        self.models['pelvis'].update_acs()

        # self.models['pelvis'].transformRigidAboutPoint(*pelvis_rigid, p=self.pelvis_origin)

    def update_femur(self, hip_rot):
        """Update position and orientation of the tibiafibula segment.
        inputs:
        hip_rot [list]: list of 3 floats describing hip flexion, rotation,
                        and adduction in radians.
        """

        # align femur to neutral hip CS
        self._reset_femur_hip()

        # get joint cs
        o, flex, rot, abd = self._get_hip_cs()

        # apply hip rotations:
        # flexion (pelvis-z)
        HJC = self.models['pelvis'].landmarks[self.HJC]
        self.models['femur'].gf.transformRotateAboutAxis(
            self.hip_flex_coeff * hip_rot[0],
            o, o + flex,
        )
        self.models['femur'].update_landmarks()
        self.models['femur'].update_acs()

        # abduction (floating)
        self.models['femur'].gf.transformRotateAboutAxis(
            self.hip_abd_coeff * hip_rot[2],
            o, o + abd
        )
        self.models['femur'].update_landmarks()
        self.models['femur'].update_acs()

        # rotations (femur-z)
        self.models['femur'].gf.transformRotateAboutAxis(
            self.hip_rot_coeff * hip_rot[1],
            o, o + rot,
        )
        self.models['femur'].update_landmarks()
        self.models['femur'].update_acs()

    def update_tibiafibula(self, knee_rot):
        """Update position and orientation of the tibiafibula segment.
        inputs:
        knee_rot [list of floats]: knee flexion and adduction angle in radians.
        """

        # align tibia to neutral knee CS
        self._reset_tibia_knee()

        # get joint cs
        o, flex, rot, abd = self._get_knee_cs()

        # apply knee rotations: flexion(femur-z)
        self.models['tibiafibula'].gf.transformRotateAboutAxis(
            self.knee_flex_coeff * knee_rot[0],
            o,
            o + flex,
        )
        self.models['tibiafibula'].update_landmarks()
        self.models['tibiafibula'].update_acs()

        # apply knee rotation: adduction(floating x)
        if self._allow_knee_adduction_dof:
            x = self._get_knee_floating_x()
            self.models['tibiafibula'].gf.transformRotateAboutAxis(
                self.knee_abd_coeff * knee_rot[1],
                o,
                o + abd
            )
            self.models['tibiafibula'].update_landmarks()
            self.models['tibiafibula'].update_acs()

        # reset knee gap
        self._reset_tibia_kneegap()

    def update_patella(self):
        """Update patella position and orientation. Patella is placed a fixed distance
        from the tibia in the tibial Y direction"""

        # align patella acs to tibial acs, centred at midpoint of tibial epicondyles
        ot = 0.5 * (self.models['tibiafibula'].landmarks['tibiafibula-LC'] +
                    self.models['tibiafibula'].landmarks['tibiafibula-MC'])
        cs_targ = np.array([ot,
                            ot + self.models['tibiafibula'].acs.x,
                            ot + self.models['tibiafibula'].acs.y,
                            ot + self.models['tibiafibula'].acs.z,
                            ])

        cs_source = self.models['patella'].acs.unit_array

        T1 = transform3D.directAffine(cs_source, cs_targ)
        self.models['patella'].gf.transformAffine(T1)
        self.models['patella'].update_landmarks()
        self.models['patella'].update_acs()

        # apply patella shift
        T2 = self.models['tibiafibula'].acs.x * self.patella_shift[0] + \
             self.models['tibiafibula'].acs.y * self.patella_shift[1] + \
             self.models['tibiafibula'].acs.z * self.patella_shift[2]

        self.models['patella'].gf.transformTranslate(T2)
        self.models['patella'].update_landmarks()
        self.models['patella'].update_acs()

    def update_all_models(self, pc_weights, pc_modes, pelvis_rigid, hip_rot, knee_rot):
        """Update the lower limb geometry by pc weights and rigid transformations

        Inputs:
        pc_weights [list of floats]: list of pc weights
        pc_modes [list of ints]: list of the pcs that the weights are for
        pelvis_rigid [1-d array]: an array of six elements defining the rigid body
            translation and rotation for the pelvis.
        hip_rot [1-d array]: an array of 3 radian angles for femur rotation about the
            hip joint (flexion, rotation, adduction)
        knee_rot [1-d array]: an array of radian angles for tibia-fibula rotation
            about the knee joint (flexion)
        """

        # evaluate shape model
        self.update_models_by_pcweights_sd(pc_weights, pc_modes)

        # rigid transform pelvis
        self.update_pelvis(pelvis_rigid)

        # place femur by hip rotation
        self.update_femur(hip_rot)

        # place tibia and fibula by knee_rot and default_knee_offset
        self.update_tibiafibula(knee_rot)

        # place patella relative to tibia
        self.update_patella()

    def update_all_models_multi_scaling(self, scalings, pelvis_rigid, hip_rot, knee_rot):
        """Update the lower limb geometry by isotropic scaling for each bone and
        rigid transformations.

        Inputs:
        scalings [2 lists]: a list of model names and a list of their scaling
        pelvis_rigid [1-d array]: an array of six elements defining the rigid body
            translation and rotation for the pelvis.
        hip_rot [1-d array]: an array of 3 radian angles for femur rotation about the
            hip joint (flexion, rotation, adduction)
        knee_rot [1-d array]: an array of radian angles for tibia-fibula rotation
            about the knee joint (flexion)
        """

        # scale models
        for mn, s in zip(scalings[0], scalings[1]):
            self.update_model_by_rigid_scale(mn, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, s)

        # rigid transform pelvis
        self.update_pelvis(pelvis_rigid)

        # place femur by hip rotation
        self.update_femur(hip_rot)

        # place tibia and fibula by knee_rot and default_knee_offset
        self.update_tibiafibula(knee_rot)

        # place patella relative to tibia
        self.update_patella()

    def update_all_models_uniform_scaling(self, scaling, pelvis_rigid, hip_rot, knee_rot):
        """Update the lower limb geometry by isotropic scaling for each bone and
        rigid transformations.

        Inputs:
        scaling [float]: scale factor for whole lower limb model
        pelvis_rigid [1-d array]: an array of six elements defining the rigid body
            translation and rotation for the pelvis.
        hip_rot [1-d array]: an array of 3 radian angles for femur rotation about the
            hip joint (flexion, rotation, adduction)
        knee_rot [1-d array]: an array of radian angles for tibia-fibula rotation
            about the knee joint (flexion)
        """

        # scale models
        self.update_models_by_uniform_rigid_scale(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, scaling)

        # rigid transform pelvis
        self.update_pelvis(pelvis_rigid)
        # self.models['pelvis'].transformRigidScaleAboutPoint(
        #                         *pelvis_rigid,
        #                         s=scaling,
        #                         p=self.pelvis_origin
        #                         )

        # place femur by hip rotation
        self.update_femur(hip_rot)

        # place tibia and fibula by knee_rot and default_knee_offset
        self.update_tibiafibula(knee_rot)

        # place patella relative to tibia
        self.update_patella()

    # default settings
    _reset_tibia_kneegap = _reset_tibia_kneegap_1


class LowerLimbRightAtlas(LowerLimbLeftAtlas):
    combined_pcs_filename = '../models/lower_limb/shapemodels/LLP26_right_mirrored_from_left_rigid.pc'
    bone_files = {'pelvis': ('../models/pelvis/template_meshes/pelvis_combined_cubic_mean_rigid_LLP26.geof',
                             '../models/pelvis/template_meshes/pelvis_combined_cubic_flat.ens',
                             '../models/pelvis/template_meshes/pelvis_combined_cubic_flat.mesh',
                             ),
                  'femur': ('../models/femur/template_meshes/femur_right_mirrored_from_left_mean_rigid_LLP26.geof',
                            '../models/femur/template_meshes/femur_right_quartic_flat.ens',
                            '../models/femur/template_meshes/femur_right_quartic_flat.mesh',
                            ),
                  'patella': (
                  '../models/patella/template_meshes/patella_right_mirrored_from_left_mean_rigid_LLP26.geof',
                  '../models/patella/template_meshes/patella_11_right.ens',
                  '../models/patella/template_meshes/patella_11_right.mesh',
                  ),
                  'tibiafibula': (
                  '../models/tibiafibula/template_meshes/tibia_fibula_cubic_right_mirrored_from_left_mean_rigid_LLP26.geof',
                  '../models/tibiafibula/template_meshes/tibia_fibula_right_cubic_flat.ens',
                  '../models/tibiafibula/template_meshes/tibia_fibula_right_cubic_flat.mesh',
                  ),
                  }
    HJC = 'pelvis-RHJC'
    side = 'right'

    hip_flex_coeff = 1.0
    hip_abd_coeff = 1.0
    hip_rot_coeff = 1.0
    knee_flex_coeff = 1.0
    knee_abd_coeff = 1.0
    knee_rot_coeff = -1.0

    patella_shift = np.array([50.0, 50.0, 10.0])
