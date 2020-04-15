"""
FILE: lowerlimbatlas.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: Classes and functions for the lowerlimb atlas model

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import numpy as np

from gias2.musculoskeletal.bonemodels.bonemodels import LowerLimbLeftAtlas, LowerLimbRightAtlas


def _trim_angle(a):
    if a < -np.pi:
        return a + 2 * np.pi
    elif a > np.pi:
        return a - 2 * np.pi
    else:
        return a


class LowerLimbAtlas(object):
    """
    Both sides
    """
    SHAPEMODESMAX = 100
    _allow_knee_adduction_dof = False
    _allow_knee_adduction_correction = False
    _neutral_params = [
        [0, ],
        [0, ],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, ],
        [0, ],
    ]
    N_PARAMS_PELVIS = 6
    N_PARAMS_HIP_L = 3
    N_PARAMS_HIP_R = 3
    N_PARAMS_KNEE_L = 1
    N_PARAMS_KNEE_R = 1
    N_PARAMS_RIGID = N_PARAMS_PELVIS + \
                     N_PARAMS_HIP_L + \
                     N_PARAMS_HIP_R + \
                     N_PARAMS_KNEE_L + \
                     N_PARAMS_KNEE_R

    def __init__(self, name):
        self.name = name
        self.models = {}
        self.model_elem_map = {}
        self.combined_model_gf = None
        self.combined_pcs = None
        self._map_model_params = None
        self._combined_param_map = None

        self.ll_l = LowerLimbLeftAtlas(self.name + '_l')
        self.ll_r = LowerLimbRightAtlas(self.name + '_r')

        self._pelvis_rigid = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self._hip_rot_l = np.array([0.0, 0.0, 0.0])
        self._hip_rot_r = np.array([0.0, 0.0, 0.0])
        self._knee_rot_l = np.array([0.0, 0.0, 0.0])
        self._knee_rot_r = np.array([0.0, 0.0, 0.0])
        self.shape_modes = [0, ]
        self._shape_mode_weights = np.zeros(self.SHAPEMODESMAX, dtype=float)
        self.uniform_scaling = 1.0
        self.pelvis_scaling = 1.0
        self.femur_scaling_l = 1.0
        self.femur_scaling_r = 1.0
        self.petalla_scaling_l = 1.0
        self.petalla_scaling_r = 1.0
        self.tibfib_scaling_l = 1.0
        self.tibfib_scaling_r = 1.0
        self.last_transform_set = None

        self._shape_model_x = None
        self._uniform_scaling_x = None
        self._per_bone_scaling_x = None

    @property
    def pelvis_rigid(self):
        return self._pelvis_rigid

    @pelvis_rigid.setter
    def pelvis_rigid(self, value):
        if len(value) != 6:
            raise ValueError('input vector not of length 6')
        else:
            self._pelvis_rigid = np.array([value[0], value[1], value[2],
                                           _trim_angle(value[3]),
                                           _trim_angle(value[4]),
                                           _trim_angle(value[5]),
                                           ])
            self.update_pelvis(self._pelvis_rigid)

    @property
    def hip_rot_l(self):
        return self._hip_rot_l

    @property
    def hip_rot_r(self):
        return self._hip_rot_r

    @hip_rot_l.setter
    def hip_rot_l(self, value):
        if len(value) != 3:
            raise ValueError('input vector not of length 3')
        else:
            self._hip_rot_l = np.array([_trim_angle(v) for v in value])
            self.ll_l.update_femur(self._hip_rot_l)

    @hip_rot_r.setter
    def hip_rot_r(self, value):
        if len(value) != 3:
            raise ValueError('input vector not of length 3')
        else:
            self._hip_rot_r = np.array([_trim_angle(v) for v in value])
            self.ll_r.update_femur(self._hip_rot_r)

    @property
    def knee_rot_l(self):
        if self._allow_knee_adduction_dof:
            return self._knee_rot_l[[0, 2]]
        else:
            return self._knee_rot_l[[0]]

    @property
    def knee_rot_r(self):
        if self._allow_knee_adduction_dof:
            return self._knee_rot_r[[0, 2]]
        else:
            return self._knee_rot_r[[0]]

    @knee_rot_l.setter
    def knee_rot_l(self, value):
        if self._allow_knee_adduction_dof:
            self._knee_rot_l[0] = _trim_angle(value[0])
            self._knee_rot_l[2] = _trim_angle(value[1])
        else:
            self._knee_rot_l[0] = _trim_angle(value[0])

        self.ll_l.update_tibiafibula(self._knee_rot_l[[0, 2]])
        self.ll_l.update_patella()

    @knee_rot_r.setter
    def knee_rot_r(self, value):
        if self._allow_knee_adduction_dof:
            self._knee_rot_r[0] = _trim_angle(value[0])
            self._knee_rot_r[2] = _trim_angle(value[1])
        else:
            self._knee_rot_r[0] = _trim_angle(value[0])

        self.ll_r.update_tibiafibula(self._knee_rot_r[[0, 2]])
        self.ll_r.update_patella()

    @property
    def shape_mode_weights(self):
        return self._shape_mode_weights[self.shape_modes]

    @shape_mode_weights.setter
    def shape_mode_weights(self, value):
        if len(value) != len(self.shape_modes):
            raise ValueError('Length of value does not match length of self.shape_modes ({})'.format(self.shape_modes))

        self._shape_mode_weights[self.shape_modes] = value
        self._update_models_by_pcweights_sd(value, self.shape_modes)

    # gets a flat array, sets using a list of arrays.
    @property
    def shape_model_x(self):
        self._shape_model_x = np.hstack([
            self.shape_mode_weights,
            self.pelvis_rigid,
            self.hip_rot_l,
            self.hip_rot_r,
            self.knee_rot_l,
            self.knee_rot_r,
        ])
        return self._shape_model_x

    @shape_model_x.setter
    def shape_model_x(self, value):
        # a = len(self.shape_modes)
        self._shape_model_x = value
        self.shape_mode_weights = value[0]
        self.pelvis_rigid = value[1]
        self.hip_rot_l = value[2]
        self.hip_rot_r = value[3]
        self.knee_rot_l = value[4]
        self.knee_rot_r = value[5]
        self.lastT_transform_set = self.shape_model_x
        self.update_all_models(
            self.shape_mode_weights,
            self.shape_modes,
            self.pelvis_rigid,
            self.hip_rot_l,
            self.hip_rot_r,
            self.knee_rot_l,
            self.knee_rot_r,
        )

    def enable_knee_adduction_correction(self):
        self.ll_l.enable_knee_adduction_correction()
        self.ll_r.enable_knee_adduction_correction()
        self._allow_knee_adduction_correction = True

    def disable_knee_adduction_correction(self):
        self.ll_l.disable_knee_adduction_correction()
        self.ll_r.disable_knee_adduction_correction()
        self._allow_knee_adduction_correction = False

    def enable_knee_adduction_dof(self):
        self.ll_l.enable_knee_adduction_dof()
        self.ll_r.enable_knee_adduction_dof()
        self._allow_knee_adduction_dof = True
        self._neutral_params = [[0, ], [0, ], [0, 0, 0, 0, 0, 0], [0, 0, 0], [0, 0, 0, ], [0, 0, ], [0, 0, ]]
        self.N_PARAMS_KNEE_L = 2
        self.N_PARAMS_KNEE_R = 2
        self.N_PARAMS_RIGID = self.N_PARAMS_PELVIS + \
                              self.N_PARAMS_HIP_L + self.N_PARAMS_HIP_R + \
                              self.N_PARAMS_KNEE_L + self.N_PARAMS_KNEE_R

    def disable_knee_adduction_dof(self):
        self.ll_l.disable_knee_adduction_dof()
        self.ll_r.disable_knee_adduction_dof()
        self._allow_knee_adduction_dof = False
        self._neutral_params = [[0, ], [0, ], [0, 0, 0, 0, 0, 0], [0, 0, 0], [0, 0, 0, ], [0, ], [0, ]]
        self.N_PARAMS_KNEE_L = 1
        self.N_PARAMS_KNEE_R = 1
        self.N_PARAMS_RIGID = self.N_PARAMS_PELVIS + \
                              self.N_PARAMS_HIP_L + self.N_PARAMS_HIP_R + \
                              self.N_PARAMS_KNEE_L + self.N_PARAMS_KNEE_R

    # @property
    # def uniformScalingX(self):
    #     self._uniformScalingX = np.hstack([
    #                             self.uniformScaling,
    #                             self.pelvisRigid,
    #                             self.hipRot,
    #                             self.kneeRot
    #                             ])
    #     return self._uniformScalingX

    # @uniformScalingX.setter
    # def uniformScalingX(self, value):
    #     print value
    #     a = 1
    #     self._uniformScalingX = value
    #     self.uniformScaling = value[0]
    #     self.pelvisRigid = value[1]
    #     self.hipRot = value[2]
    #     self.kneeRot = value[3]

    #     # propagate isotropic scaling to each bone
    #     self.pelvisScaling = self.uniformScaling
    #     self.femurScaling = self.uniformScaling
    #     self.patellaScaling = self.uniformScaling
    #     self.tibfibScaling = self.uniformScaling

    #     self.lastTransformSet = self.uniformScalingX

    # @property
    # def perBoneScalingX(self):
    #     self._perBoneScalingX = np.hstack([
    #                             self.pelvisScaling,
    #                             self.femurScaling,
    #                             self.patellaScaling,
    #                             self.tibfibScaling,
    #                             self.pelvisRigid,
    #                             self.hipRot,
    #                             self.kneeRot
    #                             ])
    #     return self._perBoneScalingX

    # @perBoneScalingX.setter
    # def perBoneScalingX(self, value):
    #     a = 4
    #     self._perBoneScalingX = value
    #     self.pelvisScaling = value[0][1][0]
    #     self.femurScaling = value[0][1][1]
    #     self.patellaScaling = value[0][1][2]
    #     self.tibfibScaling = value[0][1][3]
    #     self.pelvisRigid = value[1]
    #     self.hipRot = value[2]
    #     self.kneeRot = value[3]
    #     self.lastTransformSet = self.perBoneScalingX

    def _update_model_dict(self):
        for model_name, model in self.ll_l.models.items():
            self.models[model_name + '-l'] = model

        for model_name, model in self.ll_r.models.items():
            self.models[model_name + '-r'] = model

        # use the left pelvis as the reference
        self.models['pelvis'] = self.ll_l.models['pelvis']

    # def load_models_left(self, *args, **kwargs):
    #     self.ll_l.load_models(*args, **kwargs)
    #     self._update_model_dict()

    # def load_models_right(self, *args, **kwargs):
    #     self.ll_r.load_models(*args, **kwargs)
    #     self._update_model_dict()

    def load_bones(self):
        self.ll_l.load_bones()
        self.ll_r.load_bones()
        self._update_model_dict()

    # def load_combined_pcs_left(self, filename):
    #     """ 
    #     Load the combined left lower limb pca model
    #     """
    #     self.ll_l.load_combined_pcs(filename)

    # def load_combined_pcs_right(self, filename):
    #     """ 
    #     Load the combined right lower limb pca model
    #     """
    #     self.ll_r.load_combined_pcs(filename)

    def _update_models_by_pcweights_sd(self, pc_weights, pc_modes):
        self.ll_l.update_models_by_pcweights_sd(pc_weights, pc_modes)
        self.ll_r.update_models_by_pcweights_sd(pc_weights, pc_modes)
        # self._update_model_dict()

    # def update_models_by_combined_params_left(self, p):
    #     self.ll_l.update_models_by_combined_params(p)

    # def update_models_by_combined_params_right(self, p):
    #     self.ll_r.update_models_by_combined_params(p)

    # def update_models_by_uniform_rigid_scale(self, *args):
    #     """ Rotation and scaling is about model CoM
    #     """
    #     self.ll_l.update_models_by_uniform_rigid_scale(*args)
    #     self.ll_r.update_models_by_uniform_rigid_scale(*args)

    # def update_model_by_rigid_scale(self, modelname, tx, ty, tz, rx, ry, rz, s):
    #     """ Rotation and scaling is about model CoM
    #     """
    #     _modelname, side = modelname.split('-')
    #     if side=='left':
    #         self.ll_l.update_model_by_rigid_scale(_modelname, tx, ty, tz, rx, ry, rz, s)
    #     elif side=='right':
    #         self.ll_r.update_model_by_rigid_scale(_modelname, tx, ty, tz, rx, ry, rz, s)

    def update_pelvis(self, pelvis_rigid):
        """
        Update position and orientation of the pelvis.
        Inputs:
        -------
        pelvis_rigid [list]: list of 6 floats describing a rigid
                             transformation in the global coordinate system 
                             of the pelvis - [tx, ty, tz, rx, ry, rz].
                             Rotation is about the origin of the pelvis
                             anatomic coordinate system.
        """
        self.ll_l.update_pelvis(pelvis_rigid)
        self.ll_r.update_pelvis(pelvis_rigid)
        # self._update_model_dict()

    def update_all_models(
            self, pc_weights, pc_modes, pelvis_rigid, hip_rot_l, hip_rot_r,
            knee_rot_l, knee_rot_r
    ):
        """
        Update the lower limb geometry by pc weights and rigid transformations

        Inputs:
        -------
        pc_weights [list of floats]: list of pc weights
        pc_modes [list of ints]: list of the pcs that the weights are for
        pelvis_rigid [1-d array]: an array of six elements defining the rigid
            body translation and rotation for the pelvis.
        hip_rot_l, hip_rot_r [1-d array]: an array of 3 radian angles for left
            and right femur rotation about the hip joint (flexion, rotation,
            adduction)
        knee_rot_l, knee_rot_r [1-d array]: an array of radian angles for
            tibia-fibula rotation about the knee joint (flexion)
        """

        # evaluate shape model
        self.shape_modes = pc_modes
        self.shape_mode_weights = pc_weights

        # rigid transform pelvis
        self.pelvis_rigid = pelvis_rigid

        # place femur by hip rotation
        self.hip_rot_l = hip_rot_l
        self.hip_rot_r = hip_rot_r

        # place tibia and fibula by knee_rot and default_knee_offset
        self.knee_rot_l = knee_rot_l
        self.knee_rot_r = knee_rot_r

        # self._update_model_dict()
