"""
FILE: modelcore.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: Core functions and classes for bone models

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""
import logging
import xml.etree.cElementTree as ET

import numpy as np
from numpy.linalg import inv

from gias2.common import transform3D
from gias2.fieldwork.field import geometric_field
from gias2.learning import PCA
from gias2.musculoskeletal import fw_model_landmarks as model_landmarks

log = logging.getLogger(__name__)


# ===========================================#
# accessory functions                       #
# ===========================================#
class ACSCartesian(object):
    """ Cartesian anatomic coordinate system
    """

    global_cs = np.array([[0, 0, 0],
                          [1, 0, 0],
                          [0, 1, 0],
                          [0, 0, 1]])

    def __init__(self, o, x, y, z):
        """Define a 3D cartesian coordinate system by an origin point o,
        and unit vectors in each of the coordinates x, y, z.
        """
        self.o = None
        self.x = None
        self.y = None
        self.z = None
        self.unit_array = None
        self.local_transform = None  # affine matrix that transfrom from global to local frame
        self.inv_local_transform = None  # affine matrix that transfrom from local to global frame
        self.update(o, x, y, z)

    def update(self, o, x, y, z):
        """Define a 3D cartesian coordinate system by an origin point o,
        and unit vectors in each of the coordinates x, y, z.
        """
        self.o = np.array(o)
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)
        self.unit_array = np.array([self.o,
                                    self.o + self.x,
                                    self.o + self.y,
                                    self.o + self.z
                                    ])
        # self.local_transform = transform3D.calcAffineMatrixSVD(
        # self.unit_array, self.global_cs
        # )
        self.local_transform = transform3D.directAffine(
            self.unit_array,
            self.global_cs
        )
        self.inv_local_transform = inv(self.local_transform)

    def get_unit_array(self):
        return self.unit_array

    def make_unit_array_transform_rigid(self, tx, ty, tz, rx, ry, rz, ro):
        return transform3D.transformRigid3DAboutP(
            self.unit_array,
            [tx, ty, tz, rx, ry, rz],
            ro
        )

    def map_local(self, x):
        """Calculate the local coordinates of points in x with global
        coordinates.
        """
        return transform3D.transformAffine(x, self.local_transform)

    def rotate_local(self, x):
        """
        rotate a 3-vector in global space to the local space
        """
        return np.dot(self.local_transform[:3, :3], x)

    def map_global(self, x):
        """Calculate the global coordinates of points in x with local
        coordinates.
        """
        return transform3D.transformAffine(x, self.inv_local_transform)

    def rotate_global(self, x):
        """
        rotate a 3-vector in the local space to the global space
        """
        return np.dot(self.inv_local_transform[:3, :3], x)


def make_source_landmark_getter(landmark_names, side=None):
    """ Creates a function to return the coordinates of landmarks
    in a multibone atlas.

    inputs
    ------
    landmark_names : list of str
        List of landmark names
    side : str [optional]
        'l' or 'r' if the names of models will be sided
    """
    landmark_labels = []
    for ln in landmark_names:
        body_name, landmark_name = ln.split('-')
        if side:
            body_name += '_{}'.format(side)
        landmark_labels.append((body_name, ln))

    def get_source_landmarks(model, coords):
        for i, (bn, ln) in enumerate(landmark_labels):
            coords[i] = model.models[bn].landmarks[ln]

        return coords

    return get_source_landmarks


# ==================================================#
# region on bone surface                           #
# ==================================================#
class DiscretisedRegion(object):
    def __init__(self, vertices, faces):
        self.vertices = vertices
        self.faces = faces


class AttachmentRegion(DiscretisedRegion):

    def __init__(self, name=None, v=None, f=None, matpoints=None, end=None,
                 number=None):
        self.vertices = v
        if self.vertices is None:
            self.vertices = []
        self.faces = f
        if self.faces is None:
            self.faces = []
        self.matpoints = matpoints
        if self.matpoints is None:
            self.matpoints = []
        self.name = name
        self.end = end
        self.number = number

    def _to_xml(self, tree):
        attachment_el = ET.SubElement(tree, 'attachment')
        attachment_el.set('attachment_number', str(self.number))
        attachment_el.set('name', self.name)
        attachment_el.set('end', self.end)
        verts_el = ET.SubElement(attachment_el, 'vertices')
        verts_el.text = ', '.join('{}'.format(vi) for vi in self.vertices)
        faces_el = ET.SubElement(attachment_el, 'faces')
        faces_el.text = ', '.join('{}'.format(fi) for fi in self.faces)
        matpoints_el = ET.SubElement(attachment_el, 'mat_points')
        self._matpoints_to_xml(matpoints_el)

    def _matpoints_to_xml(self, matpoints_el):
        for mi, mp in enumerate(self.matpoints):
            matpoint_el = ET.SubElement(matpoints_el, 'mat_point')
            matpoint_el.set('id', str(mi))
            elem_el = ET.SubElement(matpoint_el, 'element')
            elem_el.text = str(mp[0])
            xi_el = ET.SubElement(matpoint_el, 'xi')
            xi_el.text = '{:12.9f}, {:12.9f}'.format(*mp[1])

    def _from_xml(self, region):
        self.name = region.attrib['name']
        self.end = region.attrib['end']
        self.number = region.attrib['attachment_number']
        properties = dict([(p.tag, p) for p in region])
        try:
            self.vertices = np.array([int(x) for x in properties['vertices'].text.split(',')])
        except KeyError:
            pass
        try:
            self.faces = np.array([int(x) for x in properties['faces'].text.split(',')])
        except KeyError:
            pass
        try:
            self.matpoints = self._matpoints_from_xml(properties['mat_points'])
        except KeyError:
            pass

    def _matpoints_from_xml(self, matpoints_el):
        mps = []
        ids = []
        for mp_el in matpoints_el:
            ids.append(int(mp_el.attrib['id']))
            properties = dict([(p.tag, p) for p in mp_el])
            elem = int(properties['element'].text)
            xi = [float(x) for x in properties['xi'].text.split(',')]
            mps.append((elem, xi))
        return mps


class BoneAttachmentRegions(object):

    def __init__(self, bone=None, disc=None, geof=None, stl=None,
                 attachmentsource=None, regions=None):
        self.bone = bone
        self.atlas_disc = disc
        self.atlas_geof = geof
        self.atlas_stl = stl
        self.attachment_source = attachmentsource
        self.regions = regions

        if self.regions is None:
            self.regions = []

    def to_xml(self, filename):
        attachment_output = ET.Element('attachment_sites')
        attachment_output.set('bone', self.bone)
        attachment_output.set('source', self.attachment_source)
        attachment_output.set('atlas_geof', self.atlas_geof)
        attachment_output.set('atlas_stl', self.atlas_stl)
        attachment_output.set('atlas_disc', self.atlas_disc)
        for r in self.regions:
            r._to_xml(attachment_output)

        tree = ET.ElementTree(attachment_output)
        tree.write(filename)

    def from_xml(self, filename):

        tree = ET.parse(filename)
        root = tree.getroot()

        self.bone = root.attrib['bone']
        self.atlas_disc = root.attrib['atlas_disc']
        self.atlas_geof = root.attrib['atlas_geof']
        self.atlas_stl = root.attrib['atlas_stl']
        self.attachment_source = root.attrib['source']

        self.regions = []
        for l in root:
            if l.tag == 'attachment':
                s = AttachmentRegion()
                s._from_xml(l)
                self.regions.append(s)
            else:
                raise IOError('unrecognised attachment file object {}'.format(l.tag))


# ==================================================#
# Individual bone models                           #
# ==================================================#
class BoneModel(object):

    def __init__(self, name, gf):
        self.name = name
        self.gf = gf
        self.acs = ACSCartesian((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
        self.landmarks = {}
        self._landmark_evaluators = {}
        self._source_field_parameters = self.gf.field_parameters.copy()

    def init_landmarks(self, landmark_names, side=None):
        self.landmarks = {}
        for ln in landmark_names:
            if side is None:
                self._landmark_evaluators[ln] = model_landmarks.makeLandmarkEvaluator(
                    ln, self.gf
                )
            else:
                self._landmark_evaluators[ln] = model_landmarks.makeLandmarkEvaluator(
                    ln, self.gf, side=side
                )
        self.update_landmarks()
        self.update_acs()

    def update_gf(self, gfparams):
        # self.gf.set_field_parameters(gfparams)
        self.gf.field_parameters = gfparams
        self.update_landmarks()
        self.update_acs()

    def update_gf_params(self, gfparams):
        """ Does not actually update the gf
        """
        self._gf_params = gfparams
        self.update_landmarks()
        self.update_acs()

    def update_landmarks(self):
        for ln, leval in list(self._landmark_evaluators.items()):
            self.landmarks[ln] = leval(self.gf.field_parameters)

        return self.landmarks

    def update_acs(self):
        """ Define this in inherited class
        """
        raise NotImplementedError('This method has not been defined')

    def transformRigid(self, tx, ty, tz, rx, ry, rz):
        """Transform the original model by a rigid
        transformation
        """
        T = np.array([tx, ty, tz, rx, ry, rz])
        xT = transform3D.transformRigid3D(
            self._source_field_parameters[:, :, 0].T,
            T)
        self.gf.field_parameters = xT.T[:, :, np.newaxis]
        self.update_landmarks()
        self.update_acs()

    def transformRigidAboutPoint(self, tx, ty, tz, rx, ry, rz, p):
        """Transform the original model by a rigid
        transformation
        """
        T = np.array([tx, ty, tz, rx, ry, rz])
        xT = transform3D.transformRigid3DAboutP(
            self._source_field_parameters[:, :, 0].T,
            T, p)
        self.gf.field_parameters = xT.T[:, :, np.newaxis]
        self.update_landmarks()
        self.update_acs()

    def transformRigidScale(self, tx, ty, tz, rx, ry, rz, s):
        """Transform the original model by a rigid
        transformation and an isotropic scaling
        """
        T = np.array([tx, ty, tz, rx, ry, rz, s])
        xT = transform3D.transformRigidScale3D(
            self._source_field_parameters[:, :, 0].T,
            T)
        self.gf.field_parameters = xT.T[:, :, np.newaxis]
        self.update_landmarks()
        self.update_acs()

    def transformRigidScaleAboutPoint(self, tx, ty, tz, rx, ry, rz, s, p):
        """Transform the original model by a rigid
        transformation and an isotropic scaling
        """
        T = np.array([tx, ty, tz, rx, ry, rz, s])
        xT = transform3D.transformRigidScale3DAboutP(
            self._source_field_parameters[:, :, 0].T,
            T, p)
        self.gf.field_parameters = xT.T[:, :, np.newaxis]
        self.update_landmarks()
        self.update_acs()


# ==================================================#
# Multiple bone models                             #
# ==================================================#
class MultiBoneAtlas(object):
    combined_model_field_basis = {'tri10': 'simplex_L3_L3',
                                  'tri15': 'simplex_L4_L4',
                                  'quad44': 'quad_L3_L3',
                                  'quad55': 'quad_L4_L4',
                                  'quad54': 'quad_L4_L3',
                                  }

    def __init__(self, name):
        self.name = name
        self.models = {}
        self.model_elem_map = {}
        self.combined_model_gf = None
        self.combined_pcs = None
        self._map_model_params = None
        self._combined_param_map = None

    def load_models(self, model_names, model_classes, model_filenames,
                    combined_gf_name=None, combined_ens_name=None,
                    combined_mesh_name=None):
        """
        Load in the individual models of the lower limb
        and create the combined model and a map of model name to
        combined model element number.

        Inputs
        ------
        model_names [list of strings]: a list of the names of each model
        model_filenames [dictionary]: a dict of (modelname: [filenames]) 
        """
        if combined_gf_name is None:
            combined_gf_name = self.name
        if combined_ens_name is None:
            combined_ens_name = self.name
        if combined_mesh_name is None:
            combined_mesh_name = self.name

        self.combined_model_gf = geometric_field.geometric_field(
            combined_gf_name, 3, field_dimensions=2,
            field_basis=self.combined_model_field_basis
        )
        self.combined_model_gf.ensemble_field_function.name = combined_ens_name
        self.combined_model_gf.ensemble_field_function.mesh.name = combined_mesh_name

        self.models = {}
        self.model_elem_map = {}
        for mn in model_names:
            log.debug(('loading models: {}'.format(mn)))
            self.models[mn] = model_classes[mn](
                mn,
                geometric_field.load_geometric_field(
                    model_filenames[mn][0],
                    model_filenames[mn][1],
                    model_filenames[mn][2],
                )
            )
            elem_i = self.combined_model_gf.add_element_with_parameters(
                self.models[mn].gf.ensemble_field_function,
                self.models[mn].gf.get_field_parameters(),
                tol=0
            )

            self.model_elem_map[mn] = elem_i

        # get the global node numbers for each model
        elem2ens = self.combined_model_gf.ensemble_field_function.mapper._element_to_ensemble_map
        self._combined_param_map = {}
        for mn, combined_elem in list(self.model_elem_map.items()):
            self._combined_param_map[mn] = [elem2ens[combined_elem][x][0][0] for x in
                                            sorted(elem2ens[combined_elem].keys())]

    def load_combined_pcs(self, filename):
        """ Load the combined lower limb pca model
        """
        self.combined_pcs = PCA.loadPrincipalComponents(filename)

    def update_models_by_pcweights_sd(self, pc_weights, pc_modes):
        """Evaluate and set the parameters of each model from pc weights
        in terms of standard deviations.

        Inputs:
        pc_weights [1-d array]: array of floats of standard deviations
        pc_modes [1-d array]: array of integers corresponding to the modes

        """

        # evaluate combined model params
        combined_model_params = self.combined_pcs.reconstruct(
            self.combined_pcs.getWeightsBySD(pc_modes, pc_weights),
            pc_modes
        ).reshape((3, -1, 1))

        # separate params for each model
        for model_name, model in list(self.models.items()):
            model_params = combined_model_params[:, self._combined_param_map[model_name], :]
            model.update_gf(model_params)

    def update_models_by_combined_params(self, p):
        for model_name, model in list(self.models.items()):
            model_params = p[:, self._combined_param_map[model_name], :]
            model.update_gf(model_params)

    def update_models_by_uniform_rigid_scale(self, tx, ty, tz, rx, ry, rz, s):
        """ Rotation and scaling is about model CoM
        """
        for model in list(self.models.values()):
            model.transformRigidScale(tx, ty, tz, rx, ry, rz, s)

    def update_model_by_rigid_scale(self, modelname, tx, ty, tz, rx, ry, rz, s):
        """ Rotation and scaling is about model CoM
        """
        self.models[modelname].transformRigidScale(tx, ty, tz, rx, ry, rz, s)
