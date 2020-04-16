"""
FILE: smutils.py
LAST MODIFIED: 05-12-2016 
DESCRIPTION: Utility functions for SimpleMeshes

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""
import copy
import itertools
import logging

import numpy as np
from scipy.spatial.ckdtree import cKDTree
from scipy.stats import mode

from gias2.mesh import simplemesh

log = logging.getLogger(__name__)


def make_sub_mesh(sm, faceinds):
    """
    Create a mesh from face indices faceinds in the mesh sm
    """
    old_faces = sm.f[faceinds, :]
    unique_old_v_inds = np.unique(old_faces.ravel())
    unique_new_v_inds = np.arange(len(unique_old_v_inds), dtype=int)
    v_ind_map = dict(zip(unique_old_v_inds, unique_new_v_inds))

    # new_v = np.array([sm.v[i,:] for i in unique_old_v_inds])
    new_v = np.array(sm.v[unique_old_v_inds, :])
    new_f = np.zeros((len(faceinds), 3), dtype=int)

    for vi in unique_old_v_inds:
        new_f[old_faces == vi] = v_ind_map[vi]

    return simplemesh.SimpleMesh(new_v, new_f)


def set_1ring_faces(sm):
    """
    Create a dict of the adjacent faces of every face in sm
    """
    log.debug('setting 1-ring for faces')

    faces_1ring_faces = {}
    # share_edge_sets = [None, None, None]
    for fi, f in enumerate(sm.f):
        # find 3 adj faces that share a side with f
        shared_edge_set_0 = set(sm.faces1Ring[f[0]]).intersection(set(sm.faces1Ring[f[1]])).difference([fi])
        shared_edge_set_1 = set(sm.faces1Ring[f[0]]).intersection(set(sm.faces1Ring[f[2]])).difference([fi])
        shared_edge_set_2 = set(sm.faces1Ring[f[1]]).intersection(set(sm.faces1Ring[f[2]])).difference([fi])
        faces_1ring_faces[fi] = shared_edge_set_0.union(shared_edge_set_1).union(shared_edge_set_2)

    sm.faces1RingFaces = faces_1ring_faces


def partition_regions(sm, maxfaces):
    """
    Partition the mesh into regions of up to maxfaces connected faces.
    If maxfaces is inf, partitions the mesh into connected regions.

    returns
    -------
    label_faces: a dict of label number and the faces of that label
    face_labels: an array of the label number of each face in sm
    """

    remaining_faces = set(range(len(sm.f)))
    label_faces = {}
    face_labels = np.zeros(len(sm.f), dtype=int)
    reg_label = 0

    # while there are unpartitioned faces in sm
    while remaining_faces:
        reg_front = {min(remaining_faces)}  # region seed face
        reg_faces = [min(remaining_faces), ]
        remaining_faces.remove(min(remaining_faces))  # remove seed face from remaining
        reg_nfaces = 1

        # while current region is below max size
        while reg_nfaces < maxfaces:
            # for each face on the front
            try:
                front_f = reg_front.pop()
            except KeyError:
                # front is empty
                break

            # get remaining adjacent faces
            adj_f = [f for f in sm.faces1RingFaces[front_f] if f in remaining_faces]
            # add remaining adjacent faces to region 
            reg_faces += adj_f
            reg_nfaces += len(adj_f)

            # update remaining and front sets
            [remaining_faces.remove(f) for f in adj_f]
            reg_front = reg_front.union(adj_f)

        # record this regions faces and label number
        label_faces[reg_label] = reg_faces
        face_labels[reg_faces] = reg_label
        reg_label += 1

    return label_faces, face_labels


def make_region_meshes(sm, region_faces):
    """
    Given a mesh a list of face lists, create a mesh for each face list
    """
    meshes = []
    for reg_faces in region_faces.values():
        meshes.append(make_sub_mesh(sm, reg_faces))

    return meshes


def remove_small_regions(sm):
    """
    Return a mesh of the largest connected region in sm
    """

    # calculate adjacent faces for each face
    sm.set1RingFaces()

    # partition mesh by connected regions
    region_faces, face_labels = partition_regions(sm, np.inf)
    log.debug('found {} regions'.format(len(region_faces)))

    # get largest region
    largest_reg = None
    largest_reg_nfaces = 0
    for rn, rf in region_faces.items():
        if len(rf) > largest_reg_nfaces:
            largest_reg = rn
            largest_reg_nfaces = len(rf)

    log.debug('keeping largest region with {} faces'.format(largest_reg_nfaces))

    # create new mesh with just the largest region
    largest_region_mesh = make_sub_mesh(sm, region_faces[largest_reg])
    return largest_region_mesh


def remove_small_regions_2(sm, k):
    """
    Removes regions with less than k faces
    """

    # calculate adjacent faces for each face
    sm.set1RingFaces()

    # partition mesh by connected regions
    region_faces, face_labels = partition_regions(sm, np.inf)
    log.debug('found {} regions'.format(len(region_faces)))

    # find regions to keep
    keep_faces = []
    for rn, rf in region_faces.items():
        if len(rf) > k:
            keep_faces += rf

    # create mesh with kept regions
    new_mesh = make_sub_mesh(sm, keep_faces)
    log.debug('keeping {} faces'.format(len(keep_faces)))

    return new_mesh


def partition_mesh(sm, maxfaces, minfaces):
    """
    Partitions sm into regions with upper and lower faces bounds
    """

    region_faces, face_labels = partition_regions(sm, maxfaces)
    log.debug('merging {} regions'.format(len(region_faces)))
    merge_regions(sm, region_faces, face_labels, minfaces)
    log.debug('making {} region meshes'.format(len(region_faces)))
    region_sms = make_region_meshes(sm, region_faces)
    return region_sms


def merge_regions(sm, region_faces, face_labels, min_faces):
    """
    Given a mesh and a partitioning of its faces, merge regions with fewer
    than minfaces faces into the neighbouring region with the longest shared
    border

    inputs
    ------
    sm : Simplemesh
        The mesh with regions to merge
    region_faces : dict
        a dictionary mapping region number to the face indices of that region
    face_labels : np.Array
        an array of the region number of each face in sm
    minfaces : int
        the minimum number of faces a region can have. Any smaller regions
        will be merged
    """

    face_centre_tree = cKDTree(sm.faceBarycenters)

    def find_adj_by_distance(r_faces):
        rb_centres = sm.faceBarycenters[r_faces]
        _d, _i = face_centre_tree.query(rb_centres, k=min_faces * 2)
        ext_faces = set(itertools.chain.from_iterable(_i)).difference(r_faces)
        return ext_faces

    is_done = False
    while not is_done:
        # get regions that are too small
        regions_to_merge = [ri for ri in region_faces.keys() if len(region_faces[ri]) < min_faces]

        # get smallest region above face min
        regions_to_keep = [ri for ri in region_faces.keys() if ri not in regions_to_merge]
        regions_to_keep_sizes = np.array([len(region_faces[ri]) for ri in regions_to_keep])
        smallest_keep_reg_label = regions_to_keep[regions_to_keep_sizes.argmin()]

        if regions_to_merge:
            for ri in regions_to_merge:
                reg_faces = region_faces[ri]
                # get external adjacents faces for this region
                reg_faces_set = set(reg_faces)
                reg_adj_faces = [sm.faces1RingFaces[f] for f in reg_faces]
                reg_adj_faces_set = set(itertools.chain.from_iterable(reg_adj_faces))
                ext_adj_faces_set = reg_adj_faces_set.difference(reg_faces_set)
                if len(ext_adj_faces_set) == 0:
                    ext_adj_faces_set = find_adj_by_distance(reg_faces)

                if len(ext_adj_faces_set) == 0:
                    # merge into smallest non-merge region
                    parent_label = smallest_keep_reg_label
                else:
                    # get most common label of external adjacent faces
                    ext_adj_labels = [face_labels[i] for i in ext_adj_faces_set]
                    parent_label = mode(ext_adj_labels)[0][0]

                # merge this region into the parent label region
                region_faces[parent_label] += reg_faces
                del region_faces[ri]
                face_labels[face_labels == ri] = parent_label
        else:
            is_done = True


def merge_sms(sms):
    new_sm = copy.deepcopy(sms[0])
    for sm in sms[1:]:
        v_offset = new_sm.v.shape[0]
        new_sm.v = np.vstack([new_sm.v, sm.v])
        new_sm.f = np.vstack([new_sm.f, np.array(sm.f) + v_offset])

    return new_sm
