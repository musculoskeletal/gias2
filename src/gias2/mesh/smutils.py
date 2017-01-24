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

import numpy as np
from gias2.mesh import simplemesh

def make_sub_mesh(sm, faceinds):
    """
    Create a mesh from face indices faceinds in the mesh sm
    """
    old_faces = sm.f[faceinds,:]
    unique_old_v_inds = np.unique(old_faces.ravel())
    unique_new_v_inds = np.arange(len(unique_old_v_inds), dtype=int)
    v_ind_map = dict(zip(unique_old_v_inds, unique_new_v_inds))

    # new_v = np.array([sm.v[i,:] for i in unique_old_v_inds])
    new_v = np.array(sm.v[unique_old_v_inds,:])
    new_f = np.zeros((len(faceinds), 3), dtype=int)

    for vi in unique_old_v_inds:
        new_f[old_faces==vi] = v_ind_map[vi]

    return simplemesh.SimpleMesh(new_v, new_f)

def set_1ring_faces(sm):
    """
    Create a dict of the adjacent faces of every face in sm
    """
    print('setting 1-ring for faces')
    
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
        reg_front = set([min(remaining_faces),])  # region seed face
        reg_faces = [min(remaining_faces),]
        remaining_faces.remove(min(remaining_faces)) # remove seed face from remaining
        reg_nfaces = 1
            

        # while current region is below max size
        while reg_nfaces<maxfaces:
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

    # partition mesh by connected regions
    region_faces, face_labels = partition_regions(sm, np.inf)
    print('found {} regions'.format(len(region_faces)))
    
    # get largest region
    largest_reg = None
    largest_reg_nfaces = 0
    for rn, rf in region_faces.items():
        if len(rf)>largest_reg_nfaces:
            largest_reg = rn
            largest_reg_nfaces = len(rf)

    print('keeping largest region with {} faces'.format(largest_reg_nfaces))
    
    # create new mesh with just the largest region
    largest_region_mesh = make_sub_mesh(sm, region_faces[largest_reg])
    return largest_region_mesh

def partition_mesh(sm, maxfaces, minfaces):
    """
    Partitions sm into regions with upper and lower faces bounds
    """

    region_faces, face_labels = partition_regions(sm, maxfaces)
    print('merging {} regions'.format(len(region_faces)))
    merge_regions(sm, region_faces, face_labels, minfaces)
    print('making {} region meshes'.format(len(region_faces)))
    region_sms = make_region_meshes(sm, region_faces)
    return region_sms