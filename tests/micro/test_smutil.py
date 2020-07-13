import unittest
from typing import List

import numpy as np
from generator import generate, generator

from gias2.mesh.simplemesh import SimpleMesh
from gias2.mesh.smutils import make_sub_mesh

# the new face indices for the new sub mesh
my_faceinds = [3, 4, 5]

my_faceinds2 = []

# the faces and vertices for our original simplemesh
my_sm_f = np.array([[0, 1, 2], [1, 2, 3], [3, 4, 5], [1, 2, 5], [2, 3, 1], [3, 2, 1]])
my_sm_v = np.array(
    [[-127.3, -195.1, -1425.3],
     [-126.8, -195.8, -1425.1],
     [-126.7, -195.2, -1425.0],
     [-127.7, -194.8, -1425.5],
     [-127.3, -195.2, -1425.4],
     [-127.2, -194.8, -1425.3]])

# the faces and vertices for our target sub mesh
my_sm_target_f = np.array([[0, 1, 3], [1, 2, 0], [2, 1, 0]])
my_sm_target_v = np.array(
    [[-126.8, -195.8, -1425.1],
     [-126.7, -195.2, -1425.],
     [-127.7, -194.8, -1425.5],
     [-127.2, -194.8, -1425.3]])


@generator
class SmUtilTestSuite(unittest.TestCase):

    @generate(
        ('test_vectorised_mapping', my_faceinds, my_sm_f, my_sm_v, my_sm_target_f, my_sm_target_v),
        ('test_vectorised_mapping_empty_list', my_faceinds2, my_sm_f, my_sm_v, my_sm_target_f, my_sm_target_v)
    )
    def test_make_sub_mesh(
            self,
            name: str,
            faceinds: List[int],
            sm_f: np.ndarray,
            sm_v: np.ndarray,
            sm_target_f: np.ndarray,
            sm_target_v: np.ndarray) -> None:
        """
        Test make_sub_mesh to ensure mapping of old vertices to new is correct

        :param name: the name of test
        :param faceinds: the face indices for the new sub mesh
        :param sm_f: the faces of our original simplemesh
        :param sm_v: the vertices of our original simplemesh
        :param sm_target_f: the faces of our target simplemesh
        :param sm_target_v: the faces of our original simplemesh
        """
        # create both the original simplemesh and the target sub mesh
        sm = SimpleMesh(sm_v, sm_f)
        sm_target = SimpleMesh(sm_target_v, sm_target_f)

        if len(faceinds) == 0:
            with self.assertRaises(ValueError) as context:
                # make our new sub mesh
                sm_new = make_sub_mesh(sm, faceinds)
                self.assertTrue('length of face_indices is zero' in str(context.exception))
        else:
            # compare new sub mesh to target sub mesh
            # make our new sub mesh
            sm_new = make_sub_mesh(sm, faceinds)
            self.assertTrue(np.all(sm_new.v == sm_target.v) and np.all(sm_new.f == sm_target.f))
