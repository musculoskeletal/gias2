import unittest

import numpy as np
from generator import generator, generate

from gias2.common.transform3D import directAffine, calcAffineMatrixSVD, transformAffine

point_set_0 = np.array([
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
])
t_matrix_0 = np.eye(4)

point_set_1 = np.array([
    [1, 2, 3],
    [1 + 1, 0 + 2, 0 + 3],
    [0 + 1, 1 + 2, 0 + 3],
    [0 + 1, 0 + 2, 1 + 3],
])
t_matrix_1 = np.array([
    [1, 0, 0, point_set_1[0, 0]],
    [0, 1, 0, point_set_1[0, 1]],
    [0, 0, 1, point_set_1[0, 2]],
    [0, 0, 0, 1],
])

point_set_2 = np.array([
    [0, 0, 0],
    [2, 0, 0],
    [0, 3, 0],
    [0, 0, 0.5],
])
t_matrix_2 = np.array([
    [2, 0, 0, 0],
    [0, 3, 0, 0],
    [0, 0, 0.5, 0],
    [0, 0, 0, 1],
])


@generator
class TestTransform3D(unittest.TestCase):

    @generate(
        ('identity', point_set_0, point_set_0, t_matrix_0),
        ('translation', point_set_0, point_set_1, t_matrix_1),
        ('scaling', point_set_0, point_set_2, t_matrix_2)
    )
    def test_direct_affine(self, name: str, source: np.ndarray, target: np.ndarray, expected_mat: np.ndarray):
        matrix: np.ndarray = directAffine(source, target)
        transformed = transformAffine(source, matrix)
        self.assertTrue(np.allclose(transformed, target, equal_nan=True))
        self.assertTrue(np.allclose(matrix, expected_mat, equal_nan=True))

    @generate(
        ('identity', point_set_0, point_set_0, t_matrix_0),
        ('translation', point_set_0, point_set_1, t_matrix_1),
    )
    def test_calc_affine_matrix_svd(self, name: str, source: np.ndarray, target: np.ndarray, expected_mat: np.ndarray):
        matrix: np.ndarray = calcAffineMatrixSVD(source, target)
        transformed = transformAffine(source, matrix)
        self.assertTrue(np.allclose(transformed, target, equal_nan=True))
        self.assertTrue(np.allclose(matrix, expected_mat, equal_nan=True))
