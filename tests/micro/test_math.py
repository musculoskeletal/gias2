import unittest
import warnings
from typing import Sequence

import numpy as np
from generator import generate, generator
from numpy import nan
from numpy.testing import assert_almost_equal

from gias2.common.math import norm


@generator
class TestMath(unittest.TestCase):
    """
        Test the 'math' module
    """

    def setUp(self) -> None:
        super().setUp()

        warnings.filterwarnings(
            'ignore',
            category=RuntimeWarning,
            module=r'gias2\.common\.mxath',
            message="invalid value encountered in true_divide")

    @generate(
        ((1, 1, 1), (0.5773503, 0.5773503, 0.5773503)),
        ((1, -1, 1), (0.5773503, -0.5773503, 0.5773503)),
        ((-1, -1, -1), (-0.5773503, -0.5773503, -0.5773503)),
        ((1, 0, 0), (1, 0, 0)),
        ((0, 1, 0), (0, 1, 0)),
        ((0, 0, 1), (0, 0, 1)),
        # WARNING: this generates a runtime warning (but not an exception)
        #      "RuntimeWarning: invalid value encountered in true_divide"
        ((0, 0, 0), (nan, nan, nan)),
        ((100, 100, 100), (0.5773503, 0.5773503, 0.5773503)),
    )
    def test_norm(self, v: Sequence[float], expected: Sequence[float]) -> None:
        result = norm(np.array(v))
        assert_almost_equal(result, np.array(expected))
