"""Fixtures shared across the test suite."""

import numpy
import pytest

from bilby.core.utils import random as bilby_random


@pytest.fixture
def seed():
    """Seed both random number generators the population samplers draw from.

    bilby's priors and numpy's shuffles both feed :meth:`Mass.sample` and
    :meth:`Spin.sample`, so seeding one of them leaves the draw dependent on
    test order. Returns a callable so a test can seed with its own value::

        def test_something(seed):
            seed(4)
    """

    def _seed(value):
        numpy.random.seed(value)
        bilby_random.seed(value)

    return _seed
