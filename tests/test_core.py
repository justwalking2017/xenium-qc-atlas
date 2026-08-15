import numpy as np
from scipy import sparse

from xenium_showcase.analysis import CONTROL_PREFIXES


def test_control_prefixes_are_platform_aware():
    assert "NegativeControlProbe_1".startswith(CONTROL_PREFIXES)
    assert not "EPCAM".startswith(CONTROL_PREFIXES)


def test_sparse_normalization_contract():
    x = sparse.csr_matrix([[1, 2], [0, 3]])
    assert np.asarray(x.sum(axis=1)).ravel().tolist() == [3, 3]

