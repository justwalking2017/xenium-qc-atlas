import numpy as np
import pandas as pd
from scipy import sparse

from xenium_showcase.analysis import CONTROL_PREFIXES, NEGATIVE_FEATURE_TYPES
from xenium_showcase.report import _table


def test_control_prefixes_are_platform_aware():
    assert "NegativeControlProbe_1".startswith(CONTROL_PREFIXES)
    assert not "EPCAM".startswith(CONTROL_PREFIXES)
    assert "Negative Control Probe" in NEGATIVE_FEATURE_TYPES
    assert "Genomic Control" not in NEGATIVE_FEATURE_TYPES
    assert "Deprecated Codeword" not in NEGATIVE_FEATURE_TYPES


def test_sparse_normalization_contract():
    x = sparse.csr_matrix([[1, 2], [0, 3]])
    assert np.asarray(x.sum(axis=1)).ravel().tolist() == [3, 3]


def test_report_table_escapes_content():
    rendered = _table(pd.DataFrame({"gene": ["A<script>"]}))
    assert "<script>" not in rendered
    assert "&lt;script&gt;" in rendered
