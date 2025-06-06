try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

_local_matrices_path = str(pkg_resources.files(__name__) / "sample_matrices")
_error_matrix_path = str(pkg_resources.files(__name__) / "float64_pinv_error_matrix.csv")