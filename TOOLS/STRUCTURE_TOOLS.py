import numpy as np


def ensure_numpy_array(vec):
    if not isinstance(vec, np.ndarray):
        vec = np.array(vec, dtype=np.float64)
    return vec



