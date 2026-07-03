"""Shared BLOSUM62 encoding for Peptide Epitope project.

ponytail: single module replaces 3 copies of BLOSUM matrix + normalization + encoding.
"""
import numpy as np

AA = list("ARNDCQEGHILKMFPSTWYV")
AA2IDX = {aa: i for i, aa in enumerate(AA)}

BLOSUM62 = np.array([
    [4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0],
    [-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3],
    [-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3],
    [-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3],
    [0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
    [-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2],
    [-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2],
    [0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3],
    [-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3],
    [-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3],
    [-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1],
    [-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2],
    [-2,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1],
    [-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1],
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2],
    [1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2],
    [0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0],
    [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3],
    [-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1],
    [0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4]
], dtype=np.float32)

# Row-normalize to [0,1]
_row_min = BLOSUM62.min(axis=1, keepdims=True)
_row_max = BLOSUM62.max(axis=1, keepdims=True)
BLOSUM_NORM = (BLOSUM62 - _row_min) / (_row_max - _row_min + 1e-8)


def encode_blosum(peptides, out=None):
    """Encode peptides as BLOSUM62 feature vectors (vectorized numpy).

    Args:
        peptides: list of 9-mer strings, length n
        out: optional pre-allocated (n, 180) float32 array

    Returns:
        (n, 180) float32 array, 9 positions × 20 amino acids
    """
    n = len(peptides)
    if out is None:
        out = np.empty((n, 180), dtype=np.float32)

    # Integer-encode all peptides at once: (n, 9)
    chars = np.array([list(p) for p in peptides])  # (n, 9)
    idx = np.vectorize(AA2IDX.get)(chars)          # (n, 9) int

    # Index BLOSUM_NORM in one shot: idx(n,9) → values(n,9,20) → reshape(n,180)
    out[:] = BLOSUM_NORM[idx].reshape(n, 180)
    return out


def encode_blosum_image(peptides, out=None):
    """Encode as (n, 9, 20, 1) image tensor for CNN."""
    flat = encode_blosum(peptides)
    return flat.reshape(len(peptides), 9, 20, 1)


def encode_blosum_lstm(peptides, out=None):
    """Encode as (n, 9, 20) for LSTM."""
    flat = encode_blosum(peptides)
    return flat.reshape(len(peptides), 9, 20)


# ── PSSM utilities ──

PSSM_A0201 = {
    "p1": dict(A=0.1, R=0.0, N=0.0, D=0.0, C=0.0, Q=0.0, E=0.0, G=0.2, H=0.0, I=0.2, L=0.1, K=0.0, M=0.1, F=0.4, P=0.0, S=0.2, T=0.2, W=0.3, Y=0.5, V=0.1),
    "p2": dict(A=0.5, R=-1, N=-1, D=-1, C=-1, Q=0.3, E=-1, G=0.0, H=-1, I=0.7, L=1.0, K=-1, M=0.9, F=0.0, P=-1, S=0.0, T=0.3, W=-1, Y=0.0, V=0.7),
    "p3": dict(A=0.0, R=0.0, N=0.0, D=0.3, E=0.3, Q=0.0, G=0.0, H=0.0, I=0.0, L=0.0, K=0.0, M=0.0, F=0.0, P=0.0, S=0.0, T=0.0, C=-0.2, W=0.0, Y=0.0, V=0.0),
    "p4": dict(A=0.1, R=0.1, N=0.0, D=0.0, E=0.0, Q=0.1, G=0.0, H=0.1, I=0.1, L=0.1, K=0.1, M=0.0, F=0.0, P=-0.2, S=0.1, T=0.1, C=-0.1, W=0.0, Y=0.0, V=0.1),
    "p5": dict(A=0.0, R=0.1, N=0.0, D=0.0, E=0.0, Q=0.0, G=0.0, H=0.1, I=0.0, L=0.1, K=0.0, M=0.0, F=0.0, P=-0.1, S=0.0, T=0.0, C=-0.1, W=0.0, Y=0.0, V=0.0),
    "p6": dict(A=0.0, R=0.0, N=0.0, D=0.0, E=0.1, Q=0.0, G=0.0, H=0.0, I=0.1, L=0.1, K=0.0, M=0.0, F=0.0, P=-0.1, S=0.0, T=0.0, C=0.0, W=0.0, Y=0.0, V=0.0),
    "p7": dict(A=0.0, R=0.0, N=0.0, D=0.0, E=0.0, Q=0.0, G=0.0, H=0.0, I=0.0, L=0.1, K=0.0, M=0.0, F=0.1, P=-0.1, S=0.0, T=0.0, C=0.0, W=0.0, Y=0.0, V=0.0),
    "p8": dict(A=0.0, R=0.0, N=0.0, D=0.0, E=0.0, Q=0.0, G=0.0, H=0.0, I=0.1, L=0.2, K=0.0, M=0.0, F=0.0, P=-0.1, S=0.0, T=0.0, C=0.0, W=0.0, Y=0.0, V=0.1),
    "p9": dict(A=0.4, R=-1, N=-1, D=-1, C=-1, Q=0.2, E=-1, G=0.0, H=-1, I=0.7, L=0.8, K=-1, M=0.5, F=0.0, P=-1, S=0.1, T=0.3, W=-1, Y=0.0, V=1.0),
}

# Pre-built (20, 9) PSSM array for vectorized scoring
_pssm_arr = np.zeros((20, 9), dtype=np.float32)
for pos in range(9):
    for aa_char, score in PSSM_A0201[f"p{pos+1}"].items():
        _pssm_arr[AA2IDX[aa_char], pos] = score


def score_pssm_vectorized(peptides):
    """Score peptides with HLA-A*02:01 PSSM (vectorized, no Python loops).

    Args:
        peptides: list of 9-mer strings

    Returns:
        (n,) float64 array of PSSM scores
    """
    n = len(peptides)
    chars = np.array([list(p) for p in peptides])        # (n, 9)
    idx = np.vectorize(AA2IDX.get)(chars)                 # (n, 9) int
    return _pssm_arr[idx[:, 0], 0] + _pssm_arr[idx[:, 1], 1] + \
           _pssm_arr[idx[:, 2], 2] + _pssm_arr[idx[:, 3], 3] + \
           _pssm_arr[idx[:, 4], 4] + _pssm_arr[idx[:, 5], 5] + \
           _pssm_arr[idx[:, 6], 6] + _pssm_arr[idx[:, 7], 7] + \
           _pssm_arr[idx[:, 8], 8]
