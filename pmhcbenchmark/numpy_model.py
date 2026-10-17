"""Tensorflow-free inference from exported weights (pmhc_models.npz)."""
import numpy as np

from .predictor import ALLELES, SAFE, encode

EPS = 1e-3  # Keras BatchNormalization default epsilon


class NumpyModel:
    """A per-allele Deep FFN running pure-numpy inference.

    Architecture: Dense(relu) -> BN -> Dense(relu) -> BN -> Dense(relu)
    -> Dense(relu) -> Dense(softmax), with dropout inactive at inference.
    """

    def __init__(self, denses, bns):
        self.denses = denses  # list of (W, b)
        self.bns = bns        # list of (gamma, beta, mean, var)

    def predict_proba(self, X):
        x = X.astype(np.float32)
        for i in range(2):
            W, b = self.denses[i]
            x = np.maximum(x @ W + b, 0)
            g, be, mean, var = self.bns[i]
            x = (x - mean) / np.sqrt(var + EPS) * g + be
        for i in (2, 3):
            W, b = self.denses[i]
            x = np.maximum(x @ W + b, 0)
        W, b = self.denses[4]
        x = x @ W + b
        x = np.exp(x - x.max(-1, keepdims=True))
        x /= x.sum(-1, keepdims=True)
        return x

    def score(self, peptides):
        """Return P(weak binder) + P(strong binder) for each peptide."""
        X = np.array([encode(p) for p in peptides], dtype=np.float32)
        p = self.predict_proba(X)
        return p[:, 1] + p[:, 2]


def load_numpy_models(path="pmhc_models.npz"):
    data = np.load(path)
    models = {}
    for allele in ALLELES:
        key = SAFE(allele)
        denses = [(data[f"{key}_d{i}_W"], data[f"{key}_d{i}_b"]) for i in range(5)]
        bns = [(data[f"{key}_bn{i}_g"], data[f"{key}_bn{i}_b"],
                data[f"{key}_bn{i}_mean"], data[f"{key}_bn{i}_var"]) for i in range(2)]
        models[allele] = NumpyModel(denses, bns)
    return models
