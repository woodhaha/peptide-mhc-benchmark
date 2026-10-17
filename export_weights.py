#!/usr/bin/env python3
"""Export the per-allele Deep FFN weights to pmhc_models.npz for
tensorflow-free inference (GUI / .exe). Run once; the .npz is all the GUI needs.

Usage: python export_weights.py   # -> pmhc_models.npz
"""
import numpy as np
import keras

from pmhcbenchmark.predictor import ALLELES, SAFE


def main():
    out = {}
    for allele in ALLELES:
        key = SAFE(allele)
        m = keras.models.load_model(f"models/deep_ffn_{key}.keras", compile=False)
        denses, bns = [], []
        for l in m.layers:
            if isinstance(l, keras.layers.Dense):
                W, b = l.get_weights()
                denses.append((W.astype(np.float32), b.astype(np.float32)))
            elif isinstance(l, keras.layers.BatchNormalization):
                g, be, mean, var = l.get_weights()
                bns.append((g.astype(np.float32), be.astype(np.float32),
                            mean.astype(np.float32), var.astype(np.float32)))
        for i, (W, b) in enumerate(denses):
            out[f"{key}_d{i}_W"] = W
            out[f"{key}_d{i}_b"] = b
        for i, (g, be, mean, var) in enumerate(bns):
            out[f"{key}_bn{i}_g"] = g
            out[f"{key}_bn{i}_b"] = be
            out[f"{key}_bn{i}_mean"] = mean
            out[f"{key}_bn{i}_var"] = var
    np.savez_compressed("pmhc_models.npz", **out)
    print(f"saved pmhc_models.npz ({len(out)} arrays for {len(ALLELES)} alleles)")


if __name__ == "__main__":
    main()
