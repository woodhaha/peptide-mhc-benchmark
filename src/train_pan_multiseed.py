#!/usr/bin/env python3
"""Train rebalanced pan model across multiple seeds; report mean+/-SD non-canonical sens.

Usage (venv TF): python src/train_pan_multiseed.py --seeds 0,1,42,123,2024
Saves models/deep_ffn_pan_seed{S}.keras, prints per-seed test acc.
"""
import argparse, os
import numpy as np

ALLELES = ["HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01", "HLA-B*07:02", "HLA-B*44:03"]
SAFE = lambda a: a.replace("*", "_").replace(":", "")
B7_FACTOR = 4

def load(seed):
    rng = np.random.RandomState(seed)
    Xp, yp, Ap = [], [], []
    for i, allele in enumerate(ALLELES):
        f = f"data/train_{SAFE(allele)}_big.npz" if allele=="HLA-B*07:02" and os.path.exists("data/train_HLA-B_0702_big.npz") else f"data/train_{SAFE(allele)}.npz"
        d = np.load(f, allow_pickle=True)
        X, y = d["X"], d["y"]
        aoh = np.zeros(len(ALLELES), dtype=np.float32); aoh[i] = 1.0
        A = np.repeat(aoh[None,:], len(X), axis=0)
        if allele == "HLA-B*07:02" and B7_FACTOR > 1:
            idx = rng.choice(len(X), len(X)*B7_FACTOR, replace=True)
            X, y, A = X[idx], y[idx], A[idx]
        Xp.append(X); yp.append(y); Ap.append(A)
    X = np.concatenate(Xp); y = np.concatenate(yp); A = np.concatenate(Ap)
    perm = rng.permutation(len(X))
    return X[perm], y[perm], A[perm]

def build(dim):
    import keras
    m = keras.Sequential([
        keras.layers.Input(shape=(dim,)),
        keras.layers.Dense(360, activation='relu'),
        keras.layers.BatchNormalization(),
        keras.layers.Dropout(0.4),
        keras.layers.Dense(180, activation='relu'),
        keras.layers.BatchNormalization(),
        keras.layers.Dropout(0.4),
        keras.layers.Dense(90, activation='relu'),
        keras.layers.Dropout(0.3),
        keras.layers.Dense(45, activation='relu'),
        keras.layers.Dense(3, activation='softmax'),
    ])
    m.compile(optimizer=keras.optimizers.Adam(0.001), loss='sparse_categorical_crossentropy', metrics=['accuracy'])
    return m

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,42,123,2024")
    ap.add_argument("--outdir", default="models")
    args = ap.parse_args()
    import keras
    keras.utils.set_random_seed(42)  # deterministic base init; per-seed variation via data RNG
    seeds = [int(s) for s in args.seeds.split(",")]
    for seed in seeds:
        X, y, A = load(seed)
        Xa = np.concatenate([X, A], axis=1).astype(np.float32)
        rng = np.random.RandomState(seed)
        idx = np.arange(len(Xa)); rng.shuffle(idx)
        ntr = int(len(Xa)*0.9)
        tr, te = idx[:ntr], idx[ntr:]
        m = build(Xa.shape[1])
        cb = [keras.callbacks.EarlyStopping(monitor='val_loss', patience=15, restore_best_weights=True),
              keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=5, min_lr=1e-6)]
        m.fit(Xa[tr], y[tr], epochs=150, batch_size=64, validation_split=0.1, callbacks=cb, verbose=0)
        acc = m.evaluate(Xa[te], y[te], verbose=0)[1]
        path = os.path.join(args.outdir, f"deep_ffn_pan_seed{seed}.keras")
        m.save(path)
        print(f"seed={seed}: test_acc={acc:.4f} -> {path}")

if __name__ == "__main__":
    main()
