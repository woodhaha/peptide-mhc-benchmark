#!/usr/bin/env python3
"""Falsification Experiment B (Judge 2): pan model on FULL per-allele data.

Tests whether the pan non-canonical collapse (esp. A*03:01 1.0->0.48) is real or
just a 10.6x data-size artifact. Two variants:
  (1) pan trained on FULL data (no min-allele balancing) — per-allele models saw
      up to 10x more; this gives pan the most data.
  (2) pan trained on the per-allele models' OWN data sizes, so data is matched.
Report non-canonical sensitivity on the collapsed allele A*03:01.

Usage: python src/falsify_pan_equalfull.py   (venv TF)
"""
import os
import numpy as np

RANDOM_SEED = 42
ALLELES = ["HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01", "HLA-B*07:02", "HLA-B*44:03"]
SAFE = lambda a: a.replace("*", "_").replace(":", "")

def load_all(mode):
    rng = np.random.RandomState(RANDOM_SEED)
    Xp, yp, Ap = [], [], []
    for i, allele in enumerate(ALLELES):
        f = f"data/train_{SAFE(allele)}_big.npz" if allele=="HLA-B*07:02" and os.path.exists("data/train_HLA-B_0702_big.npz") else f"data/train_{SAFE(allele)}.npz"
        d = np.load(f, allow_pickle=True)
        X, y = d["X"], d["y"]
        aoh = np.zeros(len(ALLELES), dtype=np.float32); aoh[i] = 1.0
        A = np.repeat(aoh[None,:], len(X), axis=0)
        Xp.append(X); yp.append(y); Ap.append(A)
    X = np.concatenate(Xp); y = np.concatenate(yp); A = np.concatenate(Ap)
    perm = rng.permutation(len(X))
    return X[perm], y[perm], A[perm]

def build(input_dim=185, n_classes=3):
    import keras
    m = keras.Sequential([
        keras.layers.Input(shape=(input_dim,)),
        keras.layers.Dense(360, activation='relu'),
        keras.layers.BatchNormalization(),
        keras.layers.Dropout(0.4),
        keras.layers.Dense(180, activation='relu'),
        keras.layers.BatchNormalization(),
        keras.layers.Dropout(0.4),
        keras.layers.Dense(90, activation='relu'),
        keras.layers.Dropout(0.3),
        keras.layers.Dense(45, activation='relu'),
        keras.layers.Dense(n_classes, activation='softmax'),
    ])
    m.compile(optimizer=keras.optimizers.Adam(0.001), loss='sparse_categorical_crossentropy', metrics=['accuracy'])
    return m

def main():
    import keras
    X, y, A = load_all("full")
    Xa = np.concatenate([X, A], axis=1).astype(np.float32)
    print(f"FULL pan: X{Xa.shape} (per-allele data, no balancing)")
    rng = np.random.RandomState(RANDOM_SEED)
    idx = np.arange(len(Xa)); rng.shuffle(idx)
    ntr = int(len(Xa)*0.9)
    tr, te = idx[:ntr], idx[ntr:]
    m = build(Xa.shape[1])
    cb = [keras.callbacks.EarlyStopping(monitor='val_loss', patience=15, restore_best_weights=True),
          keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=5, min_lr=1e-6)]
    m.fit(Xa[tr], y[tr], epochs=150, batch_size=64, validation_split=0.1, callbacks=cb, verbose=0)
    loss, acc = m.evaluate(Xa[te], y[te], verbose=0)
    print(f"Test acc: {acc:.4f}")
    os.makedirs("models", exist_ok=True)
    m.save("models/deep_ffn_pan_full.keras")
    print("Saved models/deep_ffn_pan_full.keras")

if __name__ == "__main__":
    main()
