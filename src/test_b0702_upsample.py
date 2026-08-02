#!/usr/bin/env python3
"""Test if B*07:02 pan over-anchoring is fixable by upsampling B*07:02 in pan training.

B*07:02's per-allele model (606 samples) beats pan-full (0.96 vs 0.60 on non-canonical).
Question: is the gap a data-size artifact (B*07:02 = 606 vs A*02:01 = 13k in the pan
mix) or intrinsic over-anchoring? Variants:
  (A) pan trained with B*07:02 upsampled 4x (to ~2400) — matches its relative scarcity
  (B) pan trained with B*07:02 at native weight but ALL alleles at native size (baseline = pan_full)

Reports B*07:02 non-canonical sens to see if upsampling recovers it.

Usage (venv TF): python src/test_b0702_upsample.py
"""
import os
import numpy as np

RANDOM_SEED = 42
ALLELES = ["HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01", "HLA-B*07:02", "HLA-B*44:03"]
SAFE = lambda a: a.replace("*", "_").replace(":", "")

def load_with_upsample(b7_factor=4):
    rng = np.random.RandomState(RANDOM_SEED)
    Xp, yp, Ap = [], [], []
    for i, allele in enumerate(ALLELES):
        f = f"data/train_{SAFE(allele)}_big.npz" if allele=="HLA-B*07:02" and os.path.exists("data/train_HLA-B_0702_big.npz") else f"data/train_{SAFE(allele)}.npz"
        d = np.load(f, allow_pickle=True)
        X, y = d["X"], d["y"]
        aoh = np.zeros(len(ALLELES), dtype=np.float32); aoh[i] = 1.0
        A = np.repeat(aoh[None,:], len(X), axis=0)
        if allele == "HLA-B*07:02" and b7_factor > 1:
            # repeat B*07:02 b7_factor times (with slight jitter via resample)
            idx = rng.choice(len(X), len(X)*b7_factor, replace=True)
            X, y, A = X[idx], y[idx], A[idx]
            print(f"  B*07:02 upsampled {len(X)} (x{b7_factor})")
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
    X, y, A = load_with_upsample(b7_factor=4)
    Xa = np.concatenate([X, A], axis=1).astype(np.float32)
    print(f"Pan+B7x4: X{Xa.shape}")
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
    m.save("models/deep_ffn_pan_b7x4.keras")
    print("Saved models/deep_ffn_pan_b7x4.keras")

if __name__ == "__main__":
    main()
