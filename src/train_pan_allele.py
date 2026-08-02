#!/usr/bin/env python3
"""Train a pan-specific Deep FFN with allele one-hot encoding in the input.

Hypothesis: a single model with allele identity appended to the peptide encoding
matches or beats the 5 per-allele models, at the cost of one model instead of N.

Design:
  - Input: BLOSUM62 9-mer encoding (180-dim) + 5-dim allele one-hot = 185-dim.
  - Same Deep FFN architecture as per-allele models (360-180-90-45, BN + dropout).
  - Balanced across alleles: subsample each allele to the min allele's class counts
    so no single allele dominates. Combined set = per-allele balanced sets.

Usage:
  Phase 1 (label): handled by existing per-allele data (data/train_HLA-*.npz).
  Phase 2 (train, venv): python src/train_pan_allele.py --out models/deep_ffn_pan.keras
  Phase 3 (eval): python src/evaluate_allele.py models/deep_ffn_pan.keras <allele>
"""
import argparse, os
import numpy as np

RANDOM_SEED = 42
ALLELES = ["HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01", "HLA-B*07:02", "HLA-B*44:03"]
SAFE = lambda a: a.replace("*", "_").replace(":", "")

def load_balanced():
    """Load per-allele npz, subsample each allele to min allele count, add one-hot."""
    rng = np.random.RandomState(RANDOM_SEED)
    sets = []
    min_allele_n = None
    for allele in ALLELES:
        f = f"data/train_{SAFE(allele)}_big.npz" if allele == "HLA-B*07:02" and os.path.exists("data/train_HLA-B_0702_big.npz") else f"data/train_{SAFE(allele)}.npz"
        d = np.load(f, allow_pickle=True)
        X, y = d["X"], d["y"]
        sets.append((allele, X, y))
        if min_allele_n is None or len(X) < min_allele_n:
            min_allele_n = len(X)
    # Downsample each allele to min_allele_n (per-allele stratified would be ideal;
    # these are already class-balanced per allele, so random subsample preserves balance)
    X_parts, y_parts, A_parts = [], [], []
    for i, (allele, X, y) in enumerate(sets):
        idx = rng.choice(len(X), min_allele_n, replace=False)
        X_parts.append(X[idx]); y_parts.append(y[idx])
        aoh = np.zeros(len(ALLELES), dtype=np.float32)
        aoh[i] = 1.0
        A_parts.append(np.repeat(aoh[None, :], min_allele_n, axis=0))
    X = np.concatenate(X_parts); y = np.concatenate(y_parts); A = np.concatenate(A_parts)
    # shuffle globally
    perm = rng.permutation(len(X))
    X, y, A = X[perm], y[perm], A[perm]
    return X, y, A

def build_pan_ffn(input_dim=185, n_classes=3):
    import keras
    model = keras.Sequential([
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
    model.compile(optimizer=keras.optimizers.Adam(0.001),
                  loss='sparse_categorical_crossentropy', metrics=['accuracy'])
    return model

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="models/deep_ffn_pan.keras")
    args = ap.parse_args()

    import keras
    X, y, A = load_balanced()
    print(f"Combined: X{X.shape} (per-allele min={len(X)//5}), class dist "
          f"SB/WB/NB={int((y==2).sum())}/{int((y==1).sum())}/{int((y==0).sum())}")
    Xa = np.concatenate([X, A], axis=1).astype(np.float32)

    rng = np.random.RandomState(RANDOM_SEED)
    idx = np.arange(len(Xa)); rng.shuffle(idx)
    n_train = int(len(Xa) * 0.9)
    tr, te = idx[:n_train], idx[n_train:]
    Xtr, ytr = Xa[tr], y[tr]
    Xte, yte = Xa[te], y[te]
    print(f"Train {len(tr)}, Test {len(te)}")

    model = build_pan_ffn(input_dim=Xa.shape[1])
    callbacks = [
        keras.callbacks.EarlyStopping(monitor='val_loss', patience=15, restore_best_weights=True),
        keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=5, min_lr=1e-6),
    ]
    model.fit(Xtr, ytr, epochs=150, batch_size=64, validation_split=0.1,
              callbacks=callbacks, verbose=1)
    loss, acc = model.evaluate(Xte, yte, verbose=0)
    print(f"Test accuracy: {acc:.4f}")
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    model.save(args.out)
    print(f"Saved {args.out}")

if __name__ == "__main__":
    main()
