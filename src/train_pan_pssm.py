#!/usr/bin/env python3
"""Train pan-specific Deep FFN with ALLELE PSSM conditioning (not bare one-hot).

Hypothesis (tests the one-hot-too-weak mechanism): conditioning on a per-allele
amino-acid preference profile (PSSM) instead of a 5-dim one-hot carries the
structural signal needed to retain non-canonical binding robustness that the
bare-one-hot pan model lost.

Design:
  - Input: BLOSUM62 9-mer (180-dim) + allele PSSM (5 alleles x 20-dim = 100-dim).
    The PSSM is a learned-ish representation of each allele's anchor preferences,
    derived from the same motif data used for biased training.
  - Same Deep FFN arch; total input 280-dim.

Allele PSSM construction (deterministic, from motif): for each allele, a 20-dim
profile weighting each anchor position's allowed residues 1.0 and others ~0, then
softmax-normalized. This is a hand-crafted but principled conditioning signal.
"""
import argparse, os
import numpy as np

RANDOM_SEED = 42
ALLELES = ["HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01", "HLA-B*07:02", "HLA-B*44:03"]
SAFE = lambda a: a.replace("*", "_").replace(":", "")
# canonical anchor positions: {allele: {pos1indexed: allowed_residues}}
ANCHORS = {
    "HLA-A*02:01": {2: "LMIVAT", 9: "LVIM"},
    "HLA-A*01:01": {3: "DE", 9: "YF"},
    "HLA-A*03:01": {2: "LVIM", 9: "KR"},
    "HLA-B*07:02": {2: "P", 9: "LIMF"},
    "HLA-B*44:03": {2: "ED", 9: "FWYL"},
}
AA20 = "ACDEFGHIKLMNPQRSTVWY"

def allele_pssm(allele):
    """20-dim profile: 1.0 for residues appearing at any anchor, softmax-scaled."""
    allowed = set()
    for pos, res in ANCHORS[allele].items():
        allowed.update(res)
    v = np.array([1.0 if aa in allowed else 0.05 for aa in AA20], dtype=np.float32)
    v /= v.sum()
    return v

def load_balanced():
    rng = np.random.RandomState(RANDOM_SEED)
    sets = []
    min_n = None
    for allele in ALLELES:
        f = f"data/train_{SAFE(allele)}_big.npz" if allele == "HLA-B*07:02" and os.path.exists("data/train_HLA-B_0702_big.npz") else f"data/train_{SAFE(allele)}.npz"
        d = np.load(f, allow_pickle=True)
        sets.append((allele, d["X"], d["y"]))
        if min_n is None or len(d["X"]) < min_n: min_n = len(d["X"])
    Xp, yp, Pp = [], [], []
    for i, (allele, X, y) in enumerate(sets):
        idx = rng.choice(len(X), min_n, replace=False)
        Xp.append(X[idx]); yp.append(y[idx])
        pssm = allele_pssm(allele)
        Pp.append(np.repeat(pssm[None, :], min_n, axis=0))
    X = np.concatenate(Xp); y = np.concatenate(yp); P = np.concatenate(Pp)
    perm = rng.permutation(len(X))
    return X[perm], y[perm], P[perm]

def build(input_dim=280, n_classes=3):
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
    ap.add_argument("--out", default="models/deep_ffn_pan_pssm.keras")
    args = ap.parse_args()
    import keras
    X, y, P = load_balanced()
    print(f"X{X.shape}, PSSM{P.shape}, class SB/WB/NB={int((y==2).sum())}/{int((y==1).sum())}/{int((y==0).sum())}")
    Xa = np.concatenate([X, P], axis=1).astype(np.float32)
    rng = np.random.RandomState(RANDOM_SEED)
    idx = np.arange(len(Xa)); rng.shuffle(idx)
    ntr = int(len(Xa)*0.9)
    tr, te = idx[:ntr], idx[ntr:]
    model = build(input_dim=Xa.shape[1])
    callbacks = [keras.callbacks.EarlyStopping(monitor='val_loss', patience=15, restore_best_weights=True),
                 keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=5, min_lr=1e-6)]
    model.fit(Xa[tr], y[tr], epochs=150, batch_size=64, validation_split=0.1, callbacks=callbacks, verbose=1)
    loss, acc = model.evaluate(Xa[te], y[te], verbose=0)
    print(f"Test accuracy: {acc:.4f}")
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    model.save(args.out)
    print(f"Saved {args.out}")

if __name__ == "__main__":
    main()
