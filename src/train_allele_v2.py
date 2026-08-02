#!/usr/bin/env python3
"""Train Deep FFN per HLA allele with fixed data splits + mixed biased/random data.

Two-phase pipeline (matches restored toolchain):
  1. Labeling: runs under BASE ANACONDA python (mhcflurry+torch work there).
  2. Training: runs under the r-tensorflow venv python (TF/keras work there).

Fixed seed = 42 and datasets saved to disk BEFORE training (fixes the
"training variance" problem identified in the multi-allele decision log).

Usage:
  # Phase 1 (label) — base python:
  python src/train_allele_v2.py label HLA-A*01:01 --n 50000 --out data/train_HLA-A_0101.npz
  # Phase 2 (train) — venv python:
  python src/train_allele_v2.py train --data data/train_HLA-A_0101.npz --out models/deep_ffn_HLA-A_0101.keras
"""
import argparse, os, random, sys
import numpy as np

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")
RANDOM_SEED = 42
MIX_RATIO = 0.5  # fraction of biased (motif-matched) peptides in training pool

AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
BLOSUM62 = np.array([
    [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],
    [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],
    [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],
    [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],
    [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
    [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],
    [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],
    [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],
    [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3],
    [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3],
    [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1],
    [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2],
    [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1],
    [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1],
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2],
    [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2],
    [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0],
    [-3,-3,-4,-4,-2,-2,-3,-2, 2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3],
    [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1],
    [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4],
], dtype=np.float32)
BLOSUM62_NORM = np.array([(row-row.min())/(row.max()-row.min()+1e-8) for row in BLOSUM62])

# Canonical anchor positions per allele (index 0-based)
ANCHORS = {
    "HLA-A*01:01": {2: "DE", 8: "YF"},
    "HLA-A*03:01": {1: "LVIM", 8: "KR"},
    "HLA-B*44:03": {1: "ED", 8: "FWYL"},
    "HLA-B*07:02": {2: "P", 8: "LIMF"},     # P2 proline dominant, P9 aliphatic
    "HLA-A*02:01": {1: "LMIVAT", 8: "LVIM"},
}

def encode(peptide):
    vec = []
    for aa in peptide:
        if aa not in AA_ORDER: aa = 'A'
        vec.extend(BLOSUM62_NORM[AA_ORDER.index(aa)].tolist())
    return np.array(vec, dtype=np.float32)

def biased_peptide(allele, rng):
    """Generate a peptide matching the allele's canonical anchor motif."""
    anchors = ANCHORS[allele]
    p = [rng.choice(AA_LIST) for _ in range(9)]
    for pos, allowed in anchors.items():
        p[pos] = rng.choice(list(allowed))
    return ''.join(p)

def make_training_data(allele, n_target, out_path):
    rng = random.Random(RANDOM_SEED)
    np_rng = np.random.RandomState(RANDOM_SEED)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    n_biased = int(n_target * MIX_RATIO)
    n_random = n_target - n_biased
    peps = set()
    while len(peps) < n_biased:
        peps.add(biased_peptide(allele, rng))
    while len(peps) < n_target:
        s = ''.join(rng.choice(AA_LIST) for _ in range(9))
        peps.add(s)
    peps = sorted(peps)  # deterministic order
    print(f"Generated {len(peps)} peptides ({n_biased} biased, {n_random} random)")

    from mhcflurry import Class1AffinityPredictor
    predictor = Class1AffinityPredictor.load()
    print(f"Predicting {len(peps)} affinities for {allele}...")
    affinities = predictor.predict(peptides=peps, alleles=[allele]*len(peps))
    # MHCflurry affinity threshold convention: SB<=50nM, WB<=500nM, NB>500nM
    labels = np.array([2 if a <= 50 else (1 if a <= 500 else 0) for a in affinities])

    # Balance to min class count for a deterministic balanced set
    from collections import Counter
    cnt = Counter(labels.tolist())
    min_n = min(cnt[0], cnt[1], cnt[2])
    print(f"Distribution SB/WB/NB: {cnt[2]}/{cnt[1]}/{cnt[0]}, balancing to {min_n} each")
    idx_by_class = {c: [i for i,l in enumerate(labels) if l==c] for c in (0,1,2)}
    chosen = []
    for c in (0,1,2):
        chosen.extend(np_rng.choice(idx_by_class[c], min_n, replace=False).tolist())
    np_rng.shuffle(chosen)
    peps_s = [peps[i] for i in chosen]
    labs = labels[chosen]

    X = np.array([encode(p) for p in peps_s], dtype=np.float32)
    y = labs.astype(np.int64)
    np.savez(out_path, X=X, y=y, peptides=np.array(peps_s, dtype=object))
    print(f"Saved {out_path}: X{X.shape}, class dist SB/WB/NB={int((y==2).sum())}/{int((y==1).sum())}/{int((y==0).sum())}")
    return out_path

def build_deep_ffn(input_dim=180, n_classes=3):
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

def train(data_path, out_path):
    import keras
    npz = np.load(data_path, allow_pickle=True)
    X, y = npz['X'], npz['y']
    peps = npz['peptides'].tolist()
    n = len(X)
    # Fixed stratified 90/10 split (seed 42), computed deterministically
    rng = np.random.RandomState(RANDOM_SEED)
    idx = np.arange(n)
    rng.shuffle(idx)
    n_train = int(n * 0.9)
    tr, te = idx[:n_train], idx[n_train:]
    X_train, y_train = X[tr], y[tr]
    X_test, y_test = X[te], y[te]
    print(f"Train {len(tr)}, Test {len(te)}")

    model = build_deep_ffn()
    callbacks = [
        keras.callbacks.EarlyStopping(monitor='val_loss', patience=15, restore_best_weights=True),
        keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=5, min_lr=1e-6),
    ]
    model.fit(X_train, y_train, epochs=150, batch_size=64, validation_split=0.1,
              callbacks=callbacks, verbose=1)
    test_loss, test_acc = model.evaluate(X_test, y_test, verbose=0)
    print(f"Test accuracy: {test_acc:.4f}")
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    model.save(out_path)
    print(f"Model saved to {out_path}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["label", "train"])
    ap.add_argument("allele", nargs="?", default=None)
    ap.add_argument("--n", type=int, default=50000)
    ap.add_argument("--out", default=None)
    ap.add_argument("--data", default=None)
    args = ap.parse_args()

    if args.mode == "label":
        assert args.allele and args.out
        make_training_data(args.allele, args.n, args.out)
    else:
        assert args.data and args.out
        train(args.data, args.out)

if __name__ == "__main__":
    main()
