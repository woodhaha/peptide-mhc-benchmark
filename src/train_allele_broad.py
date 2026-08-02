#!/usr/bin/env python3
"""Retrain per-allele models with BROAD anchors (Type A residues in biased generator).

Tests whether including Type A tolerated residues (A*02:01 P2-Q, A*03:01 P9-H/Y,
B*07:02 P9-V) in the biased training generator recovers Type A sensitivity (0.48)
toward canonical levels.

Two-phase, mirrors train_allele_v2.py:
  label (base python, mhcflurry):  python src/train_allele_broad.py label <allele> --out ...
  train (venv TF):                 python src/train_allele_broad.py train --data ... --out ...

Broad anchors (0-indexed positions, matching train_allele_v2):
  A*02:01: {1: LMIVATQ, 8: LVIM}
  A*03:01: {1: LVIM, 8: KRYH}
  B*07:02: {1: P, 8: LIMFV}
"""
import argparse, os, random, sys
import numpy as np

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")
RANDOM_SEED = 42
MIX_RATIO = 0.5
AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
BROAD_ANCHORS = {
    "HLA-A*02:01": {1: "LMIVATQ", 8: "LVIM"},
    "HLA-A*03:01": {1: "LVIM", 8: "KRYH"},
    "HLA-B*07:02": {1: "P", 8: "LIMFV"},
}
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

def encode(pep):
    v = []
    for aa in pep:
        if aa not in AA_ORDER: aa = 'A'
        v.extend(BLOSUM62_NORM[AA_ORDER.index(aa)].tolist())
    return np.array(v, dtype=np.float32)

def biased_peptide(allele, rng):
    anchors = BROAD_ANCHORS[allele]
    p = [rng.choice(AA_LIST) for _ in range(9)]
    for pos, allowed in anchors.items():
        p[pos] = rng.choice(list(allowed))
    return ''.join(p)

def make_training_data(allele, n_target, out_path):
    rng = random.Random(RANDOM_SEED); np_rng = np.random.RandomState(RANDOM_SEED)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    n_biased = int(n_target * MIX_RATIO)
    peps = set()
    while len(peps) < n_biased: peps.add(biased_peptide(allele, rng))
    while len(peps) < n_target: peps.add(''.join(rng.choice(AA_LIST) for _ in range(9)))
    peps = sorted(peps)
    print(f"Generated {len(peps)} peptides ({n_biased} broad-biased)")
    from mhcflurry import Class1AffinityPredictor
    predictor = Class1AffinityPredictor.load()
    aff = predictor.predict(peptides=peps, alleles=[allele]*len(peps))
    labels = np.array([2 if a<=50 else (1 if a<=500 else 0) for a in aff])
    from collections import Counter
    cnt = Counter(labels.tolist())
    min_n = min(cnt[0], cnt[1], cnt[2])
    print(f"SB/WB/NB: {cnt[2]}/{cnt[1]}/{cnt[0]}, balancing to {min_n}")
    idx_by = {c:[i for i,l in enumerate(labels) if l==c] for c in (0,1,2)}
    chosen=[]
    for c in (0,1,2): chosen.extend(np_rng.choice(idx_by[c], min_n, replace=False).tolist())
    np_rng.shuffle(chosen)
    X = np.array([encode(peps[i]) for i in chosen], dtype=np.float32)
    y = labels[chosen].astype(np.int64)
    np.savez(out_path, X=X, y=y, peptides=np.array([peps[i] for i in chosen], dtype=object))
    print(f"Saved {out_path}: X{X.shape}")

def build(input_dim=180):
    import keras
    m = keras.Sequential([
        keras.layers.Input(shape=(input_dim,)),
        keras.layers.Dense(360, activation='relu'), keras.layers.BatchNormalization(), keras.layers.Dropout(0.4),
        keras.layers.Dense(180, activation='relu'), keras.layers.BatchNormalization(), keras.layers.Dropout(0.4),
        keras.layers.Dense(90, activation='relu'), keras.layers.Dropout(0.3),
        keras.layers.Dense(45, activation='relu'),
        keras.layers.Dense(3, activation='softmax')])
    m.compile(optimizer=keras.optimizers.Adam(0.001), loss='sparse_categorical_crossentropy', metrics=['accuracy'])
    return m

def train(data_path, out_path):
    import keras
    npz = np.load(data_path, allow_pickle=True); X, y = npz['X'], npz['y']
    rng = np.random.RandomState(RANDOM_SEED)
    idx = np.arange(len(X)); rng.shuffle(idx)
    ntr = int(len(X)*0.9); tr, te = idx[:ntr], idx[ntr:]
    m = build()
    cb = [keras.callbacks.EarlyStopping(monitor='val_loss', patience=15, restore_best_weights=True),
          keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=5, min_lr=1e-6)]
    m.fit(X[tr], y[tr], epochs=150, batch_size=64, validation_split=0.1, callbacks=cb, verbose=1)
    acc = m.evaluate(X[te], y[te], verbose=0)[1]
    print(f"Test acc: {acc:.4f}")
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    m.save(out_path); print(f"Saved {out_path}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["label","train"]); ap.add_argument("allele", nargs="?", default=None)
    ap.add_argument("--n", type=int, default=50000); ap.add_argument("--out", default=None); ap.add_argument("--data", default=None)
    args = ap.parse_args()
    if args.mode=="label":
        assert args.allele in BROAD_ANCHORS, f"{args.allele} not in broad set {list(BROAD_ANCHORS)}"
        make_training_data(args.allele, args.n, args.out)
    else:
        train(args.data, args.out)

if __name__ == "__main__":
    main()
