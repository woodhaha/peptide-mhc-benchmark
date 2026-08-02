#!/usr/bin/env python3
"""Scale the rebalanced pan model to a 16-allele panel.

Labels + trains per-allele models for the new alleles, then trains an expanded
rebalanced pan model. Mirrors train_allele_v2.py / train_pan_multiseed.py but for
the expanded allele set.

Two phases (base python labels, venv trains) — see CLI.

Usage:
  label:  python src/scale_panel.py label <allele> --out data/panel_<allele>.npz
  pan:    python src/scale_panel.py pan --out models/deep_ffn_pan16.keras
"""
import argparse, os, random, sys
import numpy as np

PANEL_ANCHORS = {
    # existing 5
    "HLA-A*02:01": {1: "LMIVAT", 8: "LVIM"},
    "HLA-A*01:01": {2: "DE", 8: "YF"},
    "HLA-A*03:01": {1: "LVIM", 8: "KR"},
    "HLA-B*07:02": {1: "P", 8: "LIMF"},
    "HLA-B*44:03": {1: "ED", 8: "FWYL"},
    # new 11 (canonical anchors from literature/supertype)
    "HLA-A*11:01": {1: "LVIM", 8: "KR"},
    "HLA-A*24:02": {1: "YF", 8: "FYLIM"},
    "HLA-B*08:01": {2: "RK", 8: "L"},
    "HLA-B*15:01": {1: "Q", 8: "FYWL"},
    "HLA-B*27:05": {1: "RK", 8: "FYWLIM"},
    "HLA-B*35:01": {1: "P", 8: "YFMWL"},
    "HLA-B*40:01": {1: "E", 8: "L"},
    "HLA-B*49:01": {1: "E", 8: "LIM"},
    "HLA-B*51:01": {1: "P", 8: "ILVW"},
    "HLA-B*57:01": {1: "ATS", 8: "FWYL"},
}
AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")
RANDOM_SEED = 42
MIX_RATIO = 0.5
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

def encode(pep):
    v = []
    for aa in pep:
        if aa not in AA_ORDER: aa = 'A'
        v.extend(BLOSUM62_NORM[AA_ORDER.index(aa)].tolist())
    return np.array(v, dtype=np.float32)

def biased_peptide(allele, rng):
    anchors = PANEL_ANCHORS[allele]
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
    print(f"Generated {len(peps)} peptides ({n_biased} biased)")
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

def build(input_dim, n_classes=3):
    import keras
    m = keras.Sequential([
        keras.layers.Input(shape=(input_dim,)),
        keras.layers.Dense(360, activation='relu'), keras.layers.BatchNormalization(), keras.layers.Dropout(0.4),
        keras.layers.Dense(180, activation='relu'), keras.layers.BatchNormalization(), keras.layers.Dropout(0.4),
        keras.layers.Dense(90, activation='relu'), keras.layers.Dropout(0.3),
        keras.layers.Dense(45, activation='relu'),
        keras.layers.Dense(n_classes, activation='softmax')])
    m.compile(optimizer=keras.optimizers.Adam(0.001), loss='sparse_categorical_crossentropy', metrics=['accuracy'])
    return m

def train_pan(alleles, out_path, min_per_allele=2000):
    import keras
    rng = np.random.RandomState(RANDOM_SEED)
    Xp, yp, Ap = [], [], []
    for i, allele in enumerate(alleles):
        # new alleles saved as train_HLA-A1101.npz (no underscore), existing as train_HLA-A_0101.npz
        f = f"data/train_{SAFE(allele)}.npz"
        if not os.path.exists(f):
            f = f"data/train_{allele.replace('*','').replace(':','')}.npz"
        d = np.load(f, allow_pickle=True); X, y = d["X"], d["y"]
        # upsample thin alleles to min_per_allele (avoids a tiny allele dominating min-balance)
        if len(X) < min_per_allele:
            idx = rng.choice(len(X), min_per_allele, replace=True)
            X, y = X[idx], y[idx]
        aoh = np.zeros(len(alleles), dtype=np.float32); aoh[i] = 1.0
        A = np.repeat(aoh[None,:], len(X), axis=0)
        Xp.append(X); yp.append(y); Ap.append(A)
    X=np.concatenate(Xp); y=np.concatenate(yp); A=np.concatenate(Ap)
    perm=rng.permutation(len(X)); X,y,A=X[perm],y[perm],A[perm]
    Xa=np.concatenate([X,A],axis=1).astype(np.float32)
    print(f"Pan16: X{Xa.shape}")
    rng2=np.random.RandomState(RANDOM_SEED); idx=np.arange(len(Xa)); rng2.shuffle(idx)
    ntr=int(len(Xa)*0.9); tr,te=idx[:ntr],idx[ntr:]
    m=build(Xa.shape[1])
    cb=[keras.callbacks.EarlyStopping(monitor='val_loss',patience=15,restore_best_weights=True),
        keras.callbacks.ReduceLROnPlateau(monitor='val_loss',factor=0.5,patience=5,min_lr=1e-6)]
    m.fit(Xa[tr],y[tr],epochs=150,batch_size=64,validation_split=0.1,callbacks=cb,verbose=1)
    acc=m.evaluate(Xa[te],y[te],verbose=0)[1]
    print(f"Pan16 test acc: {acc:.4f}")
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    m.save(out_path); print(f"Saved {out_path}")

SAFE = lambda a: a.replace("*", "_").replace(":", "")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["label","pan"]); ap.add_argument("allele", nargs="?", default=None)
    ap.add_argument("--n", type=int, default=50000); ap.add_argument("--out", default=None)
    args = ap.parse_args()
    if args.mode=="label":
        assert args.allele and args.out
        make_training_data(args.allele, args.n, args.out)
    else:
        train_pan(list(PANEL_ANCHORS), args.out or "models/deep_ffn_pan16.keras")

if __name__ == "__main__":
    main()
