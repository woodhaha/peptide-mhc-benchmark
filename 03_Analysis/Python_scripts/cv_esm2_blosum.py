#!/usr/bin/env python3
"""
5-fold stratified CV: Deep FFN with BLOSUM62 vs ESM-2 t6 per-position embeddings.

Uses real_peptides.csv (5,190 peptides, SB/WB/NB balanced).
Faithful to manuscript training protocol:
  - validation_split=0.2 within training data for early stopping (no test-data leakage)
  - Same architecture, hyperparameters
  - Stratified 5-fold CV + original pre-defined split replication

Cache: 03_Analysis/models/{blosum62_embeddings.npy, esm2_t6_embeddings.npy}
Output: 02_Data/cleaned/esm2_cv_results.csv
"""
import os, sys, time, warnings
import numpy as np
import pandas as pd

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
warnings.filterwarnings("ignore")
BASE = "D:/Researching/Peptide epitope"
sys.path.insert(0, os.path.join(BASE, "03_Analysis"))
from blosum_utils import encode_blosum

print("=" * 65)
print("Deep FFN CV: ESM-2 t6 vs BLOSUM62 (proper val-split protocol)")
print("=" * 65)

# ── Load data ────────────────────────────────────────────────
print("Loading data ...")
df = pd.read_csv(os.path.join(BASE, "02_Data/raw/real_peptides.csv"))
peptides = df["peptide"].values
labels = df["label_num"].values
data_type = df["data_type"].values
n = len(df)
print(f"  {n} peptides  SB={sum(labels==2)}  WB={sum(labels==1)}  NB={sum(labels==0)}")
print(f"  Pre-defined split: train={sum(data_type=='train')}  test={sum(data_type=='test')}")

blosum_cache = os.path.join(BASE, "03_Analysis/models/blosum62_embeddings.npy")
if os.path.exists(blosum_cache):
    X_blosum = np.load(blosum_cache)
else:
    t0 = time.perf_counter()
    X_blosum = encode_blosum(peptides)
    np.save(blosum_cache, X_blosum)
    print(f"BLOSUM62 encoded: {X_blosum.shape} ({time.perf_counter()-t0:.1f}s)")

esm_cache = os.path.join(BASE, "03_Analysis/models/esm2_t6_embeddings.npy")
X_esm = np.load(esm_cache)
print(f"BLOSUM62: {X_blosum.shape} | ESM-2: {X_esm.shape}")

# ── Model / training ─────────────────────────────────────────
import tensorflow as tf
import keras
from keras import layers

def make_blosum_ffn(name="BLOSUM"):
    return keras.Sequential([
        layers.Input(shape=(180,)),
        layers.Dense(360, activation="relu"), layers.BatchNormalization(),
        layers.Dropout(0.5),
        layers.Dense(180, activation="relu"), layers.BatchNormalization(),
        layers.Dropout(0.4),
        layers.Dense(90, activation="relu"), layers.Dropout(0.3),
        layers.Dense(45, activation="relu"),
        layers.Dense(3, activation="softmax"),
    ], name=name)

def make_esm_ffn(name="ESM2"):
    """Manuscript Section 3.9: input(2880)→256→128→64→3"""
    return keras.Sequential([
        layers.Input(shape=(2880,)),
        layers.Dense(256, activation="relu"), layers.BatchNormalization(),
        layers.Dropout(0.5),
        layers.Dense(128, activation="relu"), layers.BatchNormalization(),
        layers.Dropout(0.4),
        layers.Dense(64, activation="relu"), layers.Dropout(0.3),
        layers.Dense(3, activation="softmax"),
    ], name=name)

def train_and_eval(model, X_train, y_train, X_val, y_val, verbose=0):
    """Train with internal validation_split=0.2 for early stopping (no leakage)."""
    model.compile(loss="sparse_categorical_crossentropy",
                  optimizer=keras.optimizers.Adam(0.001), metrics=["accuracy"])
    model.fit(X_train, y_train, epochs=150, batch_size=50,
              validation_split=0.2,
              callbacks=[
                  keras.callbacks.EarlyStopping(monitor="val_loss", patience=10,
                                                restore_best_weights=True),
                  keras.callbacks.ReduceLROnPlateau(monitor="val_loss", factor=0.5,
                                                    patience=5),
              ], verbose=verbose)
    y_pred = np.argmax(model.predict(X_val, verbose=0), axis=1)
    from sklearn.metrics import accuracy_score, f1_score
    return {"accuracy": round(accuracy_score(y_val, y_pred) * 100, 2),
            "macro_f1": round(f1_score(y_val, y_pred, average="macro"), 4),
            "SB_F1": round(f1_score(y_val, y_pred, labels=[2], average="macro"), 4),
            "WB_F1": round(f1_score(y_val, y_pred, labels=[1], average="macro"), 4),
            "NB_F1": round(f1_score(y_val, y_pred, labels=[0], average="macro"), 4)}

def print_metrics(label, met):
    print(f"  {label:<14}  Acc={met['accuracy']:>5.2f}%  MacroF1={met['macro_f1']:.4f}  "
          f"SB={met['SB_F1']:.4f}  WB={met['WB_F1']:.4f}  NB={met['NB_F1']:.4f}")

# ── Part 1: Original pre-defined split ───────────────────────
print("\n─── Part 1: Original pre-defined split (data_type) ───")
train_mask = data_type == "train"
for enc_name, X, make_fn in [
    ("BLOSUM62", X_blosum, make_blosum_ffn),
    ("ESM2_t6",  X_esm,    make_esm_ffn),
]:
    keras.backend.clear_session()
    met = train_and_eval(make_fn(f"Orig_{enc_name}"),
                          X[train_mask], labels[train_mask],
                          X[~train_mask], labels[~train_mask], verbose=0)
    print_metrics(enc_name, met)

# ── Part 2: 5-fold Stratified CV ─────────────────────────────
print("\n─── Part 2: 5-Fold Stratified Cross-Validation ───")
from sklearn.model_selection import StratifiedKFold

skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
results = []
hdr = f"{'Fold':>5}  {'Encoding':<14}  {'Acc%':>7}  {'MacroF1':>8}  {'SB_F1':>7}  {'WB_F1':>7}  {'NB_F1':>7}"
print(hdr); print("-" * len(hdr))

t_start = time.perf_counter()
for fold, (tr, va) in enumerate(skf.split(X_blosum, labels), 1):
    fold_t = time.perf_counter()
    for enc_name, X, make_fn in [
        ("BLOSUM62", X_blosum, make_blosum_ffn),
        ("ESM2_t6",  X_esm,    make_esm_ffn),
    ]:
        keras.backend.clear_session()
        met = train_and_eval(make_fn(f"CV_{enc_name}_f{fold}"),
                              X[tr], labels[tr], X[va], labels[va])
        results.append({"model": "DeepFFN", "encoding": enc_name, "fold": fold, **met})
        print(f"  {fold:>3}  {enc_name:<14}  {met['accuracy']:>5.2f}%  {met['macro_f1']:>8.4f}  "
              f"{met['SB_F1']:>7.4f}  {met['WB_F1']:>7.4f}  {met['NB_F1']:>7.4f}")
    print(f"         fold {fold} done ({time.perf_counter()-fold_t:.0f}s)")

# ── Summary ──────────────────────────────────────────────────
print(f"\n─── Summary (total: {time.perf_counter()-t_start:.0f}s) ───")
results_df = pd.DataFrame(results)
for enc in ["BLOSUM62", "ESM2_t6"]:
    sub = results_df[results_df["encoding"] == enc]
    print(f"\n{enc}:")
    for col in ["accuracy", "macro_f1", "SB_F1", "WB_F1", "NB_F1"]:
        accs = sub[col].values
        print(f"  {col}: {accs.mean():.2f} +- {accs.std():.2f}  [{accs.min():.2f}, {accs.max():.2f}]")

out_path = os.path.join(BASE, "02_Data/cleaned/esm2_cv_results.csv")
results_df.to_csv(out_path, index=False)
print(f"\nSaved: {out_path}")
print("Done.")
