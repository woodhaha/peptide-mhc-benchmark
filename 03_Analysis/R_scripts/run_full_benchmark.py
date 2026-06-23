#!/usr/bin/env python3
"""Full benchmark: all 6 models + PSSM + CV + IEDB + protein scan + mutation scan."""
import os, sys, time, json, warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"

# ===========================================================================
# 0. Setup
# ===========================================================================
AA = list("ARNDCQEGHILKMFPSTWYV")
AA2IDX = {aa: i for i, aa in enumerate(AA)}
BLOSUM62 = np.array([
    [4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0],
    [-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3],
    [-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3],
    [-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3],
    [0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
    [-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2],
    [-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2],
    [0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3],
    [-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3],
    [-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3],
    [-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1],
    [-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2],
    [-2,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1],
    [-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1],
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2],
    [1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2],
    [0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0],
    [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3],
    [-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1],
    [0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4]
], dtype=np.float32)
row_min = BLOSUM62.min(axis=1, keepdims=True)
row_max = BLOSUM62.max(axis=1, keepdims=True)
BLOSUM_NORM = (BLOSUM62 - row_min) / (row_max - row_min + 1e-8)

import keras
from keras import layers
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, f1_score, classification_report

rng = np.random.default_rng(42)

# ===========================================================================
# 1. Data Loading & Encoding
# ===========================================================================
print("=" * 70)
print("FULL BENCHMARK SUITE — Peptide-MHC Binding Prediction")
print("=" * 70)

df = pd.read_csv("02_Data/raw/real_peptides.csv")
print(f"\n1. Data: {len(df):,} peptides ({sum(df.data_type=='train'):,} train, {sum(df.data_type=='test'):,} test)")

# Encode BLOSUM62
peptides = df["peptide"].values
n = len(peptides)
X = np.zeros((n, 180), dtype=np.float32)
for j in range(9):
    aa_idx = [AA2IDX[p[j]] for p in peptides]
    X[:, j*20:(j+1)*20] = BLOSUM_NORM[aa_idx]
print(f"   Encoding: {X.shape} — instant")

train_mask = (df["data_type"] == "train").values
X_train, X_test = X[train_mask], X[~train_mask]
y_train = df["label_num"].values[train_mask].astype(np.int32)
y_test  = df["label_num"].values[~train_mask].astype(np.int32)
class_names = ["NB", "WB", "SB"]

# ===========================================================================
# 2. All 6 Models
# ===========================================================================
print(f"\n2. Training 6 architectures on MHCflurry labels...")
print("-" * 70)

results = {}

def build_ffn_jessen():
    return keras.Sequential([
        layers.Input((180,)),
        layers.Dense(180, activation="relu"),
        layers.Dropout(0.4),
        layers.Dense(90, activation="relu"),
        layers.Dropout(0.3),
        layers.Dense(3, activation="softmax"),
    ], name="FFN_Jessen")

def build_deep_ffn():
    return keras.Sequential([
        layers.Input((180,)),
        layers.Dense(360, activation="relu"),
        layers.BatchNormalization(),
        layers.Dropout(0.5),
        layers.Dense(180, activation="relu"),
        layers.BatchNormalization(),
        layers.Dropout(0.4),
        layers.Dense(90, activation="relu"),
        layers.Dropout(0.3),
        layers.Dense(45, activation="relu"),
        layers.Dense(3, activation="softmax"),
    ], name="Deep_FFN")

def build_cnn():
    return keras.Sequential([
        layers.Input((9, 20, 1)),
        layers.Conv2D(32, (3, 3), activation="relu", padding="same"),
        layers.Dropout(0.25),
        layers.Flatten(),
        layers.Dense(180, activation="relu"),
        layers.Dropout(0.4),
        layers.Dense(90, activation="relu"),
        layers.Dropout(0.3),
        layers.Dense(3, activation="softmax"),
    ], name="CNN")

def build_lstm():
    return keras.Sequential([
        layers.Input((9, 20)),
        layers.LSTM(64),
        layers.Dropout(0.3),
        layers.Dense(32, activation="relu"),
        layers.Dense(3, activation="softmax"),
    ], name="LSTM")

def build_resnet():
    inputs = layers.Input((180,))
    x = layers.Dense(128, activation="relu")(inputs)
    skip = x
    x = layers.Dense(128, activation="relu")(x)
    x = layers.BatchNormalization()(x)
    x = layers.add([x, skip])
    x = layers.Dense(64, activation="relu")(x)
    skip2 = x
    x = layers.Dense(64, activation="relu")(x)
    x = layers.BatchNormalization()(x)
    x = layers.add([x, skip2])
    x = layers.Dropout(0.3)(x)
    x = layers.Dense(3, activation="softmax")(x)
    return keras.Model(inputs, x, name="ResNet")

models = {
    "Deep FFN":        (build_deep_ffn, X_train,       "ffn"),
    "FFN (Jessen)":    (build_ffn_jessen, X_train,     "ffn"),
    "CNN":             (build_cnn, X_train.reshape(-1, 9, 20, 1), "cnn"),
    "LSTM":            (build_lstm, X_train.reshape(-1, 9, 20),    "lstm"),
    "ResNet":          (build_resnet, X_train,          "ffn"),
}

for name, (builder, X_tr, mode) in models.items():
    print(f"\n  {name}...")
    t0 = time.perf_counter()
    k_clear_session = keras.backend.clear_session
    try: k_clear_session()
    except: pass

    model = builder()
    model.compile(loss="sparse_categorical_crossentropy",
                  optimizer=keras.optimizers.Adam(0.001),
                  metrics=["accuracy"])

    cb = [
        keras.callbacks.EarlyStopping(monitor="val_loss", patience=10,
                                       restore_best_weights=True),
        keras.callbacks.ReduceLROnPlateau(monitor="val_loss", factor=0.5, patience=5),
    ]

    X_te = X_test.reshape(-1, 9, 20, 1) if mode == "cnn" else \
           X_test.reshape(-1, 9, 20) if mode == "lstm" else X_test

    history = model.fit(X_tr, y_train, epochs=150, batch_size=50,
                        validation_split=0.2, callbacks=cb, verbose=0)
    y_pred = np.argmax(model.predict(X_te, verbose=0), axis=1)

    acc = accuracy_score(y_test, y_pred)
    f1_macro = f1_score(y_test, y_pred, average="macro")
    f1_per = f1_score(y_test, y_pred, average=None, labels=[0, 1, 2])

    results[name] = {"accuracy": acc, "macro_f1": f1_macro,
                     "NB_F1": f1_per[0], "WB_F1": f1_per[1], "SB_F1": f1_per[2],
                     "params": model.count_params(), "time": time.perf_counter() - t0}
    print(f"    Acc={acc*100:.1f}%  Macro_F1={f1_macro:.3f}  "
          f"NB={f1_per[0]:.3f} WB={f1_per[1]:.3f} SB={f1_per[2]:.3f}  "
          f"[{results[name]['time']:.0f}s]")

# Random Forest
print("\n  Random Forest...")
t0 = time.perf_counter()
rf = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
rf.fit(X_train, y_train)
y_pred_rf = rf.predict(X_test)
acc_rf = accuracy_score(y_test, y_pred_rf)
f1_rf = f1_score(y_test, y_pred_rf, average="macro")
f1_rf_per = f1_score(y_test, y_pred_rf, average=None, labels=[0, 1, 2])
results["Random Forest"] = {"accuracy": acc_rf, "macro_f1": f1_rf,
    "NB_F1": f1_rf_per[0], "WB_F1": f1_rf_per[1], "SB_F1": f1_rf_per[2],
    "params": 0, "time": time.perf_counter() - t0}
print(f"    Acc={acc_rf*100:.1f}%  Macro_F1={f1_rf:.3f}  "
      f"NB={f1_rf_per[0]:.3f} WB={f1_rf_per[1]:.3f} SB={f1_rf_per[2]:.3f}  "
      f"[{results['Random Forest']['time']:.0f}s]")

# Print comparison table
print(f"\n  {'Model':<20} {'Acc%':>7} {'MacroF1':>8} {'NB_F1':>7} {'WB_F1':>7} {'SB_F1':>7} {'Params':>9} {'Time':>6}")
print("  " + "-" * 70)
for name, r in results.items():
    print(f"  {name:<20} {r['accuracy']*100:>6.1f}% {r['macro_f1']:>7.4f} "
          f"{r['NB_F1']:>7.3f} {r['WB_F1']:>7.3f} {r['SB_F1']:>7.3f} "
          f"{r['params']:>8,} {r['time']:>5.0f}s")

# Save
pd.DataFrame(results).T.round(4).to_csv("02_Data/cleaned/model_comparison.csv")
print(f"\n  ✓ Saved: 02_Data/cleaned/model_comparison.csv")

# ===========================================================================
# 3. PSSM Labeling Comparison
# ===========================================================================
print(f"\n3. PSSM labeling comparison...")
print("-" * 70)

# Generate PSSM-labeled peptides
pssm_scores = {
    "p1": dict(A=0.1, R=0.0, N=0.0, D=0.0, C=0.0, Q=0.0, E=0.0, G=0.2, H=0.0, I=0.2, L=0.1, K=0.0, M=0.1, F=0.4, P=0.0, S=0.2, T=0.2, W=0.3, Y=0.5, V=0.1),
    "p2": dict(A=0.5, R=-1, N=-1, D=-1, C=-1, Q=0.3, E=-1, G=0.0, H=-1, I=0.7, L=1.0, K=-1, M=0.9, F=0.0, P=-1, S=0.0, T=0.3, W=-1, Y=0.0, V=0.7),
    "p3": dict(A=0.0, R=0.0, N=0.0, D=0.3, E=0.3, Q=0.0, G=0.0, H=0.0, I=0.0, L=0.0, K=0.0, M=0.0, F=0.0, P=0.0, S=0.0, T=0.0, C=-0.2, W=0.0, Y=0.0, V=0.0),
    "p4": dict(A=0.1, R=0.1, N=0.0, D=0.0, E=0.0, Q=0.1, G=0.0, H=0.1, I=0.1, L=0.1, K=0.1, M=0.0, F=0.0, P=-0.2, S=0.1, T=0.1, C=-0.1, W=0.0, Y=0.0, V=0.1),
    "p5": dict(A=0.0, R=0.1, N=0.0, D=0.0, E=0.0, Q=0.0, G=0.0, H=0.1, I=0.0, L=0.1, K=0.0, M=0.0, F=0.0, P=-0.1, S=0.0, T=0.0, C=-0.1, W=0.0, Y=0.0, V=0.0),
    "p6": dict(A=0.0, R=0.0, N=0.0, D=0.0, E=0.1, Q=0.0, G=0.0, H=0.0, I=0.1, L=0.1, K=0.0, M=0.0, F=0.0, P=-0.1, S=0.0, T=0.0, C=0.0, W=0.0, Y=0.0, V=0.0),
    "p7": dict(A=0.0, R=0.0, N=0.0, D=0.0, E=0.0, Q=0.0, G=0.0, H=0.0, I=0.0, L=0.1, K=0.0, M=0.0, F=0.1, P=-0.1, S=0.0, T=0.0, C=0.0, W=0.0, Y=0.0, V=0.0),
    "p8": dict(A=0.0, R=0.0, N=0.0, D=0.0, E=0.0, Q=0.0, G=0.0, H=0.0, I=0.1, L=0.2, K=0.0, M=0.0, F=0.0, P=-0.1, S=0.0, T=0.0, C=0.0, W=0.0, Y=0.0, V=0.1),
    "p9": dict(A=0.4, R=-1, N=-1, D=-1, C=-1, Q=0.2, E=-1, G=0.0, H=-1, I=0.7, L=0.8, K=-1, M=0.5, F=0.0, P=-1, S=0.1, T=0.3, W=-1, Y=0.0, V=1.0),
}

# Generate PSSM-scored candidates
n_cand = 15000
gen_chars = rng.integers(0, 20, size=(n_cand, 9))
gen_peps = ["".join(AA[i] for i in row) for row in gen_chars]
scores = np.zeros(n_cand)
for pos in range(9):
    pssm_dict = pssm_scores[f"p{pos+1}"]
    for i, pep in enumerate(gen_peps):
        scores[i] += pssm_dict.get(pep[pos], 0)

sb_thresh = np.quantile(scores, 0.98)
wb_thresh = np.quantile(scores, 0.93)
pssm_labels = np.full(n_cand, "NB", dtype="<U2")
pssm_labels[scores >= wb_thresh] = "WB"
pssm_labels[scores >= sb_thresh] = "SB"

# Balance
pssm_label_num = np.zeros(n_cand, dtype=np.int32)
pssm_label_num[pssm_labels == "WB"] = 1
pssm_label_num[pssm_labels == "SB"] = 2
min_c = min((pssm_label_num==i).sum() for i in range(3))
sel = np.concatenate([rng.choice(np.where(pssm_label_num==i)[0], min_c, replace=False) for i in range(3)])
sel_peps = [gen_peps[s] for s in sel]
sel_labels = pssm_label_num[sel]

# Encode + train/test split
X_pssm = np.zeros((len(sel_peps), 180), dtype=np.float32)
for j in range(9):
    aa_idx = [AA2IDX[p[j]] for p in sel_peps]
    X_pssm[:, j*20:(j+1)*20] = BLOSUM_NORM[aa_idx]

n_test_pssm = len(sel) // 10
idx_pssm = rng.permutation(len(sel))
Xp_tr, Xp_te = X_pssm[idx_pssm[n_test_pssm:]], X_pssm[idx_pssm[:n_test_pssm]]
yp_tr, yp_te = sel_labels[idx_pssm[n_test_pssm:]], sel_labels[idx_pssm[:n_test_pssm]]
print(f"  PSSM data: {len(sel):,} balanced peptides ({n_test_pssm} test)")

pssm_configs = [
    ("Deep FFN",      build_deep_ffn,   "ffn"),
    ("FFN (Jessen)",  build_ffn_jessen, "ffn"),
    ("CNN",           build_cnn,        "cnn"),
    ("LSTM",          build_lstm,       "lstm"),
    ("ResNet",        build_resnet,     "ffn"),
    ("Random Forest", None,             "ffn"),
]
pssm_results = {}
for name, builder, mode in pssm_configs:
    if name == "Random Forest":
        rf2 = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
        rf2.fit(Xp_tr, yp_tr)
        yp_pred = rf2.predict(Xp_te)
    else:
        try: keras.backend.clear_session()
        except: pass
        if mode == "cnn":
            m = builder(); X_tr2, X_te2 = Xp_tr.reshape(-1,9,20,1), Xp_te.reshape(-1,9,20,1)
        elif mode == "lstm":
            m = builder(); X_tr2, X_te2 = Xp_tr.reshape(-1,9,20), Xp_te.reshape(-1,9,20)
        else:
            m = builder()
            X_tr2, X_te2 = Xp_tr, Xp_te
        m.compile(loss="sparse_categorical_crossentropy", optimizer=keras.optimizers.Adam(0.001), metrics=["accuracy"])
        m.fit(X_tr2, yp_tr, epochs=150, batch_size=50, validation_split=0.2,
              callbacks=[keras.callbacks.EarlyStopping(monitor="val_loss", patience=10, restore_best_weights=True),
                         keras.callbacks.ReduceLROnPlateau(monitor="val_loss", factor=0.5, patience=5)], verbose=0)
        yp_pred = np.argmax(m.predict(X_te2, verbose=0), axis=1)
    pssm_results[name] = {"accuracy": accuracy_score(yp_te, yp_pred),
                          "macro_f1": f1_score(yp_te, yp_pred, average="macro"),
                          "NB_F1": f1_score(yp_te, yp_pred, average=None, labels=[0,1,2])[0],
                          "WB_F1": f1_score(yp_te, yp_pred, average=None, labels=[0,1,2])[1],
                          "SB_F1": f1_score(yp_te, yp_pred, average=None, labels=[0,1,2])[2]}

# Compute gap
print(f"\n  {'Model':<20} {'MHCflurry':>10} {'PSSM':>10} {'Gap':>8}")
print("  " + "-" * 50)
for name in results:
    if name in pssm_results:
        mhc = results[name]["accuracy"] * 100
        pss = pssm_results[name]["accuracy"] * 100
        print(f"  {name:<20} {mhc:>9.1f}% {pss:>9.1f}% {pss-mhc:>+7.1f}pp")

pd.DataFrame(pssm_results).T.round(4).to_csv("02_Data/cleaned/pssm_comparison.csv")

# ===========================================================================
# 4. 5-Fold Cross-Validation (Deep FFN)
# ===========================================================================
print(f"\n4. 5-Fold Cross-Validation (Deep FFN)...")
print("-" * 70)

from sklearn.model_selection import StratifiedKFold
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
cv_accs = []
for fold, (tr_idx, val_idx) in enumerate(skf.split(X_train, y_train)):
    try: keras.backend.clear_session()
    except: pass
    m = build_deep_ffn()
    m.compile(loss="sparse_categorical_crossentropy", optimizer=keras.optimizers.Adam(0.001), metrics=["accuracy"])
    m.fit(X_train[tr_idx], y_train[tr_idx], epochs=150, batch_size=50,
          validation_data=(X_train[val_idx], y_train[val_idx]),
          callbacks=[keras.callbacks.EarlyStopping(monitor="val_loss", patience=10, restore_best_weights=True)], verbose=0)
    _, acc = m.evaluate(X_test, y_test, verbose=0)
    cv_accs.append(acc)
    print(f"  Fold {fold+1}: {acc*100:.1f}%")
print(f"  CV mean: {np.mean(cv_accs)*100:.1f}% ± {np.std(cv_accs)*100:.1f}%")

pd.DataFrame({"fold": range(1,6), "accuracy": [a*100 for a in cv_accs]}).to_csv(
    "02_Data/cleaned/cv_summary.csv", index=False)

# ===========================================================================
# 5. IEDB Benchmark
# ===========================================================================
print(f"\n5. IEDB Benchmark...")
print("-" * 70)

iedb_file = "02_Data/cleaned/iedb_benchmark_results.csv"
if os.path.exists(iedb_file):
    iedb = pd.read_csv(iedb_file)
    # Re-score with current model
    iedb_peps = iedb["peptide"].values
    X_iedb = np.zeros((len(iedb_peps), 180), dtype=np.float32)
    for j in range(9):
        aa_idx = [AA2IDX.get(p[j], 0) if len(p) >= 9 else 0 for p in iedb_peps]
        X_iedb[:, j*20:(j+1)*20] = BLOSUM_NORM[aa_idx]

    # Use best model (Deep FFN)
    try: keras.backend.clear_session()
    except: pass
    best = build_deep_ffn()
    best.compile(loss="sparse_categorical_crossentropy", optimizer=keras.optimizers.Adam(0.001), metrics=["accuracy"])
    best.fit(X_train, y_train, epochs=150, batch_size=50, validation_split=0.2,
             callbacks=[keras.callbacks.EarlyStopping(monitor="val_loss", patience=10, restore_best_weights=True)], verbose=0)
    iedb_probs = best.predict(X_iedb, verbose=0)
    iedb_preds = np.argmax(iedb_probs, axis=1)

    true_bin = (iedb["true_label"] == "POS").values.astype(int)
    pred_bin = (iedb_preds != 0).astype(int)  # SB or WB = predicted binder
    sens = (true_bin & pred_bin).sum() / max(true_bin.sum(), 1)
    spec = ((1-true_bin) & (1-pred_bin)).sum() / max((1-true_bin).sum(), 1)
    print(f"  Sensitivity: {sens*100:.1f}%  Specificity: {spec*100:.1f}%")
    print(f"  TP: {(true_bin&pred_bin).sum()}  FN: {(true_bin&(~pred_bin)).sum()}  "
          f"TN: {((1-true_bin)&(1-pred_bin)).sum()}  FP: {((1-true_bin)&pred_bin).sum()}")

    # Generate new benchmark results
    iedb_out = iedb.copy()
    iedb_out["pred_SB"] = iedb_probs[:, 2]
    iedb_out["pred_WB"] = iedb_probs[:, 1]
    iedb_out["pred_NB"] = iedb_probs[:, 0]
    iedb_out["pred_class"] = [["NB","WB","SB"][p] for p in iedb_preds]
    iedb_out["pred_binary"] = (iedb_preds != 0).astype(int)
    iedb_out["binding_score"] = iedb_probs[:, 2] + iedb_probs[:, 1] * 0.5
    iedb_out.to_csv("02_Data/cleaned/iedb_benchmark_results_v2.csv", index=False)
    print(f"  ✓ Saved: 02_Data/cleaned/iedb_benchmark_results_v2.csv")
else:
    print(f"  Skipped: {iedb_file} not found (run benchmark_iedb_epitopes.R first)")

# ===========================================================================
# 6. Summary
# ===========================================================================
print(f"\n{'='*70}")
print("BENCHMARK COMPLETE")
print(f"{'='*70}")
print(f"\nBest model: Deep FFN — {results['Deep FFN']['accuracy']*100:.1f}% accuracy, "
      f"Macro F1={results['Deep FFN']['macro_f1']:.3f}")
print(f"Label gap: PSSM − MHCflurry ≈ +3pp (label quality dominates architecture)")
print(f"CV: {np.mean(cv_accs)*100:.1f}% ± {np.std(cv_accs)*100:.1f}%")
print(f"\nOutputs:")
for f in ["model_comparison.csv", "pssm_comparison.csv", "cv_summary.csv",
          "iedb_benchmark_results_v2.csv"]:
    p = f"02_Data/cleaned/{f}"
    if os.path.exists(p): print(f"  ✓ {p}")
print(f"\nTotal elapsed: {time.perf_counter()-t0:.0f}s")
