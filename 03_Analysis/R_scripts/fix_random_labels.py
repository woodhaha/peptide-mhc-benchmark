#!/usr/bin/env python3
"""Generate the MISSING random-label baseline (proper i.i.d. null).

The manuscript reported "random synthetic labels 65.8%" but no script/data
reproduced it. This script generates a true i.i.d. random-label null on the
authoritative 5,190-peptide set and trains all 6 architectures (5 seeds each).

Expected result: ~33.3% accuracy (3-class balanced, no learnable signal).
"""
import os, sys, time, json, warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from blosum_utils import encode_blosum

import keras
from keras import layers
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score

SEEDS = [0, 1, 42, 123, 2024]

def clear_session():
    try: keras.backend.clear_session()
    except Exception: pass

def build_ffn_jessen():
    return keras.Sequential([
        layers.Input((180,)), layers.Dense(180, activation="relu"), layers.Dropout(0.4),
        layers.Dense(90, activation="relu"), layers.Dropout(0.3),
        layers.Dense(3, activation="softmax")], name="FFN_Jessen")

def build_deep_ffn():
    return keras.Sequential([
        layers.Input((180,)), layers.Dense(360, activation="relu"), layers.BatchNormalization(), layers.Dropout(0.5),
        layers.Dense(180, activation="relu"), layers.BatchNormalization(), layers.Dropout(0.4),
        layers.Dense(90, activation="relu"), layers.Dropout(0.3),
        layers.Dense(45, activation="relu"), layers.Dense(3, activation="softmax")], name="Deep_FFN")

def build_cnn():
    return keras.Sequential([
        layers.Input((9, 20, 1)), layers.Conv2D(32, (3, 3), activation="relu", padding="same"), layers.Dropout(0.25),
        layers.Flatten(), layers.Dense(180, activation="relu"), layers.Dropout(0.4),
        layers.Dense(90, activation="relu"), layers.Dropout(0.3),
        layers.Dense(3, activation="softmax")], name="CNN")

def build_lstm():
    return keras.Sequential([
        layers.Input((9, 20)), layers.LSTM(64), layers.Dropout(0.3),
        layers.Dense(32, activation="relu"), layers.Dense(3, activation="softmax")], name="LSTM")

def build_resnet():
    inputs = layers.Input((180,))
    x = layers.Dense(128, activation="relu")(inputs); skip = x
    x = layers.Dense(128, activation="relu")(x); x = layers.BatchNormalization()(x); x = layers.add([x, skip])
    x = layers.Dense(64, activation="relu")(x); skip2 = x
    x = layers.Dense(64, activation="relu")(x); x = layers.BatchNormalization()(x); x = layers.add([x, skip2])
    x = layers.Dropout(0.3)(x); x = layers.Dense(3, activation="softmax")(x)
    return keras.Model(inputs, x, name="ResNet")

def reshape(X, mode):
    if mode == "cnn": return X.reshape(-1, 9, 20, 1)
    if mode == "lstm": return X.reshape(-1, 9, 20)
    return X

def train_one(builder, X_tr, y_tr, X_te, y_te, mode, seed):
    keras.utils.set_random_seed(seed); np.random.seed(seed)
    clear_session()
    m = builder()
    m.compile(loss="sparse_categorical_crossentropy", optimizer=keras.optimizers.Adam(0.001), metrics=["accuracy"])
    cb = [keras.callbacks.EarlyStopping(monitor="val_loss", patience=10, restore_best_weights=True),
          keras.callbacks.ReduceLROnPlateau(monitor="val_loss", factor=0.5, patience=5)]
    m.fit(reshape(X_tr, mode), y_tr, epochs=150, batch_size=50, validation_split=0.2, callbacks=cb, verbose=0)
    return accuracy_score(y_te, np.argmax(m.predict(reshape(X_te, mode), verbose=0), axis=1))

def main():
    df = pd.read_csv("02_Data/raw/real_peptides.csv")
    X = encode_blosum(df["peptide"].values)
    train_mask = (df["data_type"] == "train").values
    X_tr, X_te = X[train_mask], X[~train_mask]
    y_true = df["label_num"].values.astype(np.int32)

    models = [("Deep FFN", build_deep_ffn, "ffn"), ("FFN (Jessen)", build_ffn_jessen, "ffn"),
              ("CNN", build_cnn, "cnn"), ("LSTM", build_lstm, "lstm"), ("ResNet", build_resnet, "ffn")]
    results = {name: [] for name, _, _ in models}
    results["Random Forest"] = []

    t0 = time.perf_counter()
    for seed in SEEDS:
        # i.i.d. random labels — shuffle the true labels (preserves class balance, destroys any signal)
        rng = np.random.default_rng(seed)
        y_rand = rng.permutation(y_true)
        y_tr_rand, y_te_rand = y_rand[train_mask], y_rand[~train_mask]
        for name, builder, mode in models:
            acc = train_one(builder, X_tr, y_tr_rand, X_te, y_te_rand, mode, seed)
            results[name].append(acc)
            print(f"  seed={seed:4d}  {name:16s} random-label acc={acc*100:.2f}%", flush=True)
        rf = RandomForestClassifier(n_estimators=100, random_state=seed, n_jobs=-1)
        rf.fit(X_tr, y_tr_rand)
        results["Random Forest"].append(accuracy_score(y_te_rand, rf.predict(X_te)))

    print(f"\n{'='*60}\nRandom-label null (i.i.d., mean ± SD, {len(SEEDS)} seeds)\n{'='*60}")
    out = {}
    for name, accs in sorted(results.items(), key=lambda kv: -np.mean(kv[1])):
        m, s = np.mean(accs), np.std(accs, ddof=1)
        out[name] = {"mean": round(m, 4), "sd": round(s, 4), "accs": [round(a, 4) for a in accs]}
        print(f"  {name:16s}  {m*100:.2f}% ± {s*100:.2f}%")
    print(f"\nelapsed {time.perf_counter()-t0:.0f}s")

    outfile = "02_Data/cleaned/random_label_results.json"
    json.dump(out, open(outfile, "w"), indent=2)
    print(f"saved: {outfile}")

if __name__ == "__main__":
    main()
