#!/usr/bin/env python3
"""Optimized Deep FFN training with Keras 3 (TF 2.21)."""
import os, sys, time, json
import numpy as np
import pandas as pd

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"

# ── Load labeled data ──
df = pd.read_csv("02_Data/raw/real_peptides.csv")
print(f"Loaded: {len(df):,} peptides ({sum(df.data_type=='train'):,} train, {sum(df.data_type=='test'):,} test)")

# ── BLOSUM62 encoding (vectorized numpy) ──
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

# Row-normalize
row_min = BLOSUM62.min(axis=1, keepdims=True)
row_max = BLOSUM62.max(axis=1, keepdims=True)
BLOSUM_NORM = (BLOSUM62 - row_min) / (row_max - row_min + 1e-8)

t0 = time.perf_counter()
peptides = df["peptide"].values
n = len(peptides)
X = np.zeros((n, 180), dtype=np.float32)
for j in range(9):
    aa_chars = [p[j] for p in peptides]
    aa_idx = [AA2IDX[c] for c in aa_chars]
    X[:, j*20:(j+1)*20] = BLOSUM_NORM[aa_idx]
print(f"Encoded {n:,} peptides in {time.perf_counter()-t0:.1f}s")

# Split
train_mask = (df["data_type"] == "train").values
X_train, X_test = X[train_mask], X[~train_mask]
y_train = df["label_num"].values[train_mask].astype(np.int32)
y_test  = df["label_num"].values[~train_mask].astype(np.int32)
print(f"X_train: {X_train.shape}  X_test: {X_test.shape}")

# ── Build Deep FFN ──
import keras
from keras import layers

model = keras.Sequential([
    layers.Input(shape=(180,)),
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

model.compile(
    loss="sparse_categorical_crossentropy",
    optimizer=keras.optimizers.Adam(0.001),
    metrics=["accuracy"],
)
model.summary()

# ── Train ──
t_train = time.perf_counter()
history = model.fit(
    X_train, y_train,
    epochs=150, batch_size=50,
    validation_split=0.2,
    callbacks=[
        keras.callbacks.EarlyStopping(
            monitor="val_loss", patience=10, restore_best_weights=True),
        keras.callbacks.ReduceLROnPlateau(
            monitor="val_loss", factor=0.5, patience=5),
    ],
    verbose=1,
)
train_time = time.perf_counter() - t_train
print(f"Training: {train_time:.0f}s ({len(history.history['loss'])} epochs)")

# ── Evaluate ──
test_loss, test_acc = model.evaluate(X_test, y_test, verbose=0)
print(f"Test accuracy: {test_acc*100:.1f}% (loss: {test_loss:.4f})")

# Per-class F1
y_pred = np.argmax(model.predict(X_test, verbose=0), axis=1)
from sklearn.metrics import classification_report
print(classification_report(y_test, y_pred, target_names=["NB", "WB", "SB"], digits=3))

# ── Save ──
model.save("03_Analysis/models/FFN_Deep_optimized.keras")
print(f"Saved: 03_Analysis/models/FFN_Deep_optimized.keras")
print(f"Total time: {time.perf_counter() - t0:.0f}s")
