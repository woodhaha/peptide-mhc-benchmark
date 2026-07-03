#!/usr/bin/env python3
"""Optimized Deep FFN training with Keras 3 (TF 2.21)."""
import os, sys, time, json
import numpy as np
import pandas as pd

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"

# Shared BLOSUM62 encoder
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from blosum_utils import AA, AA2IDX, BLOSUM_NORM, encode_blosum

# ── Load labeled data ──
df = pd.read_csv("02_Data/raw/real_peptides.csv")
print(f"Loaded: {len(df):,} peptides ({sum(df.data_type=='train'):,} train, {sum(df.data_type=='test'):,} test)")

# ── BLOSUM62 encoding (vectorized, shared module) ──
t0 = time.perf_counter()
peptides = df["peptide"].values
n = len(peptides)
X = encode_blosum(peptides)
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
