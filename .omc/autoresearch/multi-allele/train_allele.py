#!/usr/bin/env python3
"""
Train Deep FFN for a specific HLA allele.

Usage: python train_allele.py <allele> <output_model.h5> [n_peptides=50000]
"""
import sys, random, numpy as np
import tensorflow as tf
from tensorflow import keras
from collections import Counter

RANDOM_SEED = 20260619
random.seed(RANDOM_SEED); np.random.seed(RANDOM_SEED); tf.random.set_seed(RANDOM_SEED)

AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")

# BLOSUM62
BLOSUM62 = np.array([
    [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],
    [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],
    [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],
    [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],
    [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
    [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],
    [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],
    [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],
    [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2, 2, 2,-3],
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

def encode(peptide):
    vec = []
    for aa in peptide:
        if aa not in AA_ORDER: aa = 'A'
        idx = AA_ORDER.index(aa)
        vec.extend(BLOSUM62_NORM[idx].tolist())
    return np.array(vec, dtype=np.float32)

def build_deep_ffn(input_dim=180, n_classes=3):
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
                  loss='categorical_crossentropy',
                  metrics=['accuracy'])
    return model

# =============================================================================
# Main
# =============================================================================
if __name__ == "__main__":
    allele = sys.argv[1] if len(sys.argv) > 1 else "HLA-A*01:01"
    out_path = sys.argv[2] if len(sys.argv) > 2 else f"models/deep_ffn_{allele.replace('*','_').replace(':','_')}.h5"
    n_target = int(sys.argv[3]) if len(sys.argv) > 3 else 50000

    print(f"Generating {n_target} random 9-mers for {allele}...")
    peps = set()
    while len(peps) < n_target:
        s = ''.join(random.choice(AA_LIST) for _ in range(9))
        peps.add(s)
    peps = list(peps)

    print(f"Predicting with MHCflurry...")
    from mhcflurry import Class1AffinityPredictor
    predictor = Class1AffinityPredictor.load()
    affinities = predictor.predict(peptides=peps, alleles=[allele]*len(peps))

    # Label using affinity thresholds
    labels = np.array([0 if a > 500 else (1 if a >= 50 else 2) for a in affinities])
    cnt = Counter(labels)
    print(f"Distribution: SB={cnt[2]}, WB={cnt[1]}, NB={cnt[0]}")

    # Balance: downsample to min class size
    min_n = min(cnt[0], cnt[1], cnt[2])
    print(f"Balancing to {min_n} per class ({min_n*3} total)")

    balanced_idx = []
    for c in [0, 1, 2]:
        class_idx = [i for i, l in enumerate(labels) if l == c]
        sample_idx = random.sample(class_idx, min_n)
        balanced_idx.extend(sample_idx)
    random.shuffle(balanced_idx)

    peptides_subset = [peps[i] for i in balanced_idx]
    labels_subset = [labels[i] for i in balanced_idx]
    X = np.array([encode(p) for p in peptides_subset], dtype=np.float32)
    y = keras.utils.to_categorical(labels_subset, 3)

    n_train = int(len(X) * 0.9)
    X_train, X_test = X[:n_train], X[n_train:]
    y_train, y_test = y[:n_train], y[n_train:]

    print(f"Training: {n_train}, Test: {len(X_test)}")
    model = build_deep_ffn()
    model.summary()

    callbacks = [
        keras.callbacks.EarlyStopping(monitor='val_loss', patience=15, restore_best_weights=True),
        keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=5, min_lr=1e-6),
    ]

    history = model.fit(X_train, y_train, epochs=150, batch_size=64,
                        validation_split=0.1, callbacks=callbacks, verbose=1)

    test_loss, test_acc = model.evaluate(X_test, y_test, verbose=0)
    print(f"Test accuracy: {test_acc:.4f}")

    model.save(out_path)
    print(f"Model saved to {out_path}")
    print("Done.")
