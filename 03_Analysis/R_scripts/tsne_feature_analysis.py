#!/usr/bin/env python3
"""t-SNE visualization and feature attribution analysis for Deep FFN model.

Extracts penultimate-layer (45-dim) features from the trained Deep FFN model,
runs t-SNE dimensionality reduction, and computes integrated gradients to
reveal what the model learns at each peptide position.

Outputs:
  03_Analysis/figures/Figure6B_tsne_by_class.png
  03_Analysis/figures/Figure6C_tsne_by_accuracy.png
  03_Analysis/figures/FigureS6_anchor_attention_heatmap.png
"""

import os, sys, warnings
import numpy as np
import pandas as pd

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
os.environ["TF_ENABLE_ONEDNN_OPTS"] = "0"
warnings.filterwarnings("ignore")

import keras
import tensorflow as tf
from sklearn.manifold import TSNE
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from blosum_utils import encode_blosum, AA

OUTDIR = "03_Analysis/figures"
MODEL_PATH = "03_Analysis/models/Deep_FFN.keras"
DATA_PATH = "02_Data/raw/real_peptides.csv"
os.makedirs(OUTDIR, exist_ok=True)

# ── Nature NPG Palette ──
NPG = {
    "red": "#E64B35", "blue": "#4DBBD5", "green": "#00A087",
    "navy": "#3C5488", "pink": "#F39B7F", "grey": "#8491B4",
    "orange": "#DC0000", "brown": "#7E6148", "purple": "#B09C85",
}
CLASS_COLORS = [NPG["navy"], NPG["red"], NPG["green"]]  # NB, WB, SB
CLASS_LABELS = ["Non-Binder", "Weak Binder", "Strong Binder"]

# ══════════════════════════════════════════════════════════════
# Step 1: Load model & data
# ══════════════════════════════════════════════════════════════
print("Loading model ...")
model = keras.saving.load_model(MODEL_PATH)
print(f"  Model: {model.name} ({len(model.layers)} layers)")

# Penultimate layer is dense_3 (45-dim, just before softmax output).
# Strip the last softmax layer to get the 45-dim penultimate layer model.
FEAT_LAYER = "dense_3"
feat_model = keras.Sequential(model.layers[:-1], name="feature_extractor")

print("Loading data + BLOSUM62 encoding ...")
df = pd.read_csv(DATA_PATH)
peptides = df["peptide"].values
X = encode_blosum(peptides)
y = df["label_num"].values.astype(np.int32)
train_mask = (df["data_type"] == "train").values

X_test, y_test = X[~train_mask], y[~train_mask]
X_train, y_train = X[train_mask], y[train_mask]
print(f"  Train: {X_train.shape[0]}, Test: {X_test.shape[0]}")

# ══════════════════════════════════════════════════════════════
# Step 2: Extract penultimate-layer features
# ══════════════════════════════════════════════════════════════
print(f"Extracting {FEAT_LAYER} features ...")
test_feats = feat_model.predict(X_test, verbose=0)  # (n_test, 45)
print(f"  Test features shape: {test_feats.shape}")

# Predictions + correctness
y_pred_probs = model.predict(X_test, verbose=0)
y_pred = np.argmax(y_pred_probs, axis=1)
correct = (y_pred == y_test).astype(int)

# ══════════════════════════════════════════════════════════════
# Step 3: t-SNE on test-set features
# ══════════════════════════════════════════════════════════════
print("Running t-SNE (perplexity=30, 1000 iters) ...")
tsne = TSNE(n_components=2, perplexity=30, max_iter=1000,
            random_state=42, verbose=0)
coords = tsne.fit_transform(test_feats)
print(f"  t-SNE done: {coords.shape}")

# ══════════════════════════════════════════════════════════════
# Step 4: Plot – t-SNE colored by true class
# ══════════════════════════════════════════════════════════════
print("Plotting Figure6B_tsne_by_class ...")
fig, ax = plt.subplots(figsize=(8, 6))
for cls in range(3):
    mask = y_test == cls
    ax.scatter(coords[mask, 0], coords[mask, 1],
               c=CLASS_COLORS[cls], label=CLASS_LABELS[cls],
               s=14, alpha=0.7, edgecolors="none")

ax.set_title("t-SNE of Penultimate-Layer Features\nColored by True Class (SB / WB / NB)",
             fontsize=13, fontweight="bold", pad=12)
ax.set_xlabel("t-SNE Dimension 1", fontsize=11)
ax.set_ylabel("t-SNE Dimension 2", fontsize=11)
ax.legend(loc="best", fontsize=10, markerscale=2)
ax.set_facecolor("white")
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "Figure6B_tsne_by_class.png"), dpi=300)
fig.savefig(os.path.join(OUTDIR, "Figure6B_tsne_by_class.pdf"))
plt.close(fig)
print("  Saved Figure6B_tsne_by_class.{png,pdf}")

# ══════════════════════════════════════════════════════════════
# Step 5: Plot – t-SNE colored by prediction correctness
# ══════════════════════════════════════════════════════════════
print("Plotting Figure6C_tsne_by_accuracy ...")
fig, ax = plt.subplots(figsize=(8, 6))
colors_correct = [NPG["green"], NPG["red"]]
labels_correct = ["Correct", "Incorrect"]
for val in [1, 0]:
    mask = correct == val
    ax.scatter(coords[mask, 0], coords[mask, 1],
               c=colors_correct[val], label=labels_correct[val],
               s=14, alpha=0.7, edgecolors="none")

# Also mark by symbol for class
ax.set_title("t-SNE of Penultimate-Layer Features\nColored by Prediction Correctness",
             fontsize=13, fontweight="bold", pad=12)
ax.set_xlabel("t-SNE Dimension 1", fontsize=11)
ax.set_ylabel("t-SNE Dimension 2", fontsize=11)
ax.legend(loc="best", fontsize=10, markerscale=2)
ax.set_facecolor("white")
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "Figure6C_tsne_by_accuracy.png"), dpi=300)
fig.savefig(os.path.join(OUTDIR, "Figure6C_tsne_by_accuracy.pdf"))
plt.close(fig)
print("  Saved Figure6C_tsne_by_accuracy.{png,pdf}")

# ══════════════════════════════════════════════════════════════
# Step 6: Integrated Gradients – Anchor Position Attention
# ══════════════════════════════════════════════════════════════
print("Computing Integrated Gradients (input × gradient via tf.GradientTape) ...")

n_test = min(X_test.shape[0], 200)  # subset for speed
rng = np.random.RandomState(42)
sample_idx = rng.choice(X_test.shape[0], n_test, replace=False)

ig_attributions = np.zeros((n_test, 9, 20), dtype=np.float32)
for i, idx in enumerate(sample_idx):
    x_i = X_test[idx:idx+1].astype(np.float32)
    y_i = int(y_test[idx])
    x_var = tf.Variable(x_i)
    with tf.GradientTape() as tape:
        preds = model(x_var, training=False)
        loss = preds[0, y_i]
    grads = tape.gradient(loss, x_var).numpy()[0]
    ig_attributions[i] = (x_i[0] * grads).reshape(9, 20)

# Average by true class
class_attn = np.zeros((3, 9, 20), dtype=np.float32)
for cls in range(3):
    cls_mask = y_test[sample_idx] == cls
    if cls_mask.sum() > 0:
        class_attn[cls] = ig_attributions[cls_mask].mean(axis=0)

# Normalize per class for visualization
for cls in range(3):
    vals = class_attn[cls]
    vmin, vmax = vals.min(), vals.max()
    if vmax > vmin:
        class_attn[cls] = (vals - vmin) / (vmax - vmin)

# ── Plot: per-position × amino-acid heatmap ──
print("Plotting FigureS6_anchor_attention_heatmap ...")
fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)

for cls, (ax, title) in enumerate(zip(axes, CLASS_LABELS)):
    im = ax.imshow(class_attn[cls].T, aspect="auto", cmap="YlOrRd",
                   interpolation="nearest")
    ax.set_xticks(range(9))
    ax.set_xticklabels([f"p{i+1}" for i in range(9)], fontsize=10)
    ax.set_yticks(range(20))
    ax.set_yticklabels(AA, fontsize=8)
    ax.set_title(title, fontsize=12, fontweight="bold", color=CLASS_COLORS[cls])
    if cls == 0:
        ax.set_ylabel("Amino Acid", fontsize=11)
    ax.set_xlabel("Peptide Position", fontsize=11)

fig.suptitle("Integrated Gradients — Per-Position × AA Attention by Class",
             fontsize=14, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "FigureS6_anchor_attention_heatmap.png"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(OUTDIR, "FigureS6_anchor_attention_heatmap.pdf"), bbox_inches="tight")
plt.close(fig)
print("  Saved FigureS6_anchor_attention_heatmap.{png,pdf}")

# ── Bonus: Position-wise summary barplot ──
print("Plotting position-wise importance barplot ...")
fig, ax = plt.subplots(figsize=(9, 4))
pos_imp = np.abs(class_attn).sum(axis=2)  # (3, 9) — sum abs over AA
x = np.arange(9)
width = 0.25
for cls in range(3):
    ax.bar(x + cls * width, pos_imp[cls], width,
           label=CLASS_LABELS[cls], color=CLASS_COLORS[cls], alpha=0.85)
ax.set_xticks(x + width)
ax.set_xticklabels([f"p{i+1}" for i in range(9)], fontsize=10)
ax.set_title("Per-Position Feature Importance by Class (Integrated Gradients)",
             fontsize=13, fontweight="bold")
ax.set_ylabel("Mean |Attribution| (sum over AAs)", fontsize=11)
ax.set_xlabel("Peptide Position", fontsize=11)
ax.legend(fontsize=9)
ax.set_facecolor("white")
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "FigureS6_position_importance_bar.png"), dpi=300)
fig.savefig(os.path.join(OUTDIR, "FigureS6_position_importance_bar.pdf"))
plt.close(fig)
print("  Saved FigureS6_position_importance_bar.{png,pdf}")

print("\n═══ All t-SNE feature analysis figures complete ═══")
print(f"  Output: {OUTDIR}/")
