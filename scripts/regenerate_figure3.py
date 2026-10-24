#!/usr/bin/env python3
"""Regenerate Figure 3 (IEDB benchmark) with the corrected 42-epitope set.

AUC 0.970 (was 0.953); confusion matrix [[15,5],[1,41]] on 42 A*02:01-verified
positives + 20 homopolymer negatives. Publication-style matplotlib.
"""
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import roc_curve, auc, confusion_matrix

plt.rcParams.update({
    "font.family": "Arial", "font.size": 9, "axes.linewidth": 0.8,
    "axes.labelsize": 9, "axes.titlesize": 9, "xtick.labelsize": 8,
    "ytick.labelsize": 8, "figure.dpi": 350,
})

ver = {r["peptide"]: r["status"] for r in csv.DictReader(open("02_Data/cleaned/iedb_epitope_verification.csv", encoding="utf-8"))}
bench = {r["peptide"]: r for r in csv.DictReader(open("02_Data/cleaned/iedb_benchmark_results_v2.csv", encoding="utf-8"))}
pos = [p for p in bench if bench[p]["true_label"] == "POS" and ver[p] == "A2-POS"]
neg = [p for p in bench if bench[p]["true_label"] == "NEG"]

scores, y = [], []
for p in pos:
    scores.append(float(bench[p]["binding_score"])); y.append(1)
for p in neg:
    scores.append(float(bench[p]["binding_score"])); y.append(0)
fpr, tpr, _ = roc_curve(y, scores)
a = auc(fpr, tpr)
pred = [1 if bench[p]["pred_binary"] == "1" else 0 for p in (pos + neg)]
cm = confusion_matrix(y, pred)  # [[TN,FP],[FN,TP]]

fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.2))

# (a) ROC curve
ax = axes[0]
ax.plot(fpr, tpr, color="#4C72B0", lw=1.6, label=f"AUC = {a:.3f}")
ax.plot([0, 1], [0, 1], ls="--", lw=0.8, color="gray")
ax.set_xlabel("False positive rate")
ax.set_ylabel("True positive rate")
ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
ax.legend(frameon=False, loc="lower right")
ax.set_title("(a)", loc="left", fontsize=9)

# (b) confusion matrix
ax = axes[1]
im = ax.imshow(cm, cmap="Blues", vmin=0, vmax=cm.max(), origin="lower")
labels = [["TN\n%d" % cm[0, 0], "FP\n%d" % cm[0, 1]],
          ["FN\n%d" % cm[1, 0], "TP\n%d" % cm[1, 1]]]
ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
ax.set_xticklabels(["Predicted −", "Predicted +"])
ax.set_yticklabels(["Actual −", "Actual +"])
for i in range(2):
    for j in range(2):
        ax.text(j, i, labels[i][j], ha="center", va="center",
                color="white" if cm[i, j] > cm.max() / 2 else "black", fontsize=8)
ax.set_title("(b)", loc="left", fontsize=9)

for ax in axes:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig.tight_layout()
for out in ["06_Bioinformatics_Submission/figures/Figure3.pdf",
            "06_Bioinformatics_Submission/latex_template/figures/Figure3.pdf"]:
    fig.savefig(out)
    print("saved", out)
