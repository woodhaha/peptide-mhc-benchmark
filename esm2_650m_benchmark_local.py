#!/usr/bin/env python3
"""ESM-2 650M (t33) embedding benchmark — local RTX 2000 Ada 8GB (fp16).

Replicates manuscript Methods 2.3 / Table 3 (encoding comparison):
per-position embeddings -> FFN 256->128->64 -> 5-fold stratified CV.
Runs the largest ESM-2 checkpoint locally at fp16 (~1.5GB VRAM), adding the
t33 650M point the small checkpoints (t6/t12) could not cover.
"""
import csv
import os

os.environ.setdefault("HF_ENDPOINT", "https://hf-mirror.com")

import numpy as np
import torch
import torch.nn as nn
from transformers import AutoTokenizer, EsmModel
from sklearn.model_selection import StratifiedKFold

torch.manual_seed(42)
np.random.seed(42)

# 1. data (local real_peptides.csv: peptide, label_chr SB/NB/WB)
peptides, labels = [], []
with open("02_Data/raw/real_peptides.csv") as f:
    for row in csv.DictReader(f):
        peptides.append(row["peptide"])
        labels.append(row["label_chr"])
lab2i = {"NB": 0, "WB": 1, "SB": 2}
y = np.array([lab2i[l] for l in labels])
print(f"n={len(peptides)}  class counts (NB/WB/SB)={np.bincount(y)}", flush=True)

# 2. ESM-2 650M per-position embeddings (9 x 1280 = 11520 dim)
name = "facebook/esm2_t33_650M_UR50D"
tok = AutoTokenizer.from_pretrained(name)
model = EsmModel.from_pretrained(name, torch_dtype=torch.float16).cuda().eval()

X = []
bs = 16
with torch.no_grad():
    for i in range(0, len(peptides), bs):
        b = peptides[i:i + bs]
        t = tok(b, return_tensors="pt", padding=True).to("cuda")
        o = model(**t).last_hidden_state[:, 1:-1, :]  # strip <cls>/<eos>
        X.append(o.reshape(o.shape[0], -1).float().cpu().numpy())
        if (i // bs) % 25 == 0:
            print(f"  embedded {i}/{len(peptides)}", flush=True)
X = np.vstack(X)
print(f"X shape {X.shape}", flush=True)

# 3. FFN 256->128->64, 5-fold CV (matches manuscript Methods 2.3)
def make_ffn(din):
    return nn.Sequential(
        nn.Linear(din, 256), nn.ReLU(),
        nn.Linear(256, 128), nn.ReLU(),
        nn.Linear(128, 64), nn.ReLU(),
        nn.Linear(64, 3),
    )


def train_eval(tr, va):
    net = make_ffn(X.shape[1]).cuda()
    opt = torch.optim.Adam(net.parameters(), lr=0.001)
    lossf = nn.CrossEntropyLoss()
    Xtr = torch.tensor(X[tr], dtype=torch.float32).cuda()
    ytr = torch.tensor(y[tr], dtype=torch.long).cuda()
    Xva = torch.tensor(X[va], dtype=torch.float32).cuda()
    yva = torch.tensor(y[va], dtype=torch.long).cuda()
    for _ in range(150):
        net.train()
        opt.zero_grad()
        loss = lossf(net(Xtr), ytr)
        loss.backward()
        opt.step()
    net.eval()
    with torch.no_grad():
        pred = net(Xva).argmax(1)
    return (pred == yva).float().mean().item()


skf = StratifiedKFold(5, shuffle=True, random_state=42)
accs = [train_eval(tr, va) for tr, va in skf.split(X, y)]
print("ESM-2 650M (t33) 5-fold CV acc:", [f"{a:.4f}" for a in accs], flush=True)
print(f"mean {np.mean(accs):.4f} +/- {np.std(accs):.4f}", flush=True)
