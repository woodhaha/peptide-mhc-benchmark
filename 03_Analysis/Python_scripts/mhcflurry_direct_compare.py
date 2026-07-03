"""
Compare MHCflurry 2.2.x direct predictions vs trained models
on the held-out test set.

Thresholds (from MHCflurry docs):
  SB = %rank < 0.5,  WB = 0.5 <= %rank < 2.0,  NB = %rank >= 2.0
"""

import csv, time
import numpy as np
from collections import Counter
from sklearn.metrics import accuracy_score, f1_score, classification_report

# --- 1. Load test peptides ---
test_data = []
with open(r'D:\Researching\Peptide epitope\02_Data\cleaned\feature_resnet.csv', 'r') as f:
    for row in csv.DictReader(f):
        if row['data_type'] == 'test':
            test_data.append({
                'peptide': row['peptide'],
                'label': row['label_chr'],
                'model_pred': int(row['pred_class']),
            })

test_peptides = [d['peptide'] for d in test_data]
label_map = {0: 'NB', 1: 'WB', 2: 'SB'}
model_pred_chr = np.array([label_map[d['model_pred']] for d in test_data])
y_true = np.array([d['label'] for d in test_data])

print(f"Test set: {len(test_data)} peptides")
print(f"Label distribution: {dict(Counter(y_true))}")
print(f"ResNet pred distribution: {dict(Counter(model_pred_chr))}")

# --- 2. Run MHCflurry ---
print("\nLoading MHCflurry predictor...")
import mhcflurry
predictor = mhcflurry.Class1AffinityPredictor.load()

print("Running predictions...")
t0 = time.time()
df = predictor.predict_to_dataframe(peptides=test_peptides, allele='HLA-A*02:01')
elapsed = time.time() - t0
print(f"Done in {elapsed:.1f}s")

# df columns: peptide, allele, prediction(=affinity_nM), prediction_low, prediction_high, prediction_percentile(=rank%)
ranks = df['prediction_percentile'].values
affinities = df['prediction'].values

print(f"Affinity range: {affinities.min():.1f} - {affinities.max():.1f} nM")
print(f"Rank range: {ranks.min():.4f} - {ranks.max():.4f}")

# --- 3. Classify ---
def classify(rank):
    if rank < 0.5: return 'SB'
    elif rank < 2.0: return 'WB'
    else: return 'NB'

mhcflurry_pred = np.array([classify(r) for r in ranks])
print(f"MHCflurry direct pred distribution: {dict(Counter(mhcflurry_pred))}")

# --- 4. Metrics ---
mhcflurry_acc  = accuracy_score(y_true, mhcflurry_pred)
mhcflurry_f1   = f1_score(y_true, mhcflurry_pred, average='macro', labels=['NB','WB','SB'])
mhcflurry_f1pc = f1_score(y_true, mhcflurry_pred, average=None, labels=['NB','WB','SB'])

model_acc  = accuracy_score(y_true, model_pred_chr)
model_f1   = f1_score(y_true, model_pred_chr, average='macro', labels=['NB','WB','SB'])
model_f1pc = f1_score(y_true, model_pred_chr, average=None, labels=['NB','WB','SB'])

print(f"\n{'='*55}")
print(f"{'Metric':<18} {'ResNet':>8} {'MHCflurry 2.2.1':>12} {'Delta':>8}")
print(f"{'='*55}")
for name, m, f, idx in [('Accuracy', model_acc, mhcflurry_acc, None),
                          ('Macro F1', model_f1, mhcflurry_f1, None),
                          ('NB F1', model_f1pc[0], mhcflurry_f1pc[0], 0),
                          ('WB F1', model_f1pc[1], mhcflurry_f1pc[1], 1),
                          ('SB F1', model_f1pc[2], mhcflurry_f1pc[2], 2)]:
    delta = f - m if idx is None else mhcflurry_f1pc[idx] - model_f1pc[idx]
    print(f"{name:<18} {m:8.4f} {f:12.4f} {delta:+8.4f}")
print(f"{'='*55}")

print("\n--- MHCflurry direct Classification Report ---")
print(classification_report(y_true, mhcflurry_pred, digits=4))

# --- 5. Affinity-based (50/500nM) for reference ---
aff_pred = np.where(affinities < 50, 'SB', np.where(affinities < 500, 'WB', 'NB'))
print("\n--- Affinity-based (50/500nM) ---")
print(classification_report(y_true, aff_pred, digits=4))

# --- 6. Save CSV ---
rows = []
for i, d in enumerate(test_data):
    rows.append({
        'peptide': d['peptide'], 'true_label': d['label'],
        'resnet_pred': model_pred_chr[i], 'mhcflurry_pred': mhcflurry_pred[i],
        'mhcflurry_rank': round(float(ranks[i]), 6),
        'mhcflurry_affinity_nM': round(float(affinities[i]), 2),
        'resnet_correct': model_pred_chr[i] == y_true[i],
        'mhcflurry_correct': mhcflurry_pred[i] == y_true[i],
    })

def fmt(v):
    return f"{v:.4f}" if isinstance(v, (int, float)) else str(v)

summary = [
    {'peptide': k, 'true_label': '', 'resnet_pred': fmt(v_r),
     'mhcflurry_pred': fmt(v_m), 'mhcflurry_rank': '', 'mhcflurry_affinity_nM': '',
     'resnet_correct': '', 'mhcflurry_correct': ''}
    for k, v_r, v_m in [
        ('===SUMMARY===', '', ''),
        ('Accuracy', model_acc, mhcflurry_acc),
        ('Macro_F1', model_f1, mhcflurry_f1),
        ('NB_F1', model_f1pc[0], mhcflurry_f1pc[0]),
        ('WB_F1', model_f1pc[1], mhcflurry_f1pc[1]),
        ('SB_F1', model_f1pc[2], mhcflurry_f1pc[2]),
    ]
]

path = r'D:\Researching\Peptide epitope\02_Data\cleaned\mhcflurry_direct_comparison.csv'
with open(path, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys())
    w.writeheader(); w.writerows(rows); w.writerows(summary)
print(f"\nSaved {path}")

# --- 7. Disagreements ---
disagree = mhcflurry_pred != model_pred_chr
print(f"\nResNet vs MHCflurry disagreements: {disagree.sum()}/{len(test_data)} ({disagree.sum()/len(test_data)*100:.1f}%)")
for i in np.where(disagree)[0]:
    d = test_data[i]
    print(f"  {d['peptide']}: true={d['label']}  ResNet={model_pred_chr[i]}  "
          f"MHCflurry={mhcflurry_pred[i]}  (rank={ranks[i]:.3f}, aff={affinities[i]:.0f}nM)")
