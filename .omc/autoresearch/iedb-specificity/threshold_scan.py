import csv

results = []
with open('Data/cleaned/iedb_benchmark_results.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        results.append(row)

print('=== FALSE POSITIVES (homopolymers) ===')
for r in results:
    if r['true_label'] == 'NEG' and r['pred_class'] in ('SB','WB'):
        print(f"  {r['peptide']:12s} pred={r['pred_class']:3s} SB={float(r['prob_SB']):.4f} WB={float(r['prob_WB']):.4f} NB={float(r['prob_NB']):.4f}")

print('\n=== FALSE NEGATIVES ===')
for r in results:
    if r['true_label'] == 'POS' and r['pred_class'] == 'NB':
        print(f"  {r['peptide']:12s} src={r['source_protein']:20s} SB={float(r['prob_SB']):.4f} WB={float(r['prob_WB']):.4f} NB={float(r['prob_NB']):.4f} p9={r['p9']}")

print('\n=== Threshold Grid Search (SB >= threshold; WB >= 0.3) ===')
for th in [0.50, 0.60, 0.70, 0.75, 0.80, 0.85, 0.90, 0.92, 0.95, 0.97, 0.99]:
    tp = fp = tn = fn = 0
    for r in results:
        prob_sb = float(r['prob_SB'])
        prob_wb = float(r['prob_WB'])
        new_class = 'SB' if prob_sb >= th else ('WB' if prob_wb >= 0.3 else 'NB')
        is_binder = new_class in ('SB','WB')
        if r['true_label'] == 'POS' and is_binder: tp += 1
        if r['true_label'] == 'NEG' and is_binder: fp += 1
        if r['true_label'] == 'NEG' and not is_binder: tn += 1
        if r['true_label'] == 'POS' and not is_binder: fn += 1

    sens = tp/(tp+fn) if (tp+fn)>0 else 0
    spec = tn/(tn+fp) if (tn+fp)>0 else 0
    pmark = '  ** PASS **' if spec >= 0.85 and sens >= 0.90 else ''
    print(f'  SB>={th:.2f}: Sens={sens:.4f} ({tp}/{tp+fn})  Spec={spec:.4f} ({tn}/{tn+fp})  FP={fp}  FN={fn}{pmark}')
