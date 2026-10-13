"""
Run NetMHCpan-4.1 via IEDB API on 519 test peptides and score against MHCflurry labels.

Strategy:
- Batch peptides in FASTA format (50 per API call → ~11 total requests)
- Match results back by seq_num (IEDB uses 1-based index in the FASTA order)
- Classify by %rank: <0.5 => SB(2), 0.5-2.0 => WB(1), >=2.0 => NB(0)
"""
import pandas as pd
import requests
import time
import io
import os
from sklearn.metrics import confusion_matrix, accuracy_score

# ─── 1. Load data ────────────────────────────────────────────────────────────
DATA_PATH = r"D:\Researching\Peptide epitope\02_Data\raw\real_peptides.csv"
OUT_PATH  = r"D:\Researching\Peptide epitope\02_Data\cleaned\netmhcpan_results.csv"

df = pd.read_csv(DATA_PATH)
test_df = df[df["data_type"] == "test"].copy().reset_index(drop=True)
print(f"Test rows: {len(test_df)}")
print("Label distribution:")
print(test_df["label_num"].value_counts().sort_index())
print()

# ─── 2. IEDB API helper (FASTA batch) ────────────────────────────────────────
API_URL  = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
ALLELE   = "HLA-A*02:01"
LENGTH   = 9
BATCH_SZ = 50   # peptides per API call
TIMEOUT  = 300  # seconds per request (server can be slow)

def query_iedb_batch(peptides: list[str], allele: str = ALLELE,
                     length: int = LENGTH, retries: int = 5) -> dict[int, float]:
    """
    Query IEDB for a batch of peptides via FASTA.
    Returns {1-based_seq_num: percentile_rank} dict.
    """
    fasta_lines = []
    for i, pep in enumerate(peptides, 1):
        fasta_lines.append(f">pep{i}")
        fasta_lines.append(pep)
    fasta_text = "\n".join(fasta_lines)

    data = {
        "method":        "netmhcpan_ba",
        "sequence_text": fasta_text,
        "allele":        allele,
        "length":        str(length),
    }

    for attempt in range(retries):
        try:
            resp = requests.post(API_URL, data=data, timeout=TIMEOUT)
            resp.raise_for_status()
            text = resp.text.strip()
            # Parse TSV; header starts with 'allele'
            result_dict = {}
            for line in text.splitlines():
                if not line.strip() or line.startswith("allele"):
                    continue
                parts = line.split("\t")
                # cols: allele, seq_num, start, end, length, peptide, core, icore, ic50, percentile_rank
                seq_num      = int(parts[1])
                pct_rank     = float(parts[-1])
                result_dict[seq_num] = pct_rank
            return result_dict
        except Exception as e:
            wait = 2 ** attempt
            print(f"  [attempt {attempt+1}] Error: {e}. Retrying in {wait}s …")
            time.sleep(wait)
    raise RuntimeError(f"Failed after {retries} retries for batch starting with {peptides[0]!r}")

# ─── 3. Run batched API calls ─────────────────────────────────────────────────
peptides = test_df["peptide"].tolist()
ranks = [float("nan")] * len(peptides)

total_batches = (len(peptides) + BATCH_SZ - 1) // BATCH_SZ
print(f"Running {len(peptides)} peptides in {total_batches} batches of {BATCH_SZ} …\n")

for batch_idx in range(total_batches):
    start = batch_idx * BATCH_SZ
    end   = min(start + BATCH_SZ, len(peptides))
    batch = peptides[start:end]

    print(f"  Batch {batch_idx+1}/{total_batches} (peptides {start+1}–{end}) …", end=" ", flush=True)
    t0 = time.time()
    try:
        result = query_iedb_batch(batch)
        # result keys are 1-based within the batch
        for local_seqnum, rank in result.items():
            global_idx = start + local_seqnum - 1  # convert to 0-based global index
            if 0 <= global_idx < len(peptides):
                ranks[global_idx] = rank
        ok = sum(1 for k in result if 0 <= start + k - 1 < len(peptides))
        print(f"done ({ok} results, {time.time()-t0:.1f}s)")
    except RuntimeError as e:
        print(f"\n  FAILED: {e}")

    # Polite pause between batches
    if batch_idx < total_batches - 1:
        time.sleep(2)

# ─── 4. Classify & evaluate ──────────────────────────────────────────────────
test_df = test_df.copy()
test_df["netmhcpan_rank"] = ranks

def classify_rank(r):
    if pd.isna(r):
        return -1    # failed
    if r < 0.5:
        return 2     # SB
    elif r < 2.0:
        return 1     # WB
    else:
        return 0     # NB

test_df["netmhcpan_class"] = test_df["netmhcpan_rank"].apply(classify_rank)

failed = (test_df["netmhcpan_class"] == -1).sum()
if failed > 0:
    print(f"\n⚠ {failed} peptides failed — excluded from accuracy calculation.")

valid = test_df[test_df["netmhcpan_class"] >= 0]
y_true = valid["label_num"].values
y_pred = valid["netmhcpan_class"].values

acc = accuracy_score(y_true, y_pred)
cm  = confusion_matrix(y_true, y_pred, labels=[0, 1, 2])

print("\n" + "="*58)
print(f"  Accuracy (NetMHCpan vs label_num): {acc:.4f}  ({acc*100:.2f}%)")
print(f"  Valid peptides evaluated:          {len(valid)} / {len(test_df)}")
print()
print("  Confusion matrix (rows=true label, cols=predicted):")
print("             pred_NB  pred_WB  pred_SB")
for label, row_vals in zip(["true_NB(0)", "true_WB(1)", "true_SB(2)"], cm):
    print(f"  {label}  " + "  ".join(f"{v:6d}" for v in row_vals))
print("="*58)

# Per-class accuracy
for cls, name in [(0,"NB"), (1,"WB"), (2,"SB")]:
    mask = y_true == cls
    if mask.sum() > 0:
        cls_acc = (y_pred[mask] == cls).mean()
        print(f"  Class {name} accuracy: {cls_acc:.4f}  (n={mask.sum()})")

# ─── 5. Save results ─────────────────────────────────────────────────────────
out_dir = os.path.dirname(OUT_PATH)
os.makedirs(out_dir, exist_ok=True)

out_cols = ["peptide", "label_num", "netmhcpan_rank", "netmhcpan_class"]
test_df[out_cols].to_csv(OUT_PATH, index=False)
print(f"\n✓ Results saved to: {OUT_PATH}")
print(f"\n>>> FINAL ACCURACY: {acc:.4f} ({acc*100:.2f}%) <<<")
