#!/usr/bin/env python3
"""
Generate MHCflurry binding predictions for random 9-mer peptides.
Labels: SB (%rank < 0.5), WB (0.5-2.0), NB (>2.0) — netMHCpan convention.

Usage:
    python label_peptides_mhcflurry.py --n 100000 --allele HLA-A*02:01 --out data/real_peptides.csv
"""

import argparse, csv, random, sys
import numpy as np

AA = list("ARNDCQEGHILKMFPSTWYV")

def gen_peptides(n, length=9):
    return ["".join(random.choices(AA, k=length)) for _ in range(n)]

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=100000)
    p.add_argument("--allele", type=str, default="HLA-A*02:01")
    p.add_argument("--out", type=str, default="data/real_peptides.csv")
    p.add_argument("--batch", type=int, default=10000)
    args = p.parse_args()

    random.seed(42)
    np.random.seed(42)

    print(f"Generating {args.n:,} random 9-mer peptides...")
    peptides = gen_peptides(args.n)

    print("Loading MHCflurry...")
    import warnings; warnings.filterwarnings('ignore')
    from mhcflurry import Class1AffinityPredictor
    predictor = Class1AffinityPredictor.load()

    print(f"Predicting {args.allele} binding ({args.n:,} peptides, batch={args.batch})...")
    total_batches = (len(peptides) + args.batch - 1) // args.batch
    all_pct = []
    all_aff = []

    for b in range(total_batches):
        start = b * args.batch
        batch = peptides[start:start + args.batch]
        df = predictor.predict_to_dataframe(
            peptides=batch,
            allele=args.allele,
            include_percentile_ranks=True
        )
        all_pct.extend(df["prediction_percentile"].values)
        all_aff.extend(df["prediction"].values)

        if (b + 1) % 5 == 0 or b == total_batches - 1:
            pct_done = min(100, (start + len(batch)) / len(peptides) * 100)
            print(f"  {start + len(batch):,}/{args.n:,} ({pct_done:.0f}%)")

    all_pct = np.array(all_pct)
    all_aff = np.array(all_aff)

    # netMHCpan thresholds
    labels = np.full(len(peptides), "NB", dtype=object)
    labels[all_pct < 2.0] = "WB"
    labels[all_pct < 0.5] = "SB"

    label_nums = np.zeros(len(peptides), dtype=int)
    label_nums[labels == "WB"] = 1
    label_nums[labels == "SB"] = 2

    n_sb = int(np.sum(label_nums == 2))
    n_wb = int(np.sum(label_nums == 1))
    n_nb = int(np.sum(label_nums == 0))
    print(f"Distribution: SB={n_sb:,} ({100*n_sb/args.n:.1f}%) "
          f"WB={n_wb:,} ({100*n_wb/args.n:.1f}%) "
          f"NB={n_nb:,} ({100*n_nb/args.n:.1f}%)")

    # Balance
    min_count = min(n_sb, n_wb, n_nb)
    print(f"Balancing to {min_count:,} per class ({min_count*3:,} total)...")

    sb_idx = np.where(label_nums == 2)[0]
    wb_idx = np.where(label_nums == 1)[0]
    nb_idx = np.where(label_nums == 0)[0]

    sel_sb = np.sort(np.random.choice(sb_idx, min_count, replace=False))
    sel_wb = np.sort(np.random.choice(wb_idx, min_count, replace=False))
    sel_nb = np.sort(np.random.choice(nb_idx, min_count, replace=False))
    selected = np.sort(np.concatenate([sel_sb, sel_wb, sel_nb]))

    rows = []
    for new_id, idx in enumerate(selected):
        rows.append({
            "peptide_id":  f"PEP_{new_id+1:05d}",
            "peptide":     peptides[idx],
            "label_chr":   labels[idx],
            "label_num":   str(label_nums[idx]),
            "pct_rank":    f"{all_pct[idx]:.4f}",
            "affinity_nM": f"{all_aff[idx]:.1f}",
            "hla_allele":  args.allele
        })

    random.shuffle(rows)
    for i, row in enumerate(rows):
        row["peptide_id"] = f"PEP_{i+1:05d}"

    random.seed(42)
    test_n = int(len(rows) * 0.1)
    test_set = set(random.sample(range(len(rows)), test_n))
    for i, row in enumerate(rows):
        row["data_type"] = "test" if i in test_set else "train"

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "peptide_id", "peptide", "label_chr", "label_num",
            "pct_rank", "affinity_nM", "hla_allele", "data_type"
        ])
        w.writeheader()
        w.writerows(rows)

    n_train = sum(1 for r in rows if r["data_type"] == "train")
    n_test = len(rows) - n_train
    print(f"Saved: {args.out} ({len(rows):,} rows)")
    print(f"Train: {n_train:,} | Test: {n_test:,}")

if __name__ == "__main__":
    main()
