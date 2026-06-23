#!/usr/bin/env python3
"""
Label peptides using IEDB MHC-I binding prediction API.
Zero external dependencies — uses only Python stdlib.
Produces CSV compatible with peptide_mhc_binding_study.R.

Usage:
    python label_peptides_iedb.py --n 50000 --allele HLA-A*02:01 --out data/real_peptides.csv
"""

import argparse, csv, json, random, sys, time, urllib.request, urllib.parse

IEDB_URL = "https://tools-cluster-interface.iedb.org/api/v3/mhc_i"
AA = list("ARNDCQEGHILKMFPSTWYV")

def gen_peptides(n, length=9):
    return ["".join(random.choices(AA, k=length)) for _ in range(n)]

def query_iedb(peptides, allele):
    data = urllib.parse.urlencode({
        "sequence_text": "\n".join(peptides),
        "method": "consensus",
        "allele": allele,
        "length": "9"
    }).encode()
    req = urllib.request.Request(IEDB_URL, data,
        headers={"Content-Type": "application/x-www-form-urlencoded"})
    try:
        with urllib.request.urlopen(req, timeout=300) as r:
            res = json.loads(r.read())
        ranks = [float(e.get("percentile_rank", 100)) for e in res.get("results", [])]
        affs  = [float(e.get("affinity", 99999)) for e in res.get("results", [])]
        return ranks, affs
    except Exception as e:
        print(f"  API error: {e}", file=sys.stderr)
        return None, None

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=50000)
    p.add_argument("--allele", type=str, default="HLA-A*02:01")
    p.add_argument("--out", type=str, default="data/real_peptides.csv")
    p.add_argument("--batch", type=int, default=100)
    args = p.parse_args()

    random.seed(42)

    print(f"Generating {args.n:,} random 9-mer peptides...")
    peptides = gen_peptides(args.n)

    print(f"Querying IEDB API ({args.batch}/batch)...")
    ranks, affs = [], []
    n_batches = (len(peptides) + args.batch - 1) // args.batch

    for b in range(n_batches):
        start = b * args.batch
        batch = peptides[start:start + args.batch]
        r, a = query_iedb(batch, args.allele)
        if r is None:
            print(f"  Batch {b+1}/{n_batches} failed, retrying...")
            time.sleep(5)
            r, a = query_iedb(batch, args.allele)
            if r is None:
                r = [100.0] * len(batch)
                a = [50000.0] * len(batch)
        ranks.extend(r)
        affs.extend(a)
        if (b + 1) % 10 == 0 or b == n_batches - 1:
            print(f"  {start + len(batch):,}/{args.n:,} ({(b+1)/n_batches*100:.0f}%)")

    # Label using netMHCpan %rank thresholds
    labels = ["NB"] * len(peptides)
    for i, r in enumerate(ranks):
        if r < 0.5:
            labels[i] = "SB"
        elif r < 2.0:
            labels[i] = "WB"

    n_sb = labels.count("SB")
    n_wb = labels.count("WB")
    n_nb = labels.count("NB")
    print(f"Distribution: SB={n_sb:,} ({100*n_sb/args.n:.1f}%) "
          f"WB={n_wb:,} ({100*n_wb/args.n:.1f}%) "
          f"NB={n_nb:,} ({100*n_nb/args.n:.1f}%)")

    # Balance classes
    min_count = min(n_sb, n_wb, n_nb)
    print(f"Balancing to {min_count:,} per class ({min_count*3:,} total)...")

    sb_idx = [i for i, l in enumerate(labels) if l == "SB"]
    wb_idx = [i for i, l in enumerate(labels) if l == "WB"]
    nb_idx = [i for i, l in enumerate(labels) if l == "NB"]

    sel_sb = sorted(random.sample(sb_idx, min_count))
    sel_wb = sorted(random.sample(wb_idx, min_count))
    sel_nb = sorted(random.sample(nb_idx, min_count))
    selected = sorted(sel_sb + sel_wb + sel_nb)

    # Build rows
    rows = []
    for new_id, idx in enumerate(selected, 1):
        rows.append({
            "peptide_id":  f"PEP_{new_id:05d}",
            "peptide":     peptides[idx],
            "label_chr":   labels[idx],
            "label_num":   {"NB": 0, "WB": 1, "SB": 2}[labels[idx]],
            "pct_rank":    f"{ranks[idx]:.4f}",
            "affinity_nM": f"{affs[idx]:.1f}",
            "hla_allele":  args.allele
        })

    # Shuffle
    random.shuffle(rows)
    for i, row in enumerate(rows, 1):
        row["peptide_id"] = f"PEP_{i:05d}"

    # Train/test split
    random.seed(42)
    test_n = int(len(rows) * 0.1)
    test_set = set(random.sample(range(len(rows)), test_n))
    for i, row in enumerate(rows):
        row["data_type"] = "test" if i in test_set else "train"

    # Write CSV
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "peptide_id", "peptide", "label_chr", "label_num",
            "pct_rank", "affinity_nM", "hla_allele", "data_type"
        ])
        w.writeheader()
        w.writerows(rows)

    print(f"Saved: {args.out} ({len(rows):,} rows)")
    print(f"Train: {len(rows) - test_n:,} | Test: {test_n:,}")

if __name__ == "__main__":
    main()
