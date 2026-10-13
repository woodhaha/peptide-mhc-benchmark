#!/usr/bin/env python3
"""NetMHCpan-4.1 SOTA comparison on the held-out test set.

Runs NetMHCpan-4.1 (academic web server) on the 519 held-out test peptides,
classifies them SB/WB/NB with the same %rank thresholds used for MHCflurry
labels, and reports accuracy vs the MHCflurry-derived labels.

This closes the "no SOTA comparison" gap flagged by all three reviewers.

Requirements: network access to services.healthtech.dtu.dk (use Clash proxy
127.0.0.1:15715 if direct fails). NetMHCpan-4.1 is free for academic use via
the web server; the local software requires a DTU license.

Usage:
    python src/netmhcpan_compare.py --fasta            # step 1: emit FASTA
    python src/netmhcpan_compare.py --parse results.xls  # step 2: parse + score
"""
import argparse, csv, io, json, os, time, urllib.request, urllib.parse
import numpy as np
import pandas as pd

ALLELE = "HLA-A*02:01"
SB_PCT = 0.5   # %rank < 0.5  -> SB
WB_PCT = 2.0   # 0.5 <= %rank < 2.0 -> WB ; else NB

# NetMHCpan-4.1 web server (DTU Health Tech)
NETMHCPAN_URL = "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi"


def load_test_peptides():
    df = pd.read_csv("02_Data/raw/real_peptides.csv")
    test = df[df["data_type"] == "test"]
    return test[["peptide", "label_num"]].values


def emit_fasta(peptides, out="data/netmhcpan_test.fasta"):
    with open(out, "w") as f:
        for i, (pep, _) in enumerate(peptides, 1):
            f.write(f">pep_{i}\n{pep}\n")
    print(f"Wrote {len(peptides)} peptides to {out}")
    print(f"Submit FASTA to NetMHCpan-4.1 web server for allele {ALLELE}, length 9,")
    print(f"then download the XLS/TSV result and run: python src/netmhcpan_compare.py --parse <file>")


def submit_web(peptides):
    """Best-effort programmatic submission to the NetMHCpan web server.
    The DTU webface2 API is version-fragile; if this fails, use emit_fasta + manual
    web submission. Proxies are respected via HTTPS_PROXY env if set."""
    seq = "\n".join(f">p{i}\n{p}" for i, (p, _) in enumerate(peptides))
    params = {
        "configfile": json.dumps({
            "version": "4.1", "method": "netmhcpan", "allele": ALLELE,
            "length": "9", "sort": "-1", "thrs": "0.5,2.0", "peplist": seq,
        }),
        "SEQPASTE": seq, "len": "9", "allele": ALLELE,
    }
    data = urllib.parse.urlencode(params).encode()
    req = urllib.request.Request(NETMHCPAN_URL, data=data,
                                 headers={"User-Agent": "Mozilla/5.0"})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            html = r.read().decode("utf-8", "ignore")
        # job id is usually in a redirect/URL; crude extraction
        m = __import__("re").search(r"jobid=([A-Za-z0-9]+)", html)
        return m.group(1) if m else None
    except Exception as e:
        print(f"Web submission failed: {e}")
        return None


def parse_and_score(rows):
    """Parse NetMHCpan output rows [(peptide, pct_rank)] and score vs MHCflurry labels."""
    test = load_test_peptides()
    lab = {p: int(l) for p, l in test}
    # NetMHCpan -> SB/WB/NB via %rank
    pred = {}
    for pep, pr in rows:
        if pr < SB_PCT: pred[pep] = 2      # SB
        elif pr < WB_PCT: pred[pep] = 1    # WB
        else: pred[pep] = 0                # NB
    ys, yp = [], []
    for p, l in test:
        if p in pred:
            ys.append(l); yp.append(pred[p])
    acc = sum(a == b for a, b in zip(ys, yp)) / len(ys)
    print(f"NetMHCpan-4.1 vs MHCflurry labels: n={len(ys)}, accuracy={acc*100:.1f}%")
    return acc


def parse_xls(path):
    """Parse a NetMHCpan-4.1 XLS/TSV download into [(peptide, pct_rank)].
    NetMHCpan output columns include 'Pos','MHC','Peptide',...,'%Rank'. This
    handles both the xls (tab-separated) and a plain TSV."""
    rows = []
    with open(path, encoding="utf-8", errors="ignore") as f:
        lines = f.read().splitlines()
    # find header row containing '%Rank'
    hdr = next(i for i, l in enumerate(lines) if "%Rank" in l or "%rank" in l)
    cols = lines[hdr].split("\t")
    try:
        pep_i = next(i for i, c in enumerate(cols) if "Peptide" in c)
        rank_i = next(i for i, c in enumerate(cols) if "Rank" in c)
    except StopIteration:
        # fallback to known NetMHCpan-4.1 column order
        pep_i, rank_i = 2, -1
    for l in lines[hdr + 1:]:
        c = l.split("\t")
        if len(c) <= max(pep_i, rank_i): continue
        try:
            rows.append((c[pep_i].strip(), float(c[rank_i])))
        except ValueError:
            continue
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", action="store_true")
    ap.add_argument("--parse")
    args = ap.parse_args()
    test = load_test_peptides()
    if args.fasta:
        emit_fasta(test)
    elif args.parse:
        parse_and_score(parse_xls(args.parse))
    else:
        # try programmatic submit; fall back to FASTA instructions
        job = submit_web(test)
        if job:
            print(f"Submitted, job id {job} — poll results then --parse")
        else:
            emit_fasta(test)


if __name__ == "__main__":
    main()
