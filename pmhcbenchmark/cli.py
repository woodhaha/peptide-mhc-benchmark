"""Command-line interface for pmhcbenchmark."""
import argparse
import sys

from .predictor import ALLELES, predict


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="pmhc-predict",
        description="Peptide-MHC-I binding prediction (per-allele Deep FFN, BLOSUM62).",
    )
    p.add_argument("peptides", nargs="?", help="comma-separated 9-mer peptides")
    p.add_argument("--file", help="file with one 9-mer per line")
    p.add_argument("--allele", default="HLA-A*02:01", help=f"one of {sorted(ALLELES)}")
    p.add_argument("--threshold", type=float, default=0.5, help="binder cutoff (default 0.5)")
    p.add_argument("--model-dir", default="models", help="directory holding the .keras models")
    args = p.parse_args(argv)

    if args.file:
        peptides = [ln.strip() for ln in open(args.file, encoding="utf-8") if ln.strip()]
    elif args.peptides:
        peptides = [x.strip() for x in args.peptides.split(",") if x.strip()]
    else:
        p.error("provide a peptide, a comma-separated list, or --file")

    for pep in peptides:
        if len(pep) != 9:
            print(f"{pep}\tSKIP (not a 9-mer)")
            continue
        s = float(predict([pep], args.allele, args.model_dir)[0])
        cls = "BINDER" if s >= args.threshold else "NON-BINDER"
        print(f"{pep}\t{args.allele}\tscore={s:.4f}\t{cls}")


if __name__ == "__main__":
    sys.exit(main())
