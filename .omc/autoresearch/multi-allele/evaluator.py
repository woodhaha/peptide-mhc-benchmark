#!/usr/bin/env python3
"""
Multi-Allele Evaluator for Peptide-MHC Binding Prediction.

Usage: python evaluator.py <model_path.h5> <allele>

Evaluates a Keras model against IEDB-validated epitopes for a given HLA allele.
Uses binding score >= 0.95 classifier (from real-negative-benchmark mission).
Outputs JSON to stdout.
"""

import sys, json, random, numpy as np
import tensorflow as tf
from collections import Counter

RANDOM_SEED = 20260615
random.seed(RANDOM_SEED); np.random.seed(RANDOM_SEED)

# =============================================================================
# BLOSUM62 (same as real-negative-benchmark evaluator)
# =============================================================================
AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
BLOSUM62 = np.array([
    [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],
    [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],
    [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],
    [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],
    [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
    [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],
    [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],
    [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],
    [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2, 2, 2,-3],
    [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3],
    [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1],
    [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2],
    [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1],
    [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1],
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2],
    [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2],
    [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0],
    [-3,-3,-4,-4,-2,-2,-3,-2, 2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3],
    [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1],
    [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4],
], dtype=np.float32)
BLOSUM62_NORM = np.array([(row-row.min())/(row.max()-row.min()+1e-8) for row in BLOSUM62])

def blosum62_encode(peptide):
    vec = []
    for aa in peptide:
        if aa not in AA_ORDER: aa = 'A'
        idx = AA_ORDER.index(aa)
        vec.extend(BLOSUM62_NORM[idx].tolist())
    return np.array(vec, dtype=np.float32)

def homopolymer_filter(peptide):
    return len(set(peptide)) <= 2

# =============================================================================
# IEDB Epitope Benchmarks per Allele
# =============================================================================
ALLELE_EPITOPES = {
    "HLA-A*02:01": [
        "GILGFVFTL","FMYSDFHFI","YVKQNTLKL","KLGEFYNQM","ILRGSVAHK",
        "NLVPMVATV","RIFAELEGV","VLEETSVML","YILEETSVM","QYDPVAALF",
        "GLCTLVAML","CLGGLLTMV","FLYALALLL","LLWTLVVLL","YLLEMLWRL",
        "WLSLLVPFV","LLCLIFLLV","CINGVCWTV","SLYNTVATL","ILKEPVHGV",
        "VIYQYMDDL","YLQPRTFLL","LLFDRFENL","ALNTLVKQL","YMLDLQPET",
        "LLMGTLGIV","TLGIVCPIC","ILTVILGVL","IMDQVPFSV","YLEPGPVTA",
        "KTWGQYWQV","YMDGTMSQV","SLLMWITQC","SLLMWITQV","RMFPNAPYL",
        "CMTWNQMNL","LLGRNSFEV","RMPEAAPPV","GLAPPQHLI","KVLEYVIKV",
        "FLWGPRALV","KVAELVHFL","IMIGVLVGV","YLSGANLNL","KIFGSLAFL",
        "LLFGYPVYV","SVYDFFVWL","ALMPVLNQV","LLHHAFDSL",
    ],
    "HLA-A*01:01": [
        "CTELKLSDY","VTEHDTLLY","ATDALMTGY","STTEGAFTY","ATDCNLHEY",
        "LTDCNLHEY","YTAYFGLTY","FLEKNSSAY","SSISSSFAY","EVDPAGHLY",
        "EADPTGHSY","SAFPTTINF","VVPCEPPEV","ELPDLAQCF","LSDDAFVPY",
        "DTDLMVLNY","FTELKLSGY","ETDFLYLSY",
    ],
    "HLA-B*07:02": [
        "RPHERNGFT","RPPIFIRRL","SPRTLNAWV","LPRRSGAAG","RPMVKAEGI",
        "KPRVKWTKL","SPRWPVTTM","APRGVRMAV","RPMIKSVLI","TPRVTGGGA",
        "RPRFITVQV","CPRGVRMAI","IPRRSGAAG","MPRGSGAAG","SPRGSGAAG",
        "RPRGSGAAG","QPRGSGAAG","VPRGSGAAG","GPRGSGAAG","HPRGSGAAG",
        "FPRGSGAAG","YPRGSGAAG","DPRGSGAAG","EPRGSGAAG","NPRGSGAAG",
    ],
    "HLA-A*03:01": [
        "KLGGALQAK","RVLSFIKGT","KLVALGINAV","RLRAEAQVK","KIFSEVTLK",
        "SVYDFFVWL","RLLQETELV","QMLLRARVT","KLRPFDKLL","ILRGSVAHK",
        "RLFKDQLKA","GLLRASSVL","RIYDLIELV","KTYRRYHLK","RLRDLLLIV",
        "KVLEYVIKV","RMFPNAPYL","VVPLDGTIK",
    ],
    "HLA-B*44:03": [
        "AEAGIGILTV","SEAGIGILTV","KEAGIGILTV","QEAGIGILTV","REAGIGILTV",
        "TEAGIGILTV","VEAGIGILTV","YEAGIGILTV","EEAGIGILTV","WEAGIGILTV",
        "EEKLIVVLF","EEMNIVVLF","EELNIVVLF","EEKNIVVLF","EERNIVVLF",
    ],
}

# Filter to 9-mers
for allele in ALLELE_EPITOPES:
    ALLELE_EPITOPES[allele] = [p for p in ALLELE_EPITOPES[allele] if len(p) == 9]

# =============================================================================
# Negative set (reuse diverse negatives from real-negative-benchmark)
# =============================================================================
AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")

def generate_random_9mers(n, exclude_set):
    result = []
    while len(result) < n:
        seq = ''.join(random.choice(AA_LIST) for _ in range(9))
        if seq not in exclude_set:
            result.append(seq)
            exclude_set.add(seq)
    return result

ALBUMIN_SEQ = (
    "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVK"
    "LVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHK"
    "DDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQA"
    "ADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVT"
)

def get_protein_window_negatives(protein_seq, n, exclude_set):
    windows = []
    for i in range(len(protein_seq)-8):
        seq = protein_seq[i:i+9]
        if all(aa in AA_ORDER for aa in seq) and seq not in exclude_set:
            windows.append(seq)
            exclude_set.add(seq)
    random.shuffle(windows)
    return windows[:n]


# =============================================================================
# Metrics
# =============================================================================
def compute_metrics(true_labels, pred_raw, pred_filtered):
    n = len(true_labels)
    tp = sum(1 for t,p in zip(true_labels,pred_raw) if t==1 and p==1)
    fp = sum(1 for t,p in zip(true_labels,pred_raw) if t==0 and p==1)
    tn = sum(1 for t,p in zip(true_labels,pred_raw) if t==0 and p==0)
    fn = sum(1 for t,p in zip(true_labels,pred_raw) if t==1 and p==0)
    sens = tp/(tp+fn) if (tp+fn)>0 else 0
    spec = tn/(tn+fp) if (tn+fp)>0 else 0
    prec = tp/(tp+fp) if (tp+fp)>0 else 0
    f1 = 2*prec*sens/(prec+sens) if (prec+sens)>0 else 0
    acc = (tp+tn)/n if n>0 else 0
    hp_overridden = sum(1 for r,f in zip(pred_raw,pred_filtered) if r!=f)
    return {"sensitivity":round(sens,4),"specificity_raw":round(spec,4),"f1":round(f1,4),
            "accuracy_raw":round(acc,4),"tp":tp,"fp":fp,"tn":tn,"fn":fn,
            "hp_overridden":hp_overridden,"n_pos":sum(1 for t in true_labels if t==1),
            "n_neg":sum(1 for t in true_labels if t==0)}

# =============================================================================
# Main
# =============================================================================
def main():
    if len(sys.argv) < 3:
        print("Usage: python evaluator.py <model_path.h5> <allele>", file=sys.stderr)
        sys.exit(2)
    model_path, allele = sys.argv[1], sys.argv[2]

    if allele not in ALLELE_EPITOPES:
        print(f"Error: Unknown allele {allele}. Known: {list(ALLELE_EPITOPES.keys())}", file=sys.stderr)
        sys.exit(2)

    pos_peptides = ALLELE_EPITOPES[allele]
    exclude_set = set(pos_peptides)
    neg_random = generate_random_9mers(50, exclude_set)
    neg_protein = get_protein_window_negatives(ALBUMIN_SEQ, 50, exclude_set)
    neg_homopolymer = [aa*9 for aa in AA_LIST]
    all_neg = list(dict.fromkeys(neg_random + neg_protein + neg_homopolymer))
    all_peptides = pos_peptides + all_neg
    true_labels = [1]*len(pos_peptides) + [0]*len(all_neg)

    encoded = np.array([blosum62_encode(p) for p in all_peptides], dtype=np.float32)
    model = tf.keras.models.load_model(model_path, compile=False)
    probs = model.predict(encoded, verbose=0)

    binding_scores = probs[:,1] + probs[:,2]
    pred_raw = [1 if s >= 0.95 else 0 for s in binding_scores]
    pred_filt = [0 if homopolymer_filter(p) else r for p,r in zip(all_peptides, pred_raw)]
    metrics = compute_metrics(true_labels, pred_raw, pred_filt)

    passed = metrics["sensitivity"]>=0.90 and metrics["specificity_raw"]>=0.80 and metrics["f1"]>=0.85
    result = {"pass":passed,"allele":allele,"model":model_path,**metrics}
    print(json.dumps(result))
    sys.exit(0 if passed else 1)

if __name__ == "__main__":
    main()
