# Peptide Binding Prediction Results — Terminal Output

> Model: Deep FFN (MHCflurry-trained, 91.9% accuracy)
> Allele: HLA-A*02:01
> Date: June 15, 2026

---

```
==============================================================
  Peptide–MHC Binding Prediction Results
  Model: Deep FFN (MHCflurry-trained, 91.9% accuracy)
  Allele: HLA-A*02:01
==============================================================

   Peptide     Anchor  Predicted  Confidence  prob_SB  prob_WB  prob_NB
   LMAFYLYEV   M-p9-V  SB         100.0       100.00   0.00     0.00
   LLTDAQRIV   L-p9-V  SB          72.1        72.10  27.80     0.03
   FQTDRISYA   Q-p9-A  SB          94.6        94.60   5.37     0.00
   SLHLTNCFV   L-p9-V  SB          97.6        97.60   2.42     0.00
   AAAAAAAAA   A-p9-A  NB          95.6         0.00   4.42    95.60
   LLLLLLLLL   L-p9-L  SB          98.9        98.90   1.11     0.00
   KGWGHSNGS   G-p9-S  NB         100.0         0.00   0.01   100.00
   YMFVILWVA   M-p9-A  SB         100.0       100.00   0.04     0.00
   DDDDDDDDD   D-p9-D  NB         100.0         0.00   0.00   100.00
   RRRRRRRRR   R-p9-R  NB         100.0         0.00   0.00   100.00
```

## Model Comparison (MHCflurry Data)

```
          Model     Accuracy  Macro_F1  NB_F1  WB_F1  SB_F1
       FFN_Deep      91.9     0.921    0.969  0.880  0.913
 FFN_Jessen2018      90.9     0.911    0.969  0.875  0.891
 CNN_Jessen2018      90.0     0.901    0.963  0.860  0.882
         ResNet      84.3     0.847    0.907  0.792  0.843
           LSTM      83.3     0.836    0.935  0.752  0.822
  Random_Forest      81.1     0.814    0.927  0.723  0.793
```

## 5-Fold Cross-Validation

```
Fold 1 accuracy: 90.86%
Fold 2 accuracy: 89.38%
Fold 3 accuracy: 89.88%
Fold 4 accuracy: 89.29%
Fold 5 accuracy: 88.79%

Mean fold accuracy: 89.64% (±0.79%)
Ensemble accuracy:   89.64%
```

## Validation Summary

| Check | Result |
|-------|--------|
| Canonical anchors (p2=L/M, p9=V/L) → SB | ✅ All correct |
| Negative controls (poly-Asp, poly-Arg) → NB | ✅ 100% confidence |
| MHCflurry concordance (LMAFYLYEV: 11nM → SB) | ✅ 100% confidence |
| Subtle case (LLTDAQRIV: 100nM → SB/WB edge) | ✅ 72%/28% split |
| Poly-Leu (anchors present → SB) | ✅ 98.9% |

## Data Source Comparison

```
Data Source       Best Acc  Macro F1  WB F1   CV Mean
PSSM              94.8%     0.948     0.925   91.1%
MHCflurry         91.9%     0.921     0.880   89.6%
Random Synthetic  65.8%     0.558     0.127   65.4%
```
