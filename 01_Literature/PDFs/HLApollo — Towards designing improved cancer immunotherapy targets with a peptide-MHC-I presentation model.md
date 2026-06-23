# Towards designing improved cancer immunotherapy targets with a peptide-MHC-I presentation model, HLApollo

**Authors**: W. Thrift, Nicolas W. Lounsbury, Quade Broadwell, A. Heidersbach, Emily C Freund, Yassan Abdolazimi, Q. Phung, Jieming Chen, A. Capietto, A. Tong, Christopher M. Rose, Craig Blanchette, J. Lill, Benjamin Haley, L. Delamarre, R. Bourgon, Kaiyuan Liu, S. Jhunjhunwala
**Year**: 2024
**Journal**: Nature Communications
**DOI**: 10.1038/s41467-024-54887-7
**Citations**: 13
**PMC**: https://pmc.ncbi.nlm.nih.gov/articles/PMC11686168

## Abstract

Based on the success of cancer immunotherapy, personalized cancer vaccines have emerged as a leading oncology treatment. Antigen presentation on MHC class I (MHC-I) is crucial for the adaptive immune response to cancer cells, necessitating highly predictive computational methods to model this phenomenon. Here, we introduce HLApollo, a transformer-based model for peptide-MHC-I (pMHC-I) presentation prediction, leveraging the language of peptides, MHC, and source proteins. HLApollo provides end-to-end treatment of MHC-I sequences and deconvolution of multi-allelic data, using a negative-set switching strategy to mitigate misassigned negatives in unlabelled ligandome data. HLApollo shows a 12.65% increase in average precision (AP) on ligandome data and a 4.1% AP increase on immunogenicity test data compared to next-best models. Incorporating protein features from protein language models yields further gains and reduces the need for gene expression measurements. Guided by clinical use, we demonstrate pan-allelic generalization which effectively captures rare alleles in underrepresented ancestries.

**Key Results**: 78.7% AP (presentation), 12.65% better than next best; negative-set switching (+11.79% AP); ESM-1b integration; pan-allelic generalization.

---

**Source**: Semantic Scholar API | **Downloaded**: 2026-06-20 via medsci
