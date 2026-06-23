# Extended protein epitope scan — 7 additional immunology targets
options(warn = -1)
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")
suppressMessages({library(keras); library(tensorflow); library(tidyverse)})
setwd('D:/Researching/Peptide epitope')

# ---- Load model ----
cat("Loading Deep FFN model...\n")
model <- load_model_hdf5("models/FFN_Deep.h5")

# ---- BLOSUM62 encoder ----
blosum62 <- matrix(c(
    4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,
   -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,
   -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3,
   -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3,
    0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,
   -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2,
   -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2,
    0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,
   -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3,
   -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,
   -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,
   -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2,
   -2,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,
   -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,
   -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,
    1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2,
    0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,
   -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,
   -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,
    0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4
), nrow=20, ncol=20, byrow=TRUE)
aa_order <- strsplit("ARNDCQEGHILKMFPSTWYV", "")[[1]]
rownames(blosum62) <- colnames(blosum62) <- aa_order
for (i in 1:20) { rng <- range(blosum62[i,]); blosum62[i,] <- (blosum62[i,]-rng[1])/(rng[2]-rng[1]+1e-8) }

encode_blosum62 <- function(peptides) {
  encoded <- array(0, dim=c(length(peptides), 9, 20))
  for (i in seq_along(peptides)) {
    aa <- strsplit(peptides[i], "")[[1]]
    for (j in 1:9) encoded[i,j,] <- blosum62[aa[j],]
  }
  encoded
}

scan_protein <- function(seq, name, model) {
  n_pep <- nchar(seq) - 8
  peptides <- sapply(1:n_pep, function(i) substr(seq, i, i+8))
  x <- encode_blosum62(peptides)
  x <- array_reshape(x, c(dim(x)[1], 9, 20, 1))
  x_flat <- array_reshape(x, c(dim(x)[1], 180))
  preds <- predict(model, x_flat, verbose=0)
  tibble(
    protein=name, start_pos=1:n_pep, end_pos=(1:n_pep)+8, peptide=peptides,
    p2=substr(peptides,2,2), p9=substr(peptides,9,9),
    prob_SB=preds[,3], prob_WB=preds[,2], prob_NB=preds[,1],
    pred_class=c("NB","WB","SB")[apply(preds,1,which.max)]
  ) %>% mutate(
    binding_score=prob_SB*1.0+prob_WB*0.5,
    anchors=p2%in%c("L","M","I","V","A","T","Q")&p9%in%c("V","L","I","A","M","T")
  ) %>% arrange(desc(binding_score))
}

# ---- Protein sequences ----
proteins <- list(

  # Cancer/Testis Antigens
  list(seq="MQAEGRGTGGSTGDADGPGGPGIPDGPGGNAGGPGEAGATGGRGPRGAGAARASGPGGGAPRGPHGGAASGLNGCCRCGARGPESRLLEFYLAMPFATPMEAELARRSLAQDAPPLPVPGVLLKEFTVSGNILTIRLTAADHRQLQLSISSCLQQLSLLMWITQCFLPVFLAQPPSGQRR",
       name="NY-ESO-1 (CTAG1B)", known_epitopes="SLLMWITQC (157-165)"),

  # Wilms Tumor 1
  list(seq="MGSDVRDLNALLPAVPSLGGGGGCALPVSGAAQWAPVLDFAPPGASAYGSLGGPAPPPAPPPPPPPPPHSFIKQEPSWGGAEPHEEQCLSAFTVHFSGQFTGTAGACRYGPFGPPPPSQASSGQARMFPNAPYLPSCLESQPAIRNQGYSTVTFDGTPSYGHTPSHHAAQFPNHSFKHEDPMGQQGSLGEQQYSVPPPVYGCHTPTDSCTGSQALLLRTPYSSDNLYQMTSQLECMTWNQMNLGATLKGVAAGSSSSVKWTEGQSNHSTGYESDNHTTPILCGAQYRIHTHGVFRGIQDVRRVPGVAPTLVRSASETSEKRPFMCAYPGCNKRYFKLSHLQMHSRKHTGEKPYQCDFKDCERRFSRSDQLKRHQRRHTGVKPFQCKTCQRKFSRSDHLKTHTRTHTGKTSEKPFSCRWPSCQKKFARSDELVRHHNMHQRNMTKLQLAL",
       name="WT1", known_epitopes="RMFPNAPYL (126-134), CMTWNQMNL (235-243)"),

  # KRAS (wild-type)
  list(seq="MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM",
       name="KRAS (WT)", known_epitopes="KLVVVGAGGV (5-14), GADGVGKSAL (10-19)"),

  # gp100 / PMEL (melanoma)
  list(seq="MDLVLKRCLLHLAVIGALLAVGATKVPRNQDWLGVSRQLRTKAWNRQLYPEWTEAQRLDCWRGGQVSLKVSNDGPTLIGANASFSIALNFPGSQKVLPDGQVIWVNNTIINGSQVWGGQPVYPQETDDACIFPDGGPCPSGSWSQKRSFVYVWKTWGQYWQVLGGPVSGLSIGTGRAMLGTHTMEVTVYHRRGSRSYVPLAHSSSAFTITDQVPFSVSVSQLRALDGGNKHFLRNQPLTFALQLHDPSGYLAEADLSYTWDFGDSSGTLISRALVVTHTYLEPGPVTAQVVLQAAIPLTSCGSSPVPGTTDGHRPTAEAPNTTAGQVPTTEVVGTTPGQAPTAEPSGTTSVQVPTTEVISTAPVQMPTAESTGMTPEKVPVSEVMGTTLAEMSTPEATGMTPAEVSIVVLSGTTAAQVTTTEWVETTARELPIPEPEGPDASSIMSTESITGSLGPLLDGTATLRLVKRQVPLDCVLYRYGSFSVTLDIVQGIESAEILQAVPSGEGDAFELTVSCQGGLPKEACMEISSPGCQPPAQRLCQPVLPSPACQLVLHQILKGGSGTYCLNVSLADTNSLAVVSTQLIMPGQEAGLGQVPLIVGILLVLMAVVLASLIYRRRLMKQDFSVPQLPHSSSHWLRLPRIFCSCPIGENSPLLSGQQV",
       name="gp100/PMEL", known_epitopes="IMDQVPFSV (209-217), YLEPGPVTA (280-288)"),

  # Tyrosinase (melanoma)
  list(seq="MLLAVLYCLLWSFQTSAGHFPRACVSSKNLMEKECCPPWSGDRSPCGQLSGRGSCQNILLSNAPLGPQFPFTGVDDRESWPSVFYNRTCQCSGNFMGFNCGNCKFGFWGPNCTERRLLVRRNIFDLSAPEKDKFFAYLTLAKHTISSDYVIPIGTYGQMKNGSTPMFNDINIYDLFVWMHYYVSMDALLGGSEIWRDIDFAHEAPAFLPWHRLFLLRWEQEIQKLTGDENFTIPYWDWRDAEKCDICTDEYMGGQHPTNPNLLSPASFFSSWQIVCSRLEEYNSHQSLCNGTPEGPLRRNPGNHDKSRTPRLPSSADVEFCLSLTQYESGSMDKAANFSFRNTLEGFASPLTGIADASQSSMHNALHIYMNGTMSQVQGSANDPIFLLHHAFVDSIFEQWLRRHRPLQEVYPEANAPIGHNRESYMVPFIPLYRNGDFFISSKDLGYDYSYLQDSDPDSFQDYIKSYLEQASRIWSWLLGAAMVGAVLTALLAGLVSLLCRHKRKQLPEEKQPLLMEKEDYHSLYQSHL",
       name="Tyrosinase", known_epitopes="YMDGTMSQV (369-377)"),

  # Influenza M1 (excellent positive control)
  list(seq="MSLLTEVETYVLSIIPSGPLKAEIAQRLEDVFAGKNTDLEVLMEWLKTRPILSPLTKGILGFVFTLTVPSERGLQRRRFVQNALNGNGDPNNMDKAVKLYRKLKREITFHGAKEISLSYSAGALASCMGLIYNRMGAVTTEVAFGLVCATCEQIADSQHRSHRQMVTTTNPLIRHENRMVLASTTAKAMEQMAGSSEQAAEAMEVASQARQMVQAMRTIGTHPSSSAGLKNDLLENLQAYQKRMGVQMQRFK",
       name="Influenza M1", known_epitopes="GILGFVFTL (58-66) - classic immunodominant"),

  # CMV pp65 (viral, immunodominant)
  list(seq="MESRGRRCPEMISVLGPISGHVLKAVFSRGDTPVLPHETRLLQTGIHVRVSQPSLILVSQYTPDSTPCHRGDNQLQVQHTYFTGSEVENVSVNVHNPTGRSICPSQEPMSIYVYALPLKMLNIPSINVHHYPSAAERKHRHLPVADAVIHASGKQMWQARLTVSGLAWTRQQNQWKEPDVYYTSAFVFPTKDVALRHVVCAHELVCSMENTRATKMQVIGDQYVKVYLESFCEDVPSGKLFMHVTLGSDVEEDLTMTRNPQPFMRPHERNGFTVLCPKNMIIKPGKISHIMLDVAFTSHEHFGLLCPKSIPGLSISGNLLMNGQQIFLEVQAIRETVELRQYDPVAALFFFDIDLLLQRGPQYSEHPTFTSQYRIQGKLEYRHTWDRHDEGAAQGDDDVWTSGSDSDEELVTTERKTPRVTGGGAMAGASTSAGRKRKSASSATACTSGVMTRGRLKAESTCRPEEDTDEDSDNEIHNPAVFTWPPWQAGILARNLVPMVATVQGQNLKYQEFFWDANDIYRIFAELEGVWQPAAQPKRRRHRQDALPGPCIASTPKKHRG",
       name="CMV pp65", known_epitopes="NLVPMVATV (495-503) - immunodominant")
)

# ---- Run scans ----
cat("\n==============================================================\n")
cat("  EXTENDED PROTEIN EPITOPE SCAN — HLA-A*02:01\n")
cat("  Model: Deep FFN (MHCflurry, 91.9% accuracy)\n")
cat("==============================================================\n\n")

all_results <- list()
for (p in proteins) {
  cat(sprintf("Scanning %s (%d aa)...\n", p$name, nchar(p$seq)))
  res <- scan_protein(p$seq, p$name, model)
  all_results[[p$name]] <- res
}

# ---- Top candidates per protein ----
cat("\n")
cat("==============================================================\n")
cat("  TOP 5 EPITOPE CANDIDATES PER PROTEIN\n")
cat("==============================================================\n\n")

for (name in names(all_results)) {
  cat(sprintf("--- %s ---\n", name))
  res <- all_results[[name]]
  top <- head(res, 5)
  top %>% mutate(
    score=sprintf("%.3f",binding_score),
    `P(SB)`=sprintf("%.1f%%",prob_SB*100),
    `P(WB)`=sprintf("%.1f%%",prob_WB*100),
    pos=sprintf("%d-%d",start_pos,end_pos),
    A=ifelse(anchors,"+","-")
  ) %>% select(pos,peptide,pred_class,score,`P(SB)`,`P(WB)`,A) %>%
    print(n=5,width=100)
  cat("\n")
}

# ---- Known epitope check ----
cat("==============================================================\n")
cat("  KNOWN EPITOPE VALIDATION\n")
cat("==============================================================\n\n")

known_df <- tribble(
  ~protein,            ~peptide,       ~description,
  "NY-ESO-1 (CTAG1B)", "SLLMWITQC",   "NY-ESO-1 157-165 (immunodominant)",
  "WT1",               "RMFPNAPYL",   "WT1 126-134 (validated A*02:01)",
  "WT1",               "CMTWNQMNL",   "WT1 235-243 (validated A*02:01)",
  "Influenza M1",      "GILGFVFTL",   "M1 58-66 (classic immunodominant)",
  "CMV pp65",          "NLVPMVATV",   "pp65 495-503 (immunodominant)",
  "gp100/PMEL",        "IMDQVPFSV",   "gp100 209-217 (validated A*02:01)",
  "gp100/PMEL",        "YLEPGPVTA",   "gp100 280-288 (validated A*02:01)",
  "Tyrosinase",        "YMDGTMSQV",   "Tyrosinase 369-377 (validated)",
  "KRAS (WT)",         "KLVVVGAGGV",  "KRAS 5-14 (reported binder)"
)

cat("Predicting known epitopes...\n")
x_k <- encode_blosum62(known_df$peptide)
x_k <- array_reshape(x_k, c(dim(x_k)[1], 9, 20, 1))
x_k_flat <- array_reshape(x_k, c(dim(x_k)[1], 180))
pk <- predict(model, x_k_flat, verbose=0)

known_results <- known_df %>%
  mutate(
    pred=c("NB","WB","SB")[apply(pk,1,which.max)],
    P_SB=round(pk[,3]*100,1),
    P_WB=round(pk[,2]*100,1),
    P_NB=round(pk[,1]*100,1),
    match=ifelse(pred %in% c("SB","WB"), "CONFIRMED", "REVIEW")
  )

known_results %>%
  select(protein, peptide, pred, P_SB, P_WB, match, description) %>%
  print(n=15, width=130)

# ---- Stats ----
cat("\n==============================================================\n")
cat("  SCAN STATISTICS\n")
cat("==============================================================\n\n")

stats <- bind_rows(all_results) %>% group_by(protein) %>% summarise(
  length=first(nchar(protein)),
  windows=n(),
  SB=sum(pred_class=="SB"),
  WB=sum(pred_class=="WB"),
  NB=sum(pred_class=="NB"),
  .groups="drop"
) %>% mutate(
  `SB%`=sprintf("%.1f%%",100*SB/windows),
  `WB%`=sprintf("%.1f%%",100*WB/windows)
)

stats %>% print(n=20, width=100)

total <- bind_rows(all_results)
cat(sprintf("\nTotal: %d proteins, %d peptides scanned\n", length(all_results), nrow(total)))
cat(sprintf("SB: %d (%.1f%%) | WB: %d (%.1f%%) | NB: %d (%.1f%%)\n",
    sum(total$pred_class=="SB"), 100*sum(total$pred_class=="SB")/nrow(total),
    sum(total$pred_class=="WB"), 100*sum(total$pred_class=="WB")/nrow(total),
    sum(total$pred_class=="NB"), 100*sum(total$pred_class=="NB")/nrow(total)))

# ---- Save ----
write_csv(bind_rows(all_results), "data/protein_epitope_scan_extended.csv")
cat(sprintf("\nSaved: data/protein_epitope_scan_extended.csv (%d rows)\n", nrow(bind_rows(all_results))))
