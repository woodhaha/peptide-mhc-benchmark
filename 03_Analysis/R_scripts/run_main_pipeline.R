# Full pipeline execution — auto-detects real vs synthetic data
options(warn = -1)
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "2")
setwd('D:/Researching/Peptide epitope')

cat('\n==============================================================\n')
cat(sprintf('  FULL PIPELINE START: %s\n', Sys.time()))
cat('==============================================================\n\n')

suppressMessages(source('peptide_mhc_binding_study.R'))

# Auto-detect data source: use real-labeled CSV if available
REAL_DATA_PATH <- "data/real_peptides.csv"
if (file.exists(REAL_DATA_PATH) && file.info(REAL_DATA_PATH)$size > 1000) {
  cat(sprintf('\n*** Using REAL MHCflurry-labeled data: %s ***\n\n', REAL_DATA_PATH))
  study_results <- main(data_source = "custom", data_path = REAL_DATA_PATH)
} else {
  cat('\n*** Using synthetic data (no real labels found) ***\n\n')
  study_results <- main(data_source = "jessen2018")
}

cat(sprintf('\nPipeline finished: %s\n', Sys.time()))
cat(sprintf('Exit status: SUCCESS\n'))
