# Quick-start verification script
options(warn = -1)
setwd('D:/Researching/Peptide epitope')
source('peptide_mhc_binding_study.R')

cat('\nRunning quickstart_jessen_ffn()...\n\n')

result <- tryCatch(
  quickstart_jessen_ffn(),
  error = function(e) { cat('ERROR:', e$message, '\n'); NULL }
)

if (!is.null(result)) {
  cat('\n=== Quick-start Results ===\n')
  cat(sprintf('Accuracy:  %.1f%%\n', result$eval$accuracy * 100))
  cat(sprintf('Macro F1:  %.3f\n', result$eval$macro_f1))
  for(cls in names(result$eval$per_class)) {
    m <- result$eval$per_class[[cls]]
    cat(sprintf('  %s: P=%.3f R=%.3f F1=%.3f\n', cls, m['Precision'], m['Recall'], m['F1']))
  }
}
