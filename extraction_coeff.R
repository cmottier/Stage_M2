library(MCMCvis)

for (modele in c(0,1,3)) {
  load(paste0("out_modele", modele, "_5km2_2010-2024.RData"))
  resum_coeffs <- MCMCsummary(out$samples)
  save(resum_coeffs, file = paste0("resume_M", modele, ".RData"))
  rm(out)
  gc()
}
