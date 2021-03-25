library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes
CHR <- as.integer(ukb$map$chromosome)

load("data/ind_val_test.RData")

library(dplyr)
files <- tibble::tibble(basename = list.files("data/sumstats")) %>%
  mutate(
    res_file = file.path("results/lassosum2", basename),
    pheno_file = file.path("data/pheno",  basename),
    gwas_file = file.path("data/sumstats", basename)) %>%
  print()
all(file.exists(files$pheno_file)) & all(file.exists(files$gwas_file))
bigassertr::assert_dir("results/lassosum2")

files_sub <- files %>%
  filter(!file.exists(res_file)) %>%
  print()

library(future.batchtools)
NCORES <- 30
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 2, mem = "100g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(files_sub[-1], function(res_file, pheno_file, gwas_file) {

  y <- readRDS(pheno_file)
  sumstats <- readRDS(gwas_file)

  ind.val2 <- ind.val[!is.na(y[ind.val])]

  tmp <- tempfile(tmpdir = "tmp-data")
  on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

  for (chr in 1:22) {

    ## indices in 'sumstats'
    ind.chr <- which(sumstats$chr == chr)
    ## indices in 'G'
    ind.chr2 <- sumstats$`_NUM_ID_`[ind.chr]
    ## indices in 'corr'
    ind.chr3 <- match(ind.chr2, which(CHR == chr))

    corr0 <- readRDS(paste0("data/corr/chr", chr, ".rds"))[ind.chr3, ind.chr3]

    if (chr == 1) {
      df_beta <- sumstats[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
      corr <- as_SFBM(corr0, tmp)
    } else {
      df_beta <- rbind(df_beta, sumstats[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")])
      corr$add_columns(corr0, nrow(corr))
    }
  }

  # lassosum2
  beta_lassosum <- snp_lassosum2(corr, df_beta, ncores = NCORES)
  params <- attr(beta_lassosum, "grid_param")

  bigparallelr::set_blas_ncores(NCORES)
  ind <- which(rowSums(beta_lassosum != 0) > 0)
  pred_lassosum <- big_prodMat(G, beta_lassosum[ind, ], ind.row = ind.val2,
                               ind.col = df_beta[["_NUM_ID_"]][ind])
  params$score <- big_univLogReg(as_FBM(pred_lassosum), y[ind.val2])$score

  library(dplyr)
  best_beta_lassosum <- params %>%
    mutate(id = row_number()) %>%
    arrange(desc(score)) %>%
    slice(1) %>%
    pull(id) %>%
    beta_lassosum[, .]

  # compute prediction for test set
  ind <- which(best_beta_lassosum != 0)
  pred_test <- big_prodVec(G, best_beta_lassosum[ind], ind.row = ind.test,
                           ind.col = df_beta[["_NUM_ID_"]][ind],
                           ncores = NCORES)

  # save results
  res <- list(pred = pred_test, params = params)
  saveRDS(res, res_file)
})
