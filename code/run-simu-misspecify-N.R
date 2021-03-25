#### Prepare data with chromsome 22 only####

library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")

load("data/ind_gwas_val_test.RData")

(NCORES <- parallelly::availableCores() - 1L)

chr <- 22
ind.chr <- which(ukb$map$chromosome == chr)
corr0 <- readRDS(paste0("data/corr/chr", chr, ".rds"))
corr <- runonce::save_run(as_SFBM(corr0, "tmp-data/corr_chr22"),
                          "tmp-data/corr_chr22.rds")
stopifnot(length(ind.chr) == ncol(corr))
dim(corr)  # 16410 x 16410

# big_copy(ukb$genotypes, ind.col = ind.chr,
#          backingfile = "tmp-data/dosage_chr22")$save()
G <- big_attach("tmp-data/dosage_chr22.rds")

# prepare data in bed format
data <- snp_fake(n = nrow(G), m = 1)
data$map <- ukb$map[ind.chr, ]
data$map$genetic.dist <- 0
data$map$marker.ID <- data$map$rsid
data$genotypes <- G
# snp_writeBed(data, ind.row = ind.val, bedfile = "tmp-data/chr22.bed")

library(future.batchtools)
NCORES <- 15
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 1, mem = "100g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("results/simu-misN")


furrr::future_walk(1:10, function(ic) {

  res_file <- paste0("results/simu-misN/res_", ic, ".rds")
  if (file.exists(res_file)) return(NULL)

  #### Run GWAS and prepare sumstats ####

  list_ind <- list(ind.gwas,
                   sort(sample(ind.gwas, length(ind.gwas) * 0.8)),
                   sort(sample(ind.gwas, length(ind.gwas) * 0.6)))

  y <- snp_simuPheno(G, h2 = 0.2, M = 2000, ncores = NCORES)$pheno

  # GWAS to get sumstats
  gwas_set <- sample(rep_len(c(1, 1, 2, 3), ncol(G)))
  df_beta <- data.frame(beta = rep(NA, ncol(G)),
                        beta_se = NA, n_eff = NA, lpval = NA)

  for (k in 1:3) {
    ind_set <- which(gwas_set == k)
    ind.gwas.sub <- list_ind[[k]]
    gwas <- big_univLinReg(G, y[ind.gwas.sub], ind.train = ind.gwas.sub,
                           ind.col = ind_set, ncores = NCORES)
    df_beta$beta[ind_set]    <- gwas$estim
    df_beta$beta_se[ind_set] <- gwas$std.err
    df_beta$n_eff[ind_set]   <- length(ind.gwas.sub)
    df_beta$lpval[ind_set]   <- -predict(gwas)
  }
  with(df_beta, min(pchisq((beta / beta_se)^2, df = 1, lower.tail = FALSE)))

  df_beta2 <- df_beta; df_beta2$n_eff <- max(df_beta$n_eff)

  # Quality control plot from the LDpred2 paper
  sd_val <- sqrt(big_colstats(G, ind.row = ind.val, ncores = NCORES)$var)
  sd_ss <- sd(y) / with(df_beta2, sqrt(n_eff * beta_se^2 + beta^2))

  qplot(sd_val, sd_ss, alpha = I(0.6), color = as.factor(df_beta$n_eff)) +
    theme_bigstatsr(0.9) +
    coord_equal() +
    geom_abline(linetype = 2, color = "red") +
    labs(x = "Standard deviations in the validation set", color = "N",
         y = "Standard deviations derived from the summary statistics") +
    theme(legend.position = c(0.25, 0.8))
  # ggsave("figures/simu-qc-plot.pdf", width = 7, height = 7)


  #### Run LDpred2 ####

  run_ldpred2 <- function(df_beta) {

    # LDSc reg
    print(ldsc <- snp_ldsc2(corr0, df_beta))
    h2_est <- ldsc[["h2"]]

    # LDpred-inf
    beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
    pred_inf <- big_prodVec(G, beta_inf, ind.row = ind.test)
    print(r21 <- cor(pred_inf, y[ind.test])**2)

    # LDpred-auto
    auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est, vec_p_init = 0.1,
                             burn_in = 200, num_iter = 200)
    print(c(auto[[1]]$h2_est, auto[[1]]$p_est))
    pred_auto <- big_prodVec(G, auto[[1]]$beta_est, ind.row = ind.test)
    print(r23 <- cor(pred_auto, y[ind.test])**2)

    ## LDpred2-grid
    (h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4))
    (p_seq <- signif(seq_log(1e-4, 1, length.out = 9), 2))
    (params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = FALSE))

    beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)
    pred_grid <- big_prodMat(G, beta_grid)
    ind.best <- which.max(apply(pred_grid[ind.val, ], 2, cor, y = y[ind.val]))
    pred_gibbs <- pred_grid[ind.test, ind.best]
    print(r22 <- cor(pred_gibbs, y[ind.test])**2)

    c(r21, r22, r23)
  }

  # ldpred2 <- run_ldpred2(df_beta)
  # mis_ldpred2 <- run_ldpred2(df_beta2)


  #### Run lassosum2 ####

  run_lassosum2 <- function(df_beta) {
    beta_lassosum <- snp_lassosum2(corr, df_beta, s = 1:5 / 5, ncores = NCORES)
    pred_grid <- big_prodMat(G, beta_lassosum, ncores = NCORES)
    score <- big_univLinReg(as_FBM(pred_grid[ind.val, ]), y[ind.val])$score
    pred_lassosum <- pred_grid[ind.test, which.max(score)]
    print(cor(pred_lassosum, y[ind.test])**2)
  }

  # lassosum2 <- run_lassosum2(df_beta)
  # mis_lassosum2 <- run_lassosum2(df_beta2)


  #### Run C+T ####

  run_CT <- function(df_beta) {

    CHR <- as.integer(data$map$chromosome)
    POS <- data$map$physical.pos
    all_keep <- snp_grid_clumping(G, CHR, POS, ind.row = ind.val,
                                  grid.base.size = 200,
                                  grid.thr.r2 = c(0.05, 0.2, 0.8),
                                  lpS = df_beta$lpval, ncores = NCORES)
    attr(all_keep, "grid")

    multi_PRS <- snp_grid_PRS(G, all_keep, df_beta$beta, df_beta$lpval,
                              ind.row = ind.val, n_thr_lpS = 50, ncores = NCORES)
    dim(multi_PRS)  ## 150 C+T scores

    library(dplyr)
    grid2 <- attr(all_keep, "grid") %>%
      mutate(thr.lp = list(attr(multi_PRS, "grid.lpS.thr")), id = row_number()) %>%
      tidyr::unnest(cols = "thr.lp")
    s <- nrow(grid2)
    grid2$score <- big_univLinReg(multi_PRS, y[ind.val])$score
    max_prs <- grid2 %>% arrange(desc(score)) %>% slice(1:10) %>% print() %>% slice(1)

    ind.keep <- unlist(purrr::map(all_keep, max_prs$id))
    pred_CT <- snp_PRS(G, df_beta$beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
                       lpS.keep = df_beta$lpval[ind.keep], thr.list = max_prs$thr.lp)
    print(cor(pred_CT, y[ind.test])[[1]]**2)
  }

  # CT <- run_CT(df_beta)


  #### Run lassosum ####

  run_lassosum <- function(df_beta) {

    t <- with(df_beta, beta / beta_se)
    n <- df_beta$n_eff

    library(lassosum)
    doParallel::registerDoParallel(cl <- parallel::makeCluster(NCORES))
    system.time(
      out <- lassosum.pipeline(
        cor = t / sqrt(n - 2 + t^2),
        snp = data$map$marker.ID,
        A1 = data$map$allele1,
        A2 = data$map$allele2,
        exclude.ambiguous = FALSE,
        test.bfile = "tmp-data/chr22",
        LDblocks = "EUR.hg19",
        cluster = cl,
        destandardize = TRUE
      )
    ) # < 3 min
    parallel::stopCluster(cl)

    all_score <- apply(do.call("cbind", out$pgs), 2, cor, y = y[ind.val])
    best_beta <- do.call("cbind", out$beta)[, which.max(all_score)]
    pred_lassosum <- big_prodVec(G, best_beta, ind.row = ind.test)
    print(cor(pred_lassosum, y[ind.test])**2)
  }

  # lassosum <- run_lassosum(df_beta)
  # mis_lassosum <- run_lassosum(df_beta2)


  #### SBayesR -> do not run because always diverged ####

  # gctb <- "../paper-ldpred2/tmp-data/gctb_2.02_Linux/gctb"
  #
  # # Compute LD
  # if (!file.exists(paste0("tmp-data/ldm_22.ldm.shrunk.bin"))) {
  #   system.time(
  #     system(glue::glue(
  #       "{gctb} --bfile tmp-data/chr22",
  #       " --make-shrunk-ldm",
  #       " --out tmp-data/ldm_22"
  #     ))
  #   ) # 18 min
  # }
  #
  # # Compute SBayesR
  # obj.bed <- bigsnpr::bed("tmp-data/chr22.bed")
  # af <- bigsnpr::bed_MAF(obj.bed, ncores = NCORES)$af
  #
  # tmp <- tempfile(tmpdir = "tmp-data", fileext = ".ma")
  # library(dplyr)
  # df_beta %>%
  #   bind_cols(data$map) %>%
  #   transmute(SNP = rsid, A1 = allele1, A2 = allele2, freq = af,
  #             b = beta, se = beta_se, p = 10^-lpval, N = n_eff) %>%
  #   bigreadr::fwrite2(tmp, sep = " ") %>%
  #   readLines(n = 5) %>%
  #   writeLines()
  #
  # res_file <- "tmp-data/sbayesr_chr22"
  #
  # system(glue::glue(
  #   gctb,
  #   " --ldm tmp-data/ldm_{chr}.ldm.shrunk",
  #   " --sbayes R --pi 0.95,0.02,0.02,0.01 --gamma 0.0,0.01,0.1,1",
  #   # " --sbayes R --pi 0.9,0.1 --gamma 0.0,0.1",
  #   # " --p-value 0.4 --rsq 0.95",
  #   " --gwas-summary {tmp}",
  #   " --chain-length 10000 --burn-in 2000",
  #   " --out {res_file} --out-freq 100"
  # ))
  #
  # file.remove(tmp)
  #
  # library(dplyr)
  # head(res_sbayesr <- bigreadr::fread2(paste0(res_file, ".snpRes")))
  #
  # ind <- match(res_sbayesr$Name, data$map$rsid)
  # stopifnot(!anyNA(ind))
  # stopifnot(all.equal(data$map$allele1[ind], res_sbayesr$A1))
  # stopifnot(all.equal(data$map$allele2[ind], res_sbayesr$A2))
  #
  # pred_sbayesr <- big_prodVec(G, res_sbayesr$A1Effect, ind.row = ind.test,
  #                             ind.col = ind)
  # cor(pred_sbayesr, y[ind.test])**2


  #### Run PRS-CS ####

  run_prscs <- function(df_beta) {

    tmp <- tempfile(tmpdir = "tmp-data", fileext = ".txt")
    library(dplyr)
    df_beta %>%
      bind_cols(data$map) %>%
      transmute(SNP = rsid, A1 = allele1, A2 = allele2,
                BETA = beta, P = 10^-lpval) %>%
      bigreadr::fwrite2(tmp, sep = "\t") %>%
      readLines(n = 5) %>%
      writeLines()

    on.exit(file.remove(tmp), add = TRUE)

    prefix <- paste0("tmp-data/TMP_PRSCS", ic)

    PHI <- c(NA, 1e-4, 1e-2, 1)
    prscs <- "PRScs/PRScs.py"
    for (phi in PHI) {
      system(glue::glue(
        "OMP_NUM_THREADS=", NCORES[[1]],
        " python3 {prscs}",
        " --ref_dir=ldblk_1kg_eur",
        " --bim_prefix=tmp-data/chr22",
        " --sst_file={tmp}",
        " --n_gwas={max(df_beta$n_eff)}",
        if (is.na(phi)) "" else " --phi={phi}",
        " --chrom=22",
        " --out_dir={prefix}",
        " --n_iter=400 --n_burnin=200"
      ))
    }
    on.exit(unlink(paste0(prefix, "*")), add = TRUE)

    get_pred <- function(phi) {
      file <- paste0(prefix, "_pst_eff_a1_b0.5_phi", phi, "_chr22.txt")
      res <- bigreadr::fread2(file)
      betas <- rep(0, ncol(G))
      betas[match(res$V2, data$map$rsid)] <- res$V6
      big_prodVec(G, betas, ncores = NCORES)
    }

    all_pred <- sapply(c("1e-04", "1e-02", "1e+00"), get_pred)
    ind.best <- which.max(apply(all_pred[ind.val, ], 2, cor, y = y[ind.val]))
    pred_prscs <- all_pred[ind.test, ind.best]
    print(r21 <- cor(pred_prscs, y[ind.test])**2)

    pred_prscs_auto <- get_pred("auto")[ind.test]
    print(r22 <- cor(pred_prscs_auto, y[ind.test])**2)

    c(r21, r22)
  }

  # prscs <- run_prscs(df_beta)
  # mis_prscs <- run_prscs(df_beta2)


  #### Run them all ####

  res <- tibble::tribble(
    ~ Method,                              ~ misspecified,     ~ r2,
    c("LDpred2-inf", "LDpred2", "LDpred2-auto"), FALSE, run_ldpred2(df_beta),
    c("LDpred2-inf", "LDpred2", "LDpred2-auto"), TRUE,  run_ldpred2(df_beta2),
    "lassosum2", FALSE, run_lassosum2(df_beta),
    "lassosum2", TRUE,  run_lassosum2(df_beta2),
    "C+T", FALSE, run_CT(df_beta),
    "C+T", TRUE,  run_CT(df_beta2),
    "lassosum", FALSE, run_lassosum(df_beta),
    "lassosum", TRUE,  run_lassosum(df_beta2),
    c("PRS-CS", "PRS-CS-auto"), FALSE, run_prscs(df_beta),
    c("PRS-CS", "PRS-CS-auto"), TRUE,  run_prscs(df_beta2),
  )

  saveRDS(tidyr::unnest(res, cols = c(Method, r2)), res_file)
})

library(dplyr)
all_res <- list.files("results/simu-misN", full.names = TRUE) %>%
  purrr::map_dfr(readRDS) %>%
  group_by(Method, misspecified) %>%
  summarise(r2 = {
    boot <- replicate(1e4, mean(sample(r2, replace = TRUE)))
    q <- quantile(boot, c(0.025, 0.975))
    list(c(mean = mean(boot), inf = q[[1]], sup = q[[2]]))
  }) %>%
  ungroup() %>%
  tidyr::unnest_wider("r2") %>%
  print(n = Inf)
#    Method       misspecified   mean    inf    sup
#  1 C+T          FALSE        0.139  0.137  0.140
#  2 C+T          TRUE         0.139  0.137  0.141
#  3 lassosum     FALSE        0.175  0.174  0.177
#  4 lassosum     TRUE         0.174  0.172  0.175
#  5 lassosum2    FALSE        0.180  0.179  0.182
#  6 lassosum2    TRUE         0.175  0.174  0.177
#  7 LDpred2      FALSE        0.178  0.177  0.180
#  8 LDpred2      TRUE         0.115  0.109  0.121
#  9 LDpred2-auto FALSE        0.176  0.175  0.178
# 10 LDpred2-auto TRUE         0.0566 0.0473 0.0667
# 11 LDpred2-inf  FALSE        0.177  0.175  0.178
# 12 LDpred2-inf  TRUE         0.0680 0.0607 0.0756
# 13 PRS-CS       FALSE        0.154  0.152  0.156
# 14 PRS-CS       TRUE         0.154  0.152  0.156
# 15 PRS-CS-auto  FALSE        0.153  0.151  0.155
# 16 PRS-CS-auto  TRUE         0.153  0.151  0.155

library(ggplot2)
ggplot(all_res,
       aes(Method, mean, fill = misspecified)) +
  bigstatsr::theme_bigstatsr(0.8) +
  scale_fill_manual(values = c("#E69F00", "#0072B2")) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = inf, ymax = sup),
                position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  scale_y_continuous(breaks = 0:4 / 20, minor_breaks = 0:20 / 100) +
  labs(x = "Method", y = "Squared correlation between PGS and phenotype",
       fill = "Misspecified\nsample size")
# ggsave("figures/simu-misN.pdf", width = 11, height = 6)
