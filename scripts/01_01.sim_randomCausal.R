library("argparser")
library("GWASBrewer")
library("TwoSampleMR")
library("susieR")
library("hapsim")
library("dplyr")
library("Matrix")
library("data.table")
library("matrixStats")
library("coloc")
library("stringr")
library("finemappingMR")
library("gsmr2")

p <- arg_parser("Script that randomly samples causal variants and effects, generates summary statistics using UKBB
                Price Lab LD for exposure and outcome at true causal effect sizes
                of 0, 0.1, and 0.2. The script then runs traditional MR
                methods (IVW, median, Egger, RAPS, SMR) using LD-clumped variants and
                GSMR using IVs selected from SuSiE credible sets.
                This script takes the number of causal variants you want to simulate
                for the exposure, as well as the number of pleiotropic effects.")

p <- add_argument(p, "--run", help = "Run number of simulation and number to use
                  as seed for the simulation")
p <- add_argument(p, "--N_exposure", help = "Size of the exposure study to simulate", default = 33687)
p <- add_argument(p, "--N_outcome", help = "Size of the outcome study to simulate", default = 420607)
p <- add_argument(p, "--num_causal", help = "Number of causal variants to simulate (1-8)")
p <- add_argument(p, "--num_pleio", help = "Number of pleiotropic causal variants to simulate (1-7)")

p <- add_argument(p, "--ld_file", help = "Name of the LD file to read in and use")
p <- add_argument(p, "--af_file", help = "Name of the file containing the allele frequency information. 
                  Must have a column af_exposure.")

p <- add_argument(p, "--output_elboConv", help = "Name of file to write ELBO output that was used for convergence checks to")
p <- add_argument(p, "--output_res", help = "Name of file to output results to")
p <- add_argument(p, "--output_tradMR", help = "Name of file to output traditional MR results to")
p <- add_argument(p, "--output_singleVar", help = "Name of file to output the single-variant summary statistics to. This file includes the true joint beta as well.")
p <- add_argument(p, "--output_MRlink2_exp_path", help = "Path to output the exposure MR-link2 formatted summary statistics to.")
p <- add_argument(p, "--output_MRlink2_out_path", help = "Path to output the outcome MR-link2 formatted summary statistics to.")
p <- add_argument(p, "--output_cojo_exp_path", help = "Path to output the exposure COJO formatted summary statistics to.")
p <- add_argument(p, "--output_cojo_out_path", help = "Path to output the outcome COJO formatted summary statistics to.")
p <- add_argument(p, "--output_priorVar", help = "Name of file to output the prior variance estimates to")
p <- add_argument(p, "--output_posteriorMean", help = "Name of file to output the posterior mean estimates to")
p <- add_argument(p, "--output_expectedLD", help = "Name of file to output the expected LD of the causal SNPs to. First
                  5 rows/columns are the causal SNPs for X; last 2 rows/columns are the causal SNPs for Y")

#Modified MR function that fixes issue for for MR RAPS
mr_modified <- function (dat, 
                         parameters = default_parameters(), 
                         method_list = subset(mr_method_list(), use_by_default)$obj){
  
  mr_raps_modified <- function (b_exp, b_out, se_exp, se_out,parameters) 
  {
    out <- try(suppressMessages(mr.raps::mr.raps(b_exp, b_out, se_exp, se_out,
                                                 over.dispersion = parameters$over.dispersion, 
                                                 loss.function = parameters$loss.function,
                                                 diagnosis = FALSE)),
               silent = T)
    
    # The estimated overdispersion parameter is very small. Consider using the simple model without overdispersion
    # When encountering such warning, change the over.dispersion as 'FASLE'
    
    if ('try-error' %in% class(out))
    {
      output = list(b = NA, se = NA, pval = NA, nsnp = NA)
    }
    else
    {
      output = list(b = out$beta.hat, se = out$beta.se, 
                    pval = pnorm(-abs(out$beta.hat/out$beta.se)) * 2, nsnp = length(b_exp))
    }
    return(output)
  }
  
  method_list_modified <- stringr::str_replace_all(method_list, "mr_raps","mr_raps_modified")
  
  mr_tab <- plyr::ddply(dat, c("id.exposure", "id.outcome"),function(x1)
  {
    x <- subset(x1, mr_keep)
    
    if (nrow(x) == 0) {
      message("No SNPs available for MR analysis of '", x1$id.exposure[1], "' on '", x1$id.outcome[1], "'")
      return(NULL)
    }
    else {
      message("Analysing '", x1$id.exposure[1], "' on '", x1$id.outcome[1], "'")
    }
    res <- lapply(method_list_modified, function(meth)
    {
      get(meth)(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome, parameters)
    }
    )
    
    methl <- mr_method_list()
    mr_tab <- data.frame(outcome = x$outcome[1], exposure = x$exposure[1], 
                         method = methl$name[match(method_list, methl$obj)], 
                         nsnp = sapply(res, function(x) x$nsnp), 
                         b = sapply(res, function(x) x$b), 
                         se = sapply(res, function(x) x$se), 
                         pval = sapply(res, function(x) x$pval))
    
    mr_tab <- subset(mr_tab, !(is.na(b) & is.na(se) & is.na(pval)))
    
    return(mr_tab)
  }
  )
  return(mr_tab)
}

argv <- parse_args(p)

#Read in the LD
ld_curr <- readRDS(argv$ld_file)

#Read in the allele frequencies
df_af <- fread(argv$af_file, header = TRUE, data.table = FALSE)
var_pos <- unlist(strsplit(df_af$vars_hg38, split = ":", fixed = TRUE))[c(FALSE, TRUE, FALSE, FALSE)]

num_total_vars <- nrow(ld_curr)

gamma_vec <- c(0, 0.1, 0.2)
gamma_est <- c()
gamma_var <- c()
gamma_var_total_var <- c()
iter_num <- c()
p_min <- c()
pph4_reg <- c()
elbo_list <- list()
elbo_list_full <- list()
trad_mr_list <- list()

V_x_mat <- data.frame(gamma1_x = rep(0, 10),
                      gamma2_x = rep(0, 10),
                      gamma3_x = rep(0, 10))

V_y_mat <- data.frame(gamma1_y = rep(0, 10),
                      gamma2_y = rep(0, 10),
                      gamma3_y = rep(0, 10))

V_x_mat_susie <- data.frame(gamma1_x_susie = rep(0, 10),
                            gamma2_x_susie = rep(0, 10),
                            gamma3_x_susie = rep(0, 10))

V_y_mat_susie <- data.frame(gamma1_y_susie = rep(0, 10),
                            gamma2_y_susie = rep(0, 10),
                            gamma3_y_susie = rep(0, 10))

b_post_mat <- data.frame(gamma1_x = rep(0, num_total_vars),
                         gamma2_x = rep(0, num_total_vars),
                         gamma3_x = rep(0, num_total_vars))

a_post_mat <- data.frame(gamma1_y = rep(0, num_total_vars),
                         gamma2_y = rep(0, num_total_vars),
                         gamma3_y = rep(0, num_total_vars))

b_post_mat_susie <- data.frame(gamma1_x_susie = rep(0, num_total_vars),
                               gamma2_x_susie = rep(0, num_total_vars),
                               gamma3_x_susie = rep(0, num_total_vars))

a_post_mat_susie <- data.frame(gamma1_y_susie = rep(0, num_total_vars),
                               gamma2_y_susie = rep(0, num_total_vars),
                               gamma3_y_susie = rep(0, num_total_vars))

N_exposure <- as.numeric(argv$N_exposure)
N_outcome <- as.numeric(argv$N_outcome)

#Get the desired number of signals
num_signals <- as.numeric(argv$num_causal)
num_pleio <- as.numeric(argv$num_pleio)

for(i in 1:length(gamma_vec)){
  gamma_true_curr <- gamma_vec[i]
  
  af_curr <- df_af$af_exposure
  
  #Specify causal relationship 
  G <- matrix(c(0, gamma_true_curr,
                0, 0), 2, 2, byrow = TRUE)
  
  colnames(G) <- row.names(G) <- c("X", "Y")
  
  ##############################################################################
  N <- c(N_exposure, N_outcome)
  set.seed(as.numeric(argv$run))
  idx_X_causal <- sample(1:num_total_vars, num_signals)
  
  if(num_pleio > 0){
    set.seed(as.numeric(argv$run))
    idx_Y_causal <- sample(c(1:num_total_vars)[-idx_X_causal], num_pleio) #Save this as well as the LD between these SNPs (probably easier to save true LD first) 
  }
  
  #Order them for LD matrix output
  idx_X_causal <- idx_X_causal[order(idx_X_causal)]
  if(num_pleio > 0){
    idx_Y_causal <- idx_Y_causal[order(idx_Y_causal)]
  }

  pi_mat <- matrix(0, nrow = num_total_vars, ncol = 2)
  pi_mat[idx_X_causal, 1] <- 1
  if(num_pleio > 0){
    pi_mat[idx_Y_causal, 2] <- 1  
  }
  
  
  #Save LD matrix for exposure causal SNPs
  LD_X <- as.matrix(ld_curr[idx_X_causal, idx_X_causal])
  rownames(LD_X) <- colnames(LD_X) <- idx_X_causal
  
  if(num_pleio > 0){
    idx_full_Y_causal <- c(idx_X_causal, idx_Y_causal) 
  }else{
    idx_full_Y_causal <- idx_X_causal
  }
  
  LD_Y <- as.matrix(ld_curr[idx_full_Y_causal, idx_full_Y_causal])
  rownames(LD_Y) <- colnames(LD_Y) <- idx_full_Y_causal #Diagonals are not all 1 ?? Might have to regenerate?; Also I projected to nearPD after simulating genotypes so won't be exactly the same as 
  
  #Simulate summary statistics
  set.seed(as.numeric(argv$run))
  N_vec <- c(N_exposure, N_outcome)
  #Want to not decrease pleiotropy as gamma increases
  if(num_pleio > 0){
    herit_y <- 0.01 + gamma_true_curr^2*0.05 #May or may not be correct 
  }else{
    herit_y <- gamma_true_curr^2*0.05
  }
  
  #herit_y <- 0.01/(1-gamma_true_curr^2*0.05) ### This is wrong right? - was herit_y <- 0.01/(1-gamma_true_curr^2) before
  sim_dat <- sim_mv(G = G,
                    J = num_total_vars,
                    N = N_vec, 
                    h2 = c(0.05, herit_y), 
                    pi = pi_mat, 
                    R_LD = list(ld_curr), 
                    af = af_curr,
                    est_s = TRUE) #Is sporadic pleiotropy on?

  
  pruned_set <- sim_ld_prune(sim_dat, 
                             pvalue = 1, 
                             R_LD = list(ld_curr), 
                             r2_thresh = 0.01,
                             pval_thresh = 1e-6) ###
  
  gwasResultsX <- data.frame(CHR = 19, BP = 1:nrow(sim_dat$beta_hat), 
                             SNP = paste0("snp", 1:nrow(sim_dat$beta_hat)), 
                             P = 2*pnorm(-abs(sim_dat$beta_hat[,1]/sim_dat$s_estimate[,1])),
                             Z = sim_dat$beta_hat[,1]/sim_dat$s_estimate[,1])
  
  gwasResultsY <- data.frame(CHR = 19, BP = 1:nrow(sim_dat$beta_hat),
                             SNP = paste0("snp", 1:nrow(sim_dat$beta_hat)), 
                             P = 2*pnorm(-abs(sim_dat$beta_hat[,2]/sim_dat$s_estimate[,2])),
                             Z = sim_dat$beta_hat[,2]/sim_dat$s_estimate[,2])
  
  #Keep gwasResultsX and gwasResultsY as is to ensure the below code works, but create
  #two new files for use in MR-link2
  gwasResultsX_mrlink2 <- data.frame(pos_name = df_af$vars_hg38,
                                     chromosome = unlist(strsplit(df_af$vars_hg38, split = ":", fixed = TRUE))[c(TRUE, FALSE, FALSE, FALSE)],
                                     position = var_pos,
                                     effect_allele = unlist(strsplit(df_af$vars_hg38, split = ":", fixed = TRUE))[c(FALSE, FALSE, TRUE, FALSE)],
                                     reference_allele = unlist(strsplit(df_af$vars_hg38, split = ":", fixed = TRUE))[c(FALSE, FALSE, FALSE, TRUE)],
                                     beta = sim_dat$beta_hat[,1],
                                     se = sim_dat$s_estimate[,1],
                                     z = sim_dat$beta_hat[,1]/sim_dat$s_estimate[,1],
                                     pval = 2*pnorm(-abs(sim_dat$beta_hat[,1]/sim_dat$s_estimate[,1])),
                                     n_iids = N_exposure)
  
  gwasResultsY_mrlink2 <- data.frame(pos_name = df_af$vars_hg38,
                                     chromosome = unlist(strsplit(df_af$vars_hg38, split = ":", fixed = TRUE))[c(TRUE, FALSE, FALSE, FALSE)],
                                     position = var_pos,
                                     effect_allele = unlist(strsplit(df_af$vars_hg38, split = ":", fixed = TRUE))[c(FALSE, FALSE, TRUE, FALSE)],
                                     reference_allele = unlist(strsplit(df_af$vars_hg38, split = ":", fixed = TRUE))[c(FALSE, FALSE, FALSE, TRUE)],
                                     beta = sim_dat$beta_hat[,2],
                                     se = sim_dat$s_estimate[,2],
                                     z = sim_dat$beta_hat[,2]/sim_dat$s_estimate[,2],
                                     pval = 2*pnorm(-abs(sim_dat$beta_hat[,2]/sim_dat$s_estimate[,2])),
                                     n_iids = N_outcome) 
  
  #Write out the MR-link2 files
  path_exp <- argv$output_MRlink2_exp_path
  file_exp <- str_glue("{path_exp}gam0{i-1}_exposureSize_{num_signals}_outcomeSize_{num_pleio}_simNum_{argv$run}.MRLink_exp")
  fwrite(gwasResultsX_mrlink2, file_exp, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  system(sprintf("bgzip -f %s", file_exp))

  path_out <- argv$output_MRlink2_out_path
  file_out <- str_glue("{path_exp}gam0{i-1}_exposureSize_{num_signals}_outcomeSize_{num_pleio}_simNum_{argv$run}.MRLink_out")
  fwrite(gwasResultsY_mrlink2, file_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  system(sprintf("bgzip -f %s", file_out))
  
  #Now the COJO formatted files to run GCTA-COJO
  gwasResultsX_cojo <- data.frame(SNP = df_af$vars_hg38,
                                  A1 = unlist(strsplit(df_af$vars_hg38, split = ":", fixed = TRUE))[c(FALSE, FALSE, TRUE, FALSE)],
                                  A2 = unlist(strsplit(df_af$vars_hg38, split = ":", fixed = TRUE))[c(FALSE, FALSE, FALSE, TRUE)],
                                  freq = af_curr,
                                  b = sim_dat$beta_hat[,1],
                                  se = sim_dat$s_estimate[,1],
                                  p = 2*pnorm(-abs(sim_dat$beta_hat[,1]/sim_dat$s_estimate[,1])),
                                  N = N_exposure)
  
  
  gwasResultsY_cojo <- data.frame(SNP = df_af$vars_hg38,
                                  A1 = unlist(strsplit(df_af$vars_hg38, split = ":", fixed = TRUE))[c(FALSE, FALSE, TRUE, FALSE)],
                                  A2 = unlist(strsplit(df_af$vars_hg38, split = ":", fixed = TRUE))[c(FALSE, FALSE, FALSE, TRUE)],
                                  freq = af_curr,
                                  b = sim_dat$beta_hat[,2],
                                  se = sim_dat$s_estimate[,2],
                                  p = 2*pnorm(-abs(sim_dat$beta_hat[,2]/sim_dat$s_estimate[,2])),
                                  N = N_outcome)
  
  #Output the COJO results
  path_exp <- argv$output_cojo_exp_path
  file_exp <- str_glue("{path_exp}gam0{i-1}_exposureSize_{num_signals}_outcomeSize_{num_pleio}_simNum_{argv$run}.cojo_exp.tsv")
  fwrite(gwasResultsX_cojo, file_exp, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  system(sprintf("bgzip -f %s", file_exp))
  
  path_out <- argv$output_MRlink2_out_path
  file_out <- str_glue("{path_out}gam0{i-1}_exposureSize_{num_signals}_outcomeSize_{num_pleio}_simNum_{argv$run}.cojo_out.tsv")
  fwrite(gwasResultsY_cojo, file_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  system(sprintf("bgzip -f %s", file_out))
  
  
  #We don't have in-sample allele frequencies for the study so use the ones we used to generate
  gwasResultsX$eaf <- gwasResultsY$eaf <- af_curr
  
  #Use TwoSampleMR package
  #This package requires at a minimum SNP id, beta, se, and effect_allele so make up the 
  #effect allele (I think this is only to match the effect size in the outcome study so 
  #if they are the same it shouldn't cause any issues)
  gwasResultsX$beta <- sim_dat$beta_hat[,1]
  gwasResultsX$se <- sim_dat$s_estimate[,1]
  gwasResultsX$effect_allele <- "A"
  
  gwasResultsY$beta <- sim_dat$beta_hat[,2]
  gwasResultsY$se <- sim_dat$s_estimate[,2]
  gwasResultsY$effect_allele <- "A"
  
  #Keep only the LD pruned set of variants
  gwasResultsX_MR <- gwasResultsX[pruned_set,]
  gwasResultsY_MR <- gwasResultsY[pruned_set,]
  
  if(length(pruned_set) > 0){
    df_exposure <- format_data(gwasResultsX_MR, type = "exposure")
    df_outcome <- format_data(gwasResultsY_MR, type = "outcome")
    
    #Put action = 1 to make it ignore palindromic SNPs
    dat_mr <- harmonise_data(
      exposure_dat = df_exposure, 
      outcome_dat = df_outcome,
      action = 1
    )
    
    if(sum(dat_mr$mr_keep) > 0){
      if(dim(dat_mr)[1] > 1){
        df_mr <- mr(dat_mr, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
        df_raps <- mr_modified(dat_mr, method_list = c("mr_raps"), 
                               parameters = list(over.dispersion = TRUE, loss.function = "huber"))
        
        df_final <- rbind(df_mr, df_raps)
        
      }else{
        df_mr <- mr_singlesnp(dat_mr)
        df_mr$nsnp <- 1
        df_mr$method = "Wald ratio"
        names(df_mr)[names(df_mr) == "p"] <- "pval"
        df_mr <- subset(df_mr, select = c("id.exposure", "id.outcome", "outcome",
                                          "exposure", "method", "nsnp", "b", "se",
                                          "pval"))
        df_raps <- mr_modified(dat_mr, method_list = c("mr_raps"),
                               parameters = list(over.dispersion = TRUE, loss.function = "huber"))
        
        df_final <- rbind(df_mr[1,], df_raps)
      }
      
      df_final <- df_final[,5:9]
      
      #Perform SMR
      max_snp <- which.min(gwasResultsX_MR$P)
      test_stat <- (gwasResultsX_MR$Z[max_snp]^2 * gwasResultsY_MR$Z[max_snp]^2)/(gwasResultsX_MR$Z[max_snp]^2 + gwasResultsY_MR$Z[max_snp]^2)
      p_smr <- pchisq(test_stat, 1, lower.tail = FALSE)
      
      smr_est <- gwasResultsY_MR$beta[max_snp]/gwasResultsX_MR$beta[max_snp]
      smr_est_var <- ((gwasResultsY_MR$beta[max_snp]^2)/(gwasResultsX_MR$beta[max_snp]^2)) * (gwasResultsY_MR$se[max_snp]^2/gwasResultsY_MR$beta[max_snp]^2 + gwasResultsX_MR$se[max_snp]^2/gwasResultsX_MR$beta[max_snp]^2)
      
      df_final <- rbind(df_final, data.frame(method = "SMR",
                                             nsnp = 1, 
                                             b = smr_est,
                                             se = sqrt(smr_est_var),
                                             pval = p_smr))
      
      df_final$gamma <- gamma_true_curr
      
    }else{
      df_final <- data.frame(method = NA, 
                             nsnp = 0,
                             b = NA,
                             se = NA,
                             pval = NA, 
                             gamma = gamma_true_curr)
    }
  }else{
    df_final <- data.frame(method = NA, 
                           nsnp = 0,
                           b = NA,
                           se = NA,
                           pval = NA, 
                           gamma = gamma_true_curr)
  } ## Make sure to update Egger
  
  p_x <- 2*pnorm(-abs(sim_dat$beta_hat[,1]/sim_dat$s_estimate[,1]))
  p_y <- 2*pnorm(-abs(sim_dat$beta_hat[,2]/sim_dat$s_estimate[,2]))
  p_min <- c(p_min, min(p_x))
  
  #Add single-variant statistics to our matrix
  df_temp <- cbind(sim_dat$beta_joint, sim_dat$beta_hat, sim_dat$s_estimate,
                   data.frame(px = p_x, py = p_y))
  names(df_temp) <- paste(names(df_temp), i, sep = "_")
  
  if(exists("df_singleVar")){
    df_singleVar <- cbind(df_singleVar, df_temp)
  }else{
    df_singleVar <- df_temp
  }
  
  #Get sufficient statistics
  exposure_xtx <- exposure_xty <-  list()
  
  exposure_xtx[[1]] <- (N_exposure - 1)*ld_curr
  exposure_xty[[1]] <- sqrt(N_exposure - 1)*(sim_dat$beta_hat[,1]/sim_dat$s_estimate[,1]) ## Check this for proper scaling
  exposure_yty <- N_exposure - 1
  
  outcome_xtx <- outcome_xty <- list()
  
  outcome_xtx[[1]] <- (N_outcome - 1)*ld_curr
  outcome_xty[[1]] <- sqrt(N_outcome - 1)*(sim_dat$beta_hat[,2]/sim_dat$s_estimate[,2])
  outcome_yty <- N_outcome - 1 
  
  #Now get SuSiE output
  susie_fit_x <- susieR::susie_suff_stat(XtX = exposure_xtx[[1]], 
                                         Xty = exposure_xty[[1]], 
                                         yty = exposure_yty,
                                         n = N_exposure,
                                         estimate_residual_variance = FALSE)
  
  susie_fit_y <- susieR::susie_suff_stat(XtX = outcome_xtx[[1]], 
                                         Xty = outcome_xty[[1]], 
                                         yty = outcome_yty,
                                         n = N_outcome,
                                         estimate_residual_variance = FALSE)
  
  V_x_mat_susie[[i]] <- susie_fit_x$V
  V_y_mat_susie[[i]] <- susie_fit_y$V
  
  b_post_susie <- susie_get_posterior_mean(susie_fit_x)
  a_post_susie <- susie_get_posterior_mean(susie_fit_y)
  
  b_post_mat_susie[[i]] <- b_post_susie
  a_post_mat_susie[[i]] <- a_post_susie
  
  ###############################################################################
  #Now run GSMR
  #First need to pick out the variants to use
  #From each credible set, take highest PIP var
  num_cs <- length(susie_fit_x$sets$cs)
  
  #Note also that the credible sets are not necessarily in order 
  #Need to gsub out the L and get the index
  idx_sets <- as.numeric(gsub("L", "", names(susie_fit_x$sets$cs)))
  
  #Keep track of which variants to include
  idx_to_use <- c()
  for(j in 1:num_cs){
    set_curr <- idx_sets[j]
    idx_curr <- susie_fit_x$sets$cs[[j]]
    pip_curr <- susie_fit_x$alpha[set_curr, idx_curr]
    idx_to_include <- idx_curr[which.max(pip_curr)]
    idx_to_use <- c(idx_to_use, idx_to_include)
  }
  
  #One iteration had the same variant as the top var in two CS so keep only unique
  idx_to_use <- unique(idx_to_use)
  
  ldrho <- ld_curr[idx_to_use, idx_to_use, drop = FALSE] ## Error for GUSB iter 21 numcausal = 5, num_pleio = 2
  
  df_exposure_gsmr <- gwasResultsX[idx_to_use,] 
  df_outcome_gsmr <- gwasResultsY[idx_to_use,] 
  
  #Get the SNPs corresponding to these
  snp_names <- df_exposure_gsmr$SNP
  
  colnames(ldrho) <- rownames(ldrho) <- snp_names
  
  #Now back-transform the Z-scores (not sure about doing this here?)
  #z_x_orig <- df_exposure_gsmr$Z
  #z_y_orig <- df_outcome_gsmr$Z
  
  #af_x_orig <- df_exposure_gsmr$eaf
  #af_y_orig <- df_outcome_gsmr$eaf
  
  #ss_x_orig <- N_exposure
  #ss_y_orig <- N_outcome
  
  #b_hat_exposure <- z_x_orig/(sqrt(2*af_x_orig*(1-af_x_orig)*(ss_x_orig + z_x_orig^2)))
  #se_exposure <- 1/(sqrt(2*af_x_orig*(1-af_x_orig)*(ss_x_orig + z_x_orig^2)))
  
  #b_hat_outcome <- z_y_orig/(sqrt(2*af_y_orig*(1-af_y_orig)*(ss_y_orig + z_y_orig^2)))
  #se_outcome <- 1/(sqrt(2*af_y_orig*(1-af_y_orig)*(ss_y_orig + z_y_orig^2)))
  
  #Z_adj_exposure <- b_hat_exposure/se_exposure
  #Z_adj_outcome <- b_hat_outcome/se_outcome 
  
  print(num_cs)
  
  #If more than 1 SNP, run GSMR
  if(num_cs > 1){
    #Now set the parameters for GSMR
    bzx <- df_exposure_gsmr$beta
    bzx_se <- df_exposure_gsmr$se
    bzx_pval <- df_exposure_gsmr$P
    
    bzy <- df_outcome_gsmr$beta
    bzy_se <- df_outcome_gsmr$se
    bzy_pval <- df_outcome_gsmr$P
    
    n_ref <- N_exposure
    gwas_thresh <- 1
    single_snp_heidi_thresh <- 1
    multi_snps_heidi_thresh <- 1
    nsnps_thresh <- 1
    heidi_outlier_flag <- F
    ld_r2_thresh <- 1
    ld_fdr_thresh <- 1.01
    gsmr2_beta <- 0
    
    snp_coeff_id <- snp_names
    
    gsmr_results <-  gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh,
                          single_snp_heidi_thresh, multi_snps_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta)
    
    num_snps_gsmr <- length(gsmr_results$used_index)
    
  }else if(num_cs == 1){
    #Run SMR instead if only 1 SNP (maybe change to adj later?)
    test_stat <- (df_exposure_gsmr$Z^2 * df_outcome_gsmr$Z^2)/(df_exposure_gsmr$Z^2 + df_outcome_gsmr$Z^2)
    p_smr <- pchisq(test_stat, 1, lower.tail = FALSE)
    
    smr_est <- df_outcome_gsmr$beta/df_exposure_gsmr$beta
    smr_est_var <- ((df_outcome_gsmr$beta^2)/(df_exposure_gsmr$beta^2)) * (df_outcome_gsmr$se^2/df_outcome_gsmr$beta^2 + df_exposure_gsmr$se^2/df_exposure_gsmr$beta^2)
    
    gsmr_results <-  list()
    
    gsmr_results$bxy <- smr_est
    gsmr_results$bxy_se <- sqrt(smr_est_var)
    gsmr_results$bxy_pval <- p_smr
    num_snps_gsmr <- 1
  }else{
    #There are no credible sets so return NA
    gsmr_results <- list()
    
    gsmr_results$bxy <- NA
    gsmr_results$bxy_se <- NA
    gsmr_results$bxy_pval <- NA
    num_snps_gsmr <- 0
  } 
  
  #Now add the GSMR results to the traditional method data frame
  df_gsmr <- data.frame(method = c("GSMR_pack"), 
                        nsnp = num_snps_gsmr,
                        b = c(gsmr_results$bxy),
                        se = c(gsmr_results$bxy_se),
                        pval =  c(gsmr_results$bxy_pval), 
                        gamma = gamma_true_curr)
  
  df_final <- rbind(df_final, df_gsmr)
  trad_mr_list[[i]] <- df_final
  
  #susie_fit_x
  #susie_fit_y
  
  #Try to modify susie fits to have large variance
  #idx_nonZero_Vx <- which(susie_fit_x$V != 0)
  #susie_fit_x$mu2[idx_nonZero_Vx,] <- 1
  #susie_fit_x$V[idx_nonZero_Vx] <- 1
  
  #idx_nonZero_Vy <- which(susie_fit_y$V != 0)
  #susie_fit_y$V[idx_nonZero_Vy] <- 1
  #susie_fit_y$mu2[idx_nonZero_Vy,] <- 1
  
  #Try manually setting it to something large
  #susie_init <- susie_init_y  <- list(V = rep(1, 10),
  #                                    mu = matrix(0, nrow = 10, ncol = nrow(outcome_xtx[[1]])),
  #                                    mu2 = matrix(1, nrow = 10, ncol = nrow(outcome_xtx[[1]])),
  #                                    alpha = matrix(0, nrow = 10, ncol = nrow(outcome_xtx[[1]])))
  

  #Now call our finemappingMR function
  res_ourMeth <- tryCatch(expr = finemappingMR_suff_stats(exposure_xtx, exposure_xty, exposure_yty,
                                                     outcome_xtx, outcome_xty, outcome_yty,
                                                     N_exposure, N_outcome,
                                                     L_x = 10, L_y = 10,
                                                     scaled_prior_variance_x = 0.2,
                                                     scaled_prior_variance_y = 0.2,
                                                     estimate_prior_variance_x = TRUE,
                                                     estimate_prior_variance_y = TRUE,
                                                     residual_variance_x = NULL,
                                                     residual_variance_y = NULL,
                                                     estimate_residual_variance_x = FALSE,
                                                     estimate_residual_variance_y = FALSE,
                                                     tol = 1e-3,
                                                     max_iter = 1000,
                                                     calc_cs_x = TRUE,
                                                     calc_cs_y = TRUE,
                                                     verbose = TRUE,
                                                     susie_init = NULL,
                                                     susie_init_y = NULL,
                                                     mu_gamma_init = 0,
                                                     mu2_gamma_init = 0,
                                                     sigma2_gamma_prior = 10^6,
                                                     num_samples = 1000,
                                                     beta_gamma_alpha = FALSE), 
                          
                          error = function(cond){
                            to_return <- list()
                            to_return$res <- data.frame(gamma = NA,
                                                        gamma_var = NA, 
                                                        gamma_total_var = NA,
                                                        iter = NA) #Note that before gamma_total_var didn't exist here and caused issues so check the others
                            to_return$elbo <- NA
                            to_return$V_x <- matrix(rep(0, 10), nrow = 10, ncol = 1)
                            to_return$V_y <- matrix(rep(0, 10), nrow = 10, ncol = 1)
                            to_return$alpha_a <- to_return$alpha_b <- to_return$mu_a <- to_return$mu_b <- list(matrix(0, nrow = 10, ncol = num_total_vars))
                            
                            return(to_return)
                          })
  
  #Store posterior estimates of gamma
  gamma_est <- c(gamma_est, res_ourMeth$res$gamma[1])
  gamma_var <- c(gamma_var, res_ourMeth$res$gamma_var[1])
  gamma_var_total_var <- c(gamma_var_total_var, res_ourMeth$res$gamma_total_var[1])
  iter_num <- c(iter_num, res_ourMeth$res$iter[1])
  elbo_list[[i]] <- res_ourMeth$elbo
  
  #Store final prior variance estimates to compare to SuSiE estimates
  V_x_mat[[i]] <- res_ourMeth$V_x
  V_y_mat[[i]] <- res_ourMeth$V_y
  
  #Store final posterior mean estimates to compare to SuSiE estimates
  prior_tol <- 1e-9
  include_idx_x <- which(res_ourMeth$V_x[,1] > prior_tol)
  b_post <- colSums((res_ourMeth$alpha_b[[1]]*res_ourMeth$mu_b[[1]])[include_idx_x,,drop = FALSE])
  
  include_idx_y <- which(res_ourMeth$V_y[,1] > prior_tol)
  a_post <- colSums((res_ourMeth$alpha_a[[1]]*res_ourMeth$mu_a[[1]])[include_idx_y,,drop = FALSE])
  
  b_post_mat[[i]] <- b_post
  a_post_mat[[i]] <- a_post
  
}


#Save the expected LD
out_file <- sub(".gz", "", argv$output_expectedLD, fixed = T)
fwrite(LD_Y, file = out_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE) #This only saves the final LD matrix
system(sprintf("bgzip -f %s", out_file))

#Write single variant results
out_file <- sub(".gz", "", argv$output_singleVar, fixed = T)
fwrite(df_singleVar, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
system(sprintf("bgzip -f %s", out_file))

#Create data frame to output 
output_df <- data.frame(gamma_true = c(0, 0.1, 0.2), 
                        gamma_est = gamma_est, 
                        gamma_var = gamma_var,
                        gamma_var_total_var = gamma_var_total_var,
                        iter = iter_num)

out_file <- sub(".gz", "", argv$output_res, fixed = T)
fwrite(output_df, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
system(sprintf("bgzip -f %s", out_file))

#Write ELBO used for convergence
max_elbo <- max(sapply(elbo_list, length))
for(i in 1:length(elbo_list)){
  len_elbo <- length(elbo_list[[i]])
  if(len_elbo < max_elbo){
    new_elbo <- c(elbo_list[[i]], rep(elbo_list[[i]][len_elbo], max_elbo - len_elbo))
    elbo_list[[i]] <- new_elbo
  }
}

output_elbo <- data.frame(run = 1:max_elbo, 
                          ELBO1 = elbo_list[[1]], 
                          ELBO2 = elbo_list[[2]], 
                          ELBO3 = elbo_list[[3]])

out_file <- sub(".gz", "", argv$output_elboConv, fixed = T)
fwrite(output_elbo, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
system(sprintf("bgzip -f %s", out_file))

#Write trad MR data frame - need to fix this to output every estimate (think it does now?)
trad_mr_df <- trad_mr_list[[1]]
for(i in 2:length(trad_mr_list)){
  trad_mr_df <- rbind(trad_mr_df, trad_mr_list[[i]])
}

out_file <- sub(".gz", "", argv$output_tradMR, fixed = T)
fwrite(trad_mr_df, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
system(sprintf("bgzip -f %s", out_file))

#Write prior variance estimates
V_comb <- cbind(V_x_mat, V_y_mat, V_x_mat_susie, V_y_mat_susie)

out_file <- sub(".gz", "", argv$output_priorVar, fixed = T)
fwrite(V_comb, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
system(sprintf("bgzip -f %s", out_file))

#Write posterior mean estimates
post_comb <- cbind(b_post_mat, a_post_mat, b_post_mat_susie, a_post_mat_susie)

out_file <- sub(".gz", "", argv$output_posteriorMean, fixed = T)
fwrite(post_comb, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
system(sprintf("bgzip -f %s", out_file))

