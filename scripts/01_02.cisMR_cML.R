library("argparser")
library("dplyr")
library("data.table")
library("stringr")
library("cisMRcML")

p <- arg_parser("Script that runs cisMR-cML on our pleiotropic simulation")

p <- add_argument(p, "--gwasX_path", help = "Path to sumstats for exposure")
p <- add_argument(p, "--gwasY_path", help = "Path to sumstats for outcome")
p <- add_argument(p, "--cojo_exp_path", help = "Path to the GCTA-COJO results for the exposure")
p <- add_argument(p, "--cojo_out_path", help = "Path to the GCTA-COJO results for the outcome")
p <- add_argument(p, "--ld_file", help = "Name of .rds file containing the LD")
p <- add_argument(p, "--af_file", help = "Name of file containing the allele frequencies of the variants")
p <- add_argument(p, "--numCausalX", help = "Number of causal variants for the exposure")
p <- add_argument(p, "--numCausalY", help = "Number of causal variants for the outcome (pleiotropic)")
p <- add_argument(p, "--N_exposure", help = "Size of the exposure")
p <- add_argument(p, "--run", help = "Seed used to generate data")

p <- add_argument(p, "--output", help = "Name of file to output to")

argv <- parse_args(p)

numCausalX <- argv$numCausalX
numCausalY <- argv$numCausalY
run <- argv$run
N_exposure <- as.numeric(argv$N_exposure)

#Read in the LD
R_full <- readRDS(argv$ld_file)

#Read in allele frequencies to match variants in the cojo file to the LD
df_af <- fread(argv$af_file, header = TRUE, data.table = FALSE)
var_pos <- as.numeric(unlist(strsplit(df_af$vars_hg38, split = ":", fixed = TRUE))[c(FALSE, TRUE, FALSE, FALSE)])

gamma_true <- c(0, 0.1, 0.2)
df_final <- data.frame()

for(i in 1:3){
  print(i)
  gamma_curr <- gamma_true[i]
  
  #Read in the COJO files
  file_exposure <- str_glue("{argv$cojo_exp_path}gam0{i-1}_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.jma.cojo")
  df_exp_cojo <- fread(file_exposure, header = TRUE, data.table = FALSE)
  
  file_outcome <- str_glue("{argv$cojo_out_path}gam0{i-1}_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.out.jma.cojo")
  df_out_cojo <- fread(file_outcome, header = TRUE, data.table = FALSE)
  
  #Get the jointly associated SNPs
  if(nrow(df_exp_cojo) > 0){
    idx_x <- which(var_pos %in% df_exp_cojo$bp)
  }else{
    idx_x <- c()
  }
  
  if(nrow(df_out_cojo) > 0){
    idx_y <- which(var_pos %in% df_out_cojo$bp)
  }else{
    idx_y <- c()
  }
  
  idx_full <- unique(c(idx_x, idx_y))
  
  if(length(idx_full) != 0){
    R_red <- R_full[idx_full, idx_full]
    
    ########################################
    #Now read in the exposure and outcome summary stats
    gwasX_file <- str_glue("{argv$gwasX_path}gam0{i-1}_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.cojo_exp.tsv.gz")
    df_exp <- fread(cmd = str_glue("zcat {gwasX_file}"), header = TRUE, data.table = FALSE)
    
    gwasY_file <- str_glue("{argv$gwasY_path}gam0{i-1}_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.cojo_out.tsv.gz")
    df_out <- fread(cmd = str_glue("zcat {gwasY_file}"), header = TRUE, data.table = FALSE)
    
    df_exp <- df_exp[idx_full,]
    df_out <- df_out[idx_full,]
    
    #Get conditional betas according to their vignette
    cor_exp <- df_exp$b / sqrt(df_exp$b^2 + (df_exp$N - 2) *
                                 df_exp$se^2)
    
    df_exp$bJ <- solve(R_red) %*% cor_exp
    
    
    cor_out <- df_out$b / sqrt(df_out$b^2 + (df_out$N - 2) *
                                 df_out$se^2)
    
    df_out$bJ <- solve(R_red) %*% cor_out
    
    #Get adjusted ses
    se_exp_cor <- df_exp$se / sqrt(df_exp$b^2 + (df_exp$N - 2) * df_exp$se^2)
    se_out_cor <- df_out$se / sqrt(df_out$b^2 + (df_out$N - 2) * df_out$se^2)
    
    #Run cisMR-cML
    b_exp_cond <- df_exp$bJ
    b_out_cond <- df_out$bJ 
    
    Sig_exp1 <- solve(R_red) %*% (se_exp_cor %o% se_exp_cor * R_red) %*% solve(R_red)
    
    Sig_out1 <- solve(R_red) %*% (se_out_cor %o% se_out_cor * R_red) %*% solve(R_red)
    
    Sig_exp_inv <- solve(Sig_exp1)
    Sig_out_inv <- solve(Sig_out1)
    
    ciscML_res <- tryCatch(expr = {cismr_cML_DP(b_exp = b_exp_cond,
                                                b_out = b_out_cond,
                                                Sig_exp_inv = Sig_exp_inv,
                                                Sig_out_inv = Sig_out_inv,
                                                maxit = 200,
                                                n = N_exposure, 
                                                random_start = 5,
                                                min_theta_range = -0.2,
                                                max_theta_range = 0.2, #min_theta - max_theta = -0.2:0.2 worked well but will do this instead as it makes more sense
                                                num_pert = 100,
                                                random_start_pert = 5,
                                                random_seed = 12345)},
                           
                           error = function(cond){
                             to_return <- list()
                             to_return$BIC_DP_theta <- NA
                             to_return$BIC_DP_se <- NA
                             to_return$BIC_DP_p <- NA
                             return(to_return)
                           } #Just returning errors now since I am probably not going to debug their code
    )
    
    #Make a data frame
    df_curr <- data.frame(gamma_true = gamma_true[i],
                          gamma_est = ciscML_res$BIC_DP_theta,
                          gamma_se = ciscML_res$BIC_DP_se,
                          gamma_p = ciscML_res$BIC_DP_p,
                          num_exp = length(idx_x),
                          num_out = length(idx_y))
    
  }else{
    df_curr <- data.frame(gamma_true = gamma_true[i],
                          gamma_est = NA,
                          gamma_se = NA,
                          gamma_p = NA,
                          num_exp = length(idx_x),
                          num_out = length(idx_y))
  }
  
  df_final <- rbind(df_final, df_curr)
  
}

#Write out the table
out_file <- sub(".gz", "", argv$output, fixed = T)
fwrite(df_final, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "NA")
system(sprintf("bgzip -f %s", out_file))

