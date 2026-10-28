# Snakemake pipeline for the finemappingMR manuscript simulations.
# Simulation 1: resampled LD from 5 proteins (GUSB, SERPINC1, F2, RBKS, BCAM), 5 random causal
# variants for the exposure and either 0 or 2 pleiotropic causal variants for the outcome.

import pandas as pd
import snakemake as _smk

REQUIRED_SNAKEMAKE = "7.12.1"
if _smk.__version__ != REQUIRED_SNAKEMAKE:
    raise WorkflowError(
        f"This workflow requires Snakemake {REQUIRED_SNAKEMAKE} (found {_smk.__version__}). "
        "Activate the environment with: conda activate finemappingMR_manuscript"
    )

gene_ranges = pd.read_csv("config/gene_ranges.tsv", sep="\t")

# gene -> "Chr:Start-End"
gene_to_region = {
    row["Gene"]: f'{row["Chr"]}:{row["Start"]}-{row["End"]}'
    for _, row in gene_ranges.iterrows()
}

MRLINK2 = "external/mrlink2/mr_link_2_standalone.py"

RUN = range(1,201)
GENES = ["GUSB", "SERPINC1", "F2", "RBKS", "BCAM"]
NUM_CAUSAL_X = [5]
NUM_CAUSAL_Y = [0, 2]
N_EXPOSURE = [5000]
N_OUTCOME = [5000]

rule all:
    input:
        expand("01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.cisMRcML.res.tsv.gz", gene = GENES, N_exposure = N_EXPOSURE, N_outcome = N_OUTCOME, numCausalX = NUM_CAUSAL_X, numCausalY = NUM_CAUSAL_Y, run = RUN)

#Rule to run finemappingMR and generate the summary statistics 
rule sim_finemappingMR_resampledLD:
    input:
        ld = "00.inputs/{gene}/{gene}_LD_resampled.rds",
        af = "00.inputs/{gene}/{gene}_af_resampled.tsv"
    output:
        res = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.BayesianRes.gz",
        elbo_conv = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.ELBORes.gz",
        trad_MR = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.tradMR.gz",
        singleVar = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.singleVar.gz",
        priorVar = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.priorVar.gz", 
        posteriorMean = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.posteriorMean.gz",
        expectedLD = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.expectedLD.gz",
        MRLink_exp_gam0 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam00_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.MRLink_exp.gz",
        MRLink_out_gam0 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam00_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.MRLink_out.gz",
        MRLink_exp_gam01 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam01_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.MRLink_exp.gz",
        MRLink_out_gam01 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam01_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.MRLink_out.gz",
        MRLink_exp_gam02 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam02_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.MRLink_exp.gz",
        MRLink_out_gam02 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam02_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.MRLink_out.gz",
        cojo_exp_gam0 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam00_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.cojo_exp.tsv.gz",
        cojo_out_gam0 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam00_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.cojo_out.tsv.gz",
        cojo_exp_gam01 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam01_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.cojo_exp.tsv.gz",
        cojo_out_gam01 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam01_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.cojo_out.tsv.gz",
        cojo_exp_gam02 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam02_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.cojo_exp.tsv.gz",
        cojo_out_gam02 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam02_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.cojo_out.tsv.gz"
    params:
        MRLink_exp_path = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/",
        MRLink_out_path = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/",
        cojo_exp_path = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/",
        cojo_out_path = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/"
    shell:
        "Rscript scripts/01_01.sim_randomCausal.R --run {wildcards.run} --N_exposure {wildcards.N_exposure} --N_outcome {wildcards.N_outcome} \
        --output_res {output.res} --output_elboConv {output.elbo_conv} --output_tradMR {output.trad_MR} --output_singleVar {output.singleVar} \
        --output_priorVar {output.priorVar} --output_posteriorMean {output.posteriorMean} --output_MRlink2_exp_path {params.MRLink_exp_path}  --output_MRlink2_out_path {params.MRLink_out_path} \
        --output_expectedLD {output.expectedLD} --num_causal {wildcards.numCausalX} --num_pleio {wildcards.numCausalY} \
        --ld_file {input.ld} --af_file {input.af} \
        --output_cojo_exp_path {params.cojo_exp_path}  --output_cojo_out_path {params.cojo_out_path}"
        
#Rule to run MR-link-2
rule python_mr_link2:
    input:
        "00.inputs/{gene}/{gene}_resampled_hg38_genio.bed", 
        "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam00_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.MRLink_exp.gz",
        "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam01_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.MRLink_exp.gz",
        "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam02_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.MRLink_exp.gz",
        "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam00_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.MRLink_out.gz",
        "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam01_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.MRLink_out.gz",
        "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam02_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.MRLink_out.gz"
    output:
        "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam00_MRLink.res.txt",
        "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam01_MRLink.res.txt",
        "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam02_MRLink.res.txt"
    params:
        bed_prefix = "00.inputs/{gene}/{gene}_resampled_hg38_genio",
        region = lambda wc: gene_to_region[wc.gene]
    shell:
        """
        for g in 0 1 2; do
            exp_gz="01.sim_resampled/{wildcards.gene}/{wildcards.N_exposure}/{wildcards.N_outcome}/{wildcards.numCausalX}/{wildcards.numCausalY}/{wildcards.run}/gam0${{g}}_exposureSize_{wildcards.numCausalX}_outcomeSize_{wildcards.numCausalY}_simNum_{wildcards.run}.MRLink_exp.gz"
            out_gz="01.sim_resampled/{wildcards.gene}/{wildcards.N_exposure}/{wildcards.N_outcome}/{wildcards.numCausalX}/{wildcards.numCausalY}/{wildcards.run}/gam0${{g}}_exposureSize_{wildcards.numCausalX}_outcomeSize_{wildcards.numCausalY}_simNum_{wildcards.run}.MRLink_out.gz"

            exp_tmp="${{exp_gz%.gz}}_tmp.tsv"
            out_tmp="${{out_gz%.gz}}_tmp.tsv"

            out_res="01.sim_resampled/{wildcards.gene}/{wildcards.N_exposure}/{wildcards.N_outcome}/{wildcards.numCausalX}/{wildcards.numCausalY}/{wildcards.run}/gam0${{g}}_MRLink.res.txt"

            zcat $exp_gz > $exp_tmp
            zcat $out_gz > $out_tmp

            python3 {MRLINK2} \
                --reference_bed {params.bed_prefix} \
                --sumstats_exposure $exp_tmp \
                --sumstats_outcome $out_tmp \
                --out $out_res \
                --maf 0.00000001 \
                --prespecified_regions {params.region} \
                --max_correlation 0.99 \
                --maf_threshold 0.000000001 

            rm -f $exp_tmp $out_tmp
        done
    """
 
rule gcta_exposure: 
    input:
        bed_file = "00.inputs/{gene}/{gene}_resampled_hg38_genio.bed",
        cojo_exp_gam0 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam00_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.cojo_exp.tsv.gz",
        cojo_exp_gam1 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam01_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.cojo_exp.tsv.gz",
        cojo_exp_gam2 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam02_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.cojo_exp.tsv.gz"
    output:
        cojo_out_gam0 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/cojo/gam00_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.jma.cojo",
        cojo_out_gam1 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/cojo/gam01_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.jma.cojo",
        cojo_out_gam2 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/cojo/gam02_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.jma.cojo"
    params:
        bed_prefix = "00.inputs/{gene}/{gene}_resampled_hg38_genio"
    shell:
        """
        for g in 0 1 2; do
            
            exp_gz="01.sim_resampled/{wildcards.gene}/{wildcards.N_exposure}/{wildcards.N_outcome}/{wildcards.numCausalX}/{wildcards.numCausalY}/{wildcards.run}/gam0${{g}}_exposureSize_{wildcards.numCausalX}_outcomeSize_{wildcards.numCausalY}_simNum_{wildcards.run}.cojo_exp.tsv"
            out_path="01.sim_resampled/{wildcards.gene}/{wildcards.N_exposure}/{wildcards.N_outcome}/{wildcards.numCausalX}/{wildcards.numCausalY}/{wildcards.run}/cojo/gam0${{g}}_exposureSize_{wildcards.numCausalX}_outcomeSize_{wildcards.numCausalY}_simNum_{wildcards.run}"
            
            exp_tmp="${{exp_gz%.gz}}_tmp.tsv"
            
            zcat $exp_gz > $exp_tmp

            gcta \
                --bfile {params.bed_prefix} \
                --cojo-file $exp_tmp \
                --cojo-slct \
                --cojo-p 1e-6 \
                --maf 0.00001 \
                --out $out_path
                
            if [ ! -f ${{out_path}}.jma.cojo ]; then
                echo -e "Chr\tSNP\tbp\tfreq\trefA\tb\tse\tp\tn\tfreq_geno\tbJ\tbJ_se\tpJ\tLD_r" > ${{out_path}}.jma.cojo
            fi
            
            rm -f $exp_tmp
        done
    """
    
rule gcta_outcome:
    input:
        bed_file = "00.inputs/{gene}/{gene}_resampled_hg38_genio.bed",
        cojo_out_gam0 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam00_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.cojo_out.tsv.gz",
        cojo_out_gam1 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam01_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.cojo_out.tsv.gz",
        cojo_out_gam2 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam02_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.cojo_out.tsv.gz"
    output:
        cojo_out_gam0 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/cojo/gam00_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.out.jma.cojo",
        cojo_out_gam1 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/cojo/gam01_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.out.jma.cojo",
        cojo_out_gam2 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/cojo/gam02_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.out.jma.cojo"
    params:
        bed_prefix = "00.inputs/{gene}/{gene}_resampled_hg38_genio"
    shell:
        """
        for g in 0 1 2; do
            
            out_gz="01.sim_resampled/{wildcards.gene}/{wildcards.N_exposure}/{wildcards.N_outcome}/{wildcards.numCausalX}/{wildcards.numCausalY}/{wildcards.run}/gam0${{g}}_exposureSize_{wildcards.numCausalX}_outcomeSize_{wildcards.numCausalY}_simNum_{wildcards.run}.cojo_out.tsv"
            out_path="01.sim_resampled/{wildcards.gene}/{wildcards.N_exposure}/{wildcards.N_outcome}/{wildcards.numCausalX}/{wildcards.numCausalY}/{wildcards.run}/cojo/gam0${{g}}_exposureSize_{wildcards.numCausalX}_outcomeSize_{wildcards.numCausalY}_simNum_{wildcards.run}.out"
            
            out_tmp="${{out_gz%.gz}}_tmp.tsv"
            
            zcat $out_gz > $out_tmp

            gcta \
                --bfile {params.bed_prefix} \
                --cojo-file $out_tmp \
                --cojo-slct \
                --cojo-p 1e-6 \
                --maf 0.00001 \
                --out $out_path
            
            if [ ! -f ${{out_path}}.jma.cojo ]; then
                echo -e "Chr\tSNP\tbp\tfreq\trefA\tb\tse\tp\tn\tfreq_geno\tbJ\tbJ_se\tpJ\tLD_r" > ${{out_path}}.jma.cojo
            fi
            
            rm -f $out_tmp
        done
    """
    
rule cisMR_cML:
    input:
        ld_file = "00.inputs/{gene}/{gene}_LD_resampled.rds",
        af_file = "00.inputs/{gene}/{gene}_af_resampled.tsv",
        cojo_exp_gam0 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/cojo/gam00_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.jma.cojo",
        cojo_exp_gam1 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/cojo/gam01_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.jma.cojo",
        cojo_exp_gam2 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/cojo/gam02_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.jma.cojo",
        cojo_out_gam0 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/cojo/gam00_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.out.jma.cojo",
        cojo_out_gam1 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/cojo/gam01_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.out.jma.cojo",
        cojo_out_gam2 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/cojo/gam02_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.out.jma.cojo",
        cojo_out_gwas_gam0 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam00_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.cojo_out.tsv.gz",
        cojo_out_gwas_gam1 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam01_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.cojo_out.tsv.gz",
        cojo_out_gwas_gam2 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam02_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.cojo_out.tsv.gz",
        cojo_exp_gwas_gam0 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam00_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.cojo_exp.tsv.gz",
        cojo_exp_gwas_gam1 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam01_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.cojo_exp.tsv.gz",
        cojo_exp_gwas_gam2 = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/gam02_exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.cojo_exp.tsv.gz"
    output:
        "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/exposureSize_{numCausalX}_outcomeSize_{numCausalY}_simNum_{run}.cisMRcML.res.tsv.gz"
    params:
        gwasX_path = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/",
        gwasY_path = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/",
        cojo_exp_path = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/cojo/",
        cojo_out_path = "01.sim_resampled/{gene}/{N_exposure}/{N_outcome}/{numCausalX}/{numCausalY}/{run}/cojo/"
    shell:
        "Rscript scripts/01_02.cisMR_cML.R  --run {wildcards.run} \
        --gwasX_path {params.gwasX_path} --gwasY_path {params.gwasY_path} --cojo_exp_path {params.cojo_exp_path} \
        --cojo_out_path {params.cojo_out_path} --ld_file {input.ld_file} --numCausalX {wildcards.numCausalX} \
        --numCausalY {wildcards.numCausalY} --af_file {input.af_file} --N_exposure {wildcards.N_exposure} --output {output}"
