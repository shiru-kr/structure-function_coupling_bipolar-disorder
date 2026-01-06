#!/bin/bash

# This script runs the second step of a REGENIE analysis for GWAS using imputed genotype data from GEL.
# It processes each chromosome from 1 to 22, extracting relevant SNPs and phenotypes, and outputs the results.
# .sample files had to be created with a job pn RAP because UKB ones had a fault (known issue)

data_file_dir="/Bulk/Imputation/Imputation_from_genotype_GEL"

loco_files=(
  "step_1_1.loco"
  "step_1_2.loco"
  "step_1_3.loco"
  "step_1_4.loco"
  "step_1_5.loco"
)

for chr in {1..22}; do

  run_regenie_step2=$(cat <<EOF 
  regenie --step 2 \
  --bgen ukb21008_c${chr}_b0_v1.bgen \
  --sample ukb21008_c${chr}_b0_v1.sample \
  --extract gel_imputed_snps_data_qc_pass.snplist \
  --phenoFile coupling_df.phe \
  --covarFile coupling_df.phe \
  --phenoColList temporal_pole,superior_frontal_gyrus,precentral_gyrus,supramarginal_gyrus,frontal_pole \
  --covarColList age,age_squared,pca_1,pca_2,pca_3,pca_4,pca_5,pca_6,pca_7,pca_8,pca_9,pca_10,pca_11,pca_12,pca_13,pca_14,pca_15,pca_16,pca_17,pca_18,pca_19,pca_20,rs_motion,d_motion,icv \
  --catCovarList sex,imaging_cite \
  --pred step_1_pred.list \
  --apply-rint \
  --bsize 500 \
  --apply-rint \
  --out step_2_chr${chr}
EOF
)
 
  dx run swiss-army-knife \
    -iin=${data_file_dir}/ukb21008_c${chr}_b0_v1.bgen \
    -iin=${data_file_dir}/ukb21008_c${chr}_b0_v1.sample \
    -iin=/Data/gel_imputed_snps_data_qc_pass.snplist \
    -iin=coupling_df.phe \
    -iin=/Data/GWAS_output_full_rint/step_1_pred.list \
    -iin=/Data/GWAS_output_full_rint/step_1_1.loco \
    -iin=/Data/GWAS_output_full_rint/step_1_2.loco \
    -iin=/Data/GWAS_output_full_rint/step_1_3.loco \
    -iin=/Data/GWAS_output_full_rint/step_1_4.loco \
    -iin=/Data/GWAS_output_full_rint/step_1_5.loco \
    -icmd="$run_regenie_step2" \
    --tag="regenie_step2_chr${chr}" \
    --destination="GWAS_output_full_rint/" \
    --instance-type "mem1_ssd1_v2_x16" \
    --brief --yes
done

