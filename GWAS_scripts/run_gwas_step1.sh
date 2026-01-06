#!/bin/sh

# Define paths and filenames
data_file_dir="/Data/lifted_arrays"

# Construct REGENIE step 1 command
run_regenie_step1="regenie \
  --step 1 \
  --bed ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged \
  --extract array_snps_qc_pass.snplist \
  --phenoFile coupling_df.phe \
  --covarFile coupling_df.phe \
  --phenoColList temporal_pole,superior_frontal_gyrus,precentral_gyrus,supramarginal_gyrus,frontal_pole \
  --covarColList age,age_squared,pca_1,pca_2,pca_3,pca_4,pca_5,pca_6,pca_7,pca_8,pca_9,pca_10,pca_11,pca_12,pca_13,pca_14,pca_15,pca_16,pca_17,pca_18,pca_19,pca_20,rs_motion,d_motion,icv \
  --catCovarList sex,imaging_cite \
  --apply-rint \
  --bsize 1000 \
  --lowmem \
  --threads 8 \
  --out step_1 \
  --apply-rint" \
  

# Submit job via Swiss Army Knife
dx run swiss-army-knife \
  -iin=${data_file_dir}/ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged.bed \
  -iin=${data_file_dir}/ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged.bim \
  -iin=${data_file_dir}/ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged.fam \
  -iin=coupling_df.phe \
  -iin=/Data/array_snps_qc_pass.snplist \
  -icmd="${run_regenie_step1}" \
  --tag="ReGENIE Step 1" \
  --instance-type "mem1_ssd1_v2_x36" \
  --priority high \
  --destination="GWAS_output_full_rint/" \
  --brief --yes

