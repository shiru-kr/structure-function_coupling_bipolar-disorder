#!/bin/sh

# This script runs the QC process using PLINK on the lifted over merged PLINK files generated
# using liftover_plinks_bed.wdl script

data_file_dir="/Data/lifted_arrays"

run_plink_qc="plink2 --bfile ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged \
 --keep coupling_df.phe --autosome \
 --max-alleles 2 --maf 0.01 --mac 100 --geno 0.1 \
 --hwe 1e-6 --mind 0.1 --write-snplist \
 --write-samples --no-id-header --out array_snps_qc_pass"
 
dx run swiss-army-knife \
  -iin=${data_file_dir}/ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged.bed \
  -iin=${data_file_dir}/ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged.bim \
  -iin=${data_file_dir}/ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged.fam \
  -iin=coupling_df.phe \
  -icmd="${run_plink_qc}" \
  --tag="Array QC" \
  --instance-type "mem1_ssd1_v2_x36" \
  --priority high \
  --destination="array_qc/" \
  --brief --yes

