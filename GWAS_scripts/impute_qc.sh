#!/bin/bash

# Directory where your BGEN and SAMPLE files are located
#make sure you create all folders here and that all paths are correct

data_file_dir="/Bulk/Imputation/Imputation_from_genotype_GEL" #I changed the name of this forlder so make sure it fits yours

# Prefix for output files
output_prefix="imputed_snps_data_qc_pass"

# Destination folder for outputs on DNAnexus
output_destination="GWAS_output/impute_qc/"

# Plink2 QC options 
plink2_qc_options="--rm-dup exclude-all list --max-alleles 2 --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-6 --mind 0.1"

# Array for storing snplist outputs for downstream merging

echo "--- Starting Plink2 QC Job Submissions ---"
echo "Submitting Plink2 QC jobs for chromosomes 1 to 22..."

for chr in {1..22}; do

  run_plink2_qc=$(cat <<EOF
  plink2 \
  --import-max-alleles 2 \
  --bgen ukb21008_c${chr}_b0_v1.bgen ref-first \
  --sample ukb21008_c${chr}_b0_v1.sample \
  --keep coupling_df.phe \
  ${plink2_qc_options} \
  --write-snplist --write-samples --no-id-header \
  --out ${output_prefix}_chr${chr}

EOF
)


  echo "Submitting QC job for chromosome ${chr}..."
  dx run swiss-army-knife \
    -iin=${data_file_dir}/ukb21008_c${chr}_b0_v1.bgen \
    -iin=${data_file_dir}/ukb21008_c${chr}_b0_v1.sample \
    -iin=coupling_df.phe \
    -icmd="$run_plink2_qc" \
    --tag="plink2_qc_chr${chr}" \
    --destination="GWAS_output/impute_qc/" \
    --instance-type "mem1_ssd1_v2_x36" \
    --priority high \
    --brief --yes
done

echo "All per-chromosome Plink2 QC jobs have been submitted."
echo "--- Script Finished ---"

