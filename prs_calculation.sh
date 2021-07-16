
# Calculation of polygenic risk scores

plink2 --bfile ${bed_file}
--score ${weightfile} column_for_variant_id column_for_effect_allele column_for_effect_size header list-variants cols=scoresums
--out ${resultfile}

# for the genotype input, depending on dataset genotype format, bfile can be replaced with --pfile/--vcf/--bgen with respective commands
