

# -+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

# LDpred workflow

# 2019

# -+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

##### STEP 1

#!/bin/bash
phenotypes=( "" "" "" "" )
samplesize=( "" "" "" "")
ldradius=( "" "" "" "")

for ((k=0;k<${#phenotypes[@]};++k));
do
  for i in {1..22};
  do
  python coord_genotypes.py \
    --gf=eur1000G_chr.acgt.chr"${i}" \
    --ssf="${phenotypes[k]}".txt \
    --N="${samplesize[k]}" \
    --out="${phenotypes[k]}"_coord_chr"${i}" \
    --ssf_format=STANDARD
  done
done

##### STEP 2

# LDpred STEP 2
#!/bin/bash -l
for ((k=0;k<${#phenotypes[@]};++k)); do
  for i in {1..22};
  do
    # Analysis performed in google cloud but cloud maintenance and startup script not shown due to privacy reasons
    ldpred --coord="${phenotypes[k]}"_coord_chr"${i}" --ld_radius="${ldradius[k]}" --PS=1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001 --local_ld_file_prefix=LD_"${phenotypes[k]}"_chr"${i}" --N="${samplesize[k]}" --out=ldpred_"${phenotypes[k]}"_chr"${i}"
  done
done


##### Merge files to one in R

library(data.table)
phenos <- c("", "", "","")

for(pheno in phenos) {
  print(paste0("Starting ", pheno))

  setwd(paste0("ldpred_", pheno, "/"))

  pb.chr <- txtProgressBar(1, 22, style = 3)

  for (chr in c(1:22)) {
    print(paste0("// Processing ", pheno, ", chromosome ", chr))

    ldpred.prefix <- paste0(pheno, "_chr", chr, "_LDpred")
    ldpred.inf <- fread(paste(paste0(ldpred.prefix, "-inf"), "txt", sep = "."))
    ldpred.full <- ldpred.inf

    for (ldpred.p.file in dir(pattern = paste0(ldpred.prefix, "_"))) {

      ldpred.p <- fread(ldpred.p.file)
      ldpred.full <- cbind(ldpred.full, ldpred.p$ldpred_beta)
      names(ldpred.full)[ncol(ldpred.full)] <-
        paste("ldpred", gsub(".txt", "", strsplit(ldpred.p.file, "_")[[1]][5]), "beta", sep = "_")
    }

    print(paste0("// Writing table"))
    write.table(ldpred.full, paste0("LDpred_", pheno, "/", ldpred.prefix, ".txt"),
              sep = "\t", quote = FALSE, na = "", row.names = FALSE)

    setTxtProgressBar(pb.chr, chr)
  }
}
