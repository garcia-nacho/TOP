library(readxl)
system("wget 'https://ngstar.canada.ca/alleles/download?lang=en&loci_name=23S' -O /home/docker/pyngSTar/pyngSTarDB_rolling/23S_alleles.fasta")
system("wget 'https://ngstar.canada.ca/alleles/download?lang=en&loci_name=penA' -O /home/docker/pyngSTar/pyngSTarDB_rolling/penA_alleles.fasta")
system("wget 'https://ngstar.canada.ca/alleles/download?lang=en&loci_name=mtrR' -O /home/docker/pyngSTar/pyngSTarDB_rolling/mtrR_alleles.fasta")
system("wget 'https://ngstar.canada.ca/alleles/download?lang=en&loci_name=porB' -O /home/docker/pyngSTar/pyngSTarDB_rolling/porB_alleles.fasta")
system("wget 'https://ngstar.canada.ca/alleles/download?lang=en&loci_name=ponA' -O /home/docker/pyngSTar/pyngSTarDB_rolling/ponA_alleles.fasta")
system("wget 'https://ngstar.canada.ca/alleles/download?lang=en&loci_name=gyrA' -O /home/docker/pyngSTar/pyngSTarDB_rolling/gyrA_alleles.fasta")
system("wget 'https://ngstar.canada.ca/alleles/download?lang=en&loci_name=parC' -O /home/docker/pyngSTar/pyngSTarDB_rolling/parC_alleles.fasta")


system("wget 'https://ngstar.canada.ca/sequence_types/download?lang=en' -O /home/docker/pyngSTar/pyngSTarDB_rolling/profiles.xlsx")

df<-read_xlsx("/home/docker/pyngSTar/pyngSTarDB_rolling/profiles.xlsx")
colnames(df)[1]<-"ST"
write.table(df, "/home/docker/pyngSTar/pyngSTarDB_rolling/ngstar_profiles.tab", row.names = FALSE, quote = FALSE, sep = "\t")
