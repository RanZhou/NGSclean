clean=/scratch/twh58439/2016AthSA/NGSclean
cd $clean
mkdir $clean/rRNA_ref
ml STAR/2.5.3a-foss-2016b
/usr/local/apps/eb/STAR/2.5.3a-foss-2016b/bin/STAR \
  --runThreadN 2  \
  --runMode genomeGenerate  \
  --genomeDir rRNA_ref  \
  --genomeChrBinNbits  10 \
  --genomeSAindexNbases 4 \
  --genomeFastaFiles $clean/rRNA_only_NR.fas 
