# Extracting the tumor bam and normal bam files
wget 'https://s3-eu-west-1.amazonaws.com/wtsi-pancancer/testdata/HCC1143_ds.tar'
tar xvf HCC1143_ds.tar

# Extracting reference files
wget 'https://dcc.icgc.org/api/v1/download?fn=/PCAWG/reference_data/pcawg-bwa-mem/genome.fa.gz'
bgzip -d genome.fa.gz

# Getting the .excl file
wget 'https://github.com/dellytools/delly/blob/master/excludeTemplates/human.hg19.excl.tsv'
mv human.hg19.excl.tsv hg19.excl
