from parsl import python_app, bash_app

@bash_app
def run_delly(tumor_bam: str, normal_bam: str , sample_num: int):
    '''
    This function takes the normal_bam and tumor_bam file strings and runs the Delly SV Caller on them.
    '''
    sv_discovery = 'delly call -x {} -o t{}.bcf -g {} {} {}'.format(excl, sample_num, ref_fa,tumor_bam, normal_bam)
    conversion = 'dos2unix samples.tsv'
    filter = 'delly filter -f somatic -o t{}.pre.bcf -s samples.tsv t{}.bcf'.format(sample_num, sample_num)
    pre_filter_call =  'delly call -g {} -v t{}.pre.bcf -o geno{}.bcf -x {} {} {}'.format(ref_fa, sample_num,
                                                                                          sample_num, excl,
                                                                                          tumor_bam, normal_bam)
    final_filter = 'delly filter -f somatic -o t{}.somatic.bcf -s sample{}.tsv geno{}.bcf'.format(sample_num, sample_num, sample_num)
    conversion = 'bcftools view t{}.somatic.bcf > t{}.vcf'.format(sample_num, sample_num)
    return '{} ; {} ; {} ; {} ; {}; {}'.format(sv_discovery, conversion, filter, pre_filter_call, final_filter, conversion)

@bash_app
def run_lumpy(tumor_bam: str, normal_bam: str , sample_num: int):
    '''
    This function takes the normal_bam and tumor_bam file strings and runs the Lumpy SV-Caller on them.
    '''
    view_tumor = 'samtools view -b -F 1294 {} > tumor{}.discordants.unsorted.bam'.format(tumor_bam, sample_num)
    view_normal = 'samtools view -b -F 1294 {} > normal{}.discordants.unsorted.bam'.format(normal_bam, sample_num)

    extract_tumor = 'samtools view -h {} \
               | scripts/extractSplitReads_BwaMem -i stdin \
               | samtools view -Sb - \
               > tumor{}.splitters.unsorted.bam'.format(tumor_bam, sample_num)

    extract_normal = 'samtools view -h {} \
               | scripts/extractSplitReads_BwaMem -i stdin \
               | samtools view -Sb - \
               > normal{}.splitters.unsorted.bam'.format(normal_bam, sample_num)

    discord_tumor = 'samtools sort tumor{}.discordants.unsorted.bam tumor{}.discordants'.format(sample_num, sample_num)
    discord_normal = 'samtools sort normal{}.discordants.unsorted.bam normal{}.discordants'.format(sample_num, sample_num)

    splitter_tumor = 'samtools sort tumor{}.splitters.unsorted.bam tumor{}.splitters'.format(sample_num, sample_num)
    splitter_normal = 'samtools sort normal{}.splitters.unsorted.bam normal{}.splitters'.format(sample_num, sample_num)


    final_call = 'lumpyexpress \
                    -B tumor{}.bam, normal{}.bam \
                    -S tumor{}.splitters.bam, normal{}.splitters.bam \
                    -D tumor{}.discordants.bam, normal{}.discordants.bam \
                    -o tumor_normal{}.vcf'.format(sample_num, sample_num, sample_num, sample_num,
                                                  sample_num, sample_num, sample_num)

    return '{}; {}; {}; {}; {}; {}; {}; {}; {}'.format(view_tumor, view_normal, extract_tumor, extract_normal,
                                                    discord_tumor, discord_normal, splitter_tumor, splitter_normal, final_call)


@python_app
def make_tsv(tumor_bam:str, normal_bam: str, sample_num: int) -> None:
    '''
    This function takes the normal_bam and tumor_bam file strings and makes the .tsv file for that sample.
    '''
    import csv

    tumor = tumor_bam.split('.bam')[0]
    control = normal_bam.split('.bam')[0]
    out_file = open('sample{}.tsv'.format(sample_num), 'wt')
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow([tumor, 'tumor'])
    tsv_writer.writerow([control, 'control'])
    out_file.close()

@bash_app
def merge_vcf_files(sample_num: int) -> None:
    import glob
    import os
    pwd = os.getcwd()

    files = glob.glob('*{}.vcf'.format(sample_num))
    os.mkdir("/{}vcf".format(sample_num))

    for file in files:
        if not os.path.exists(file):
            raise FileNotFoundError('The vcf file was expected but not found: {}'.format(file))

        shutil.move(file, os.path.join(pwd, "/{}vcf".format(sample_num), file))

    return './SURVIVOR merge {}vcf 1000 2 1 1 0 30 sample{}_merged.vcf'.format(sample_num, sample_num)
