import sys
import json
import subprocess
import csv

try:
    inputs = json.load(open(sys.argv[1]))
except:
    raise TypeError('Input JSON file not specified or incorrectly formatted.')

def initialise_delly():
    subprocess.call('make PARALLEL=1 -B src/delly', shell=True)
    subprocess.call('export OMP_NUM_THREADS={}'.format(number_of_samples), shell=True)

tumor_bam = inputs['tumor-bam']['path']
normal_bam = inputs['normal-bam']['path']
ref_fa = inputs['reference-fa']['path']
excl = inputs['exclude']['path']

number_of_samples = 1 # needs to be dynamic from inputs.json file

sv_discovery = 'delly call -x {} -o t1.bcf -g {} {} {}'.format(excl, ref_fa,tumor_bam, normal_bam)
subprocess.call(sv_discovery, shell=True)
subprocess.call('dos2unix samples.tsv', shell=True)
subprocess.call('delly filter -f somatic -o t1.pre.bcf -s samples.tsv t1.bcf', shell=True)
pre_filter_call =  'delly call -g {} -v t1.pre.bcf -o geno.bcf -x {} {} {}'.format(ref_fa, excl,tumor_bam, normal_bam)
subprocess.call(pre_filter_call, shell=True)
subprocess.call('delly filter -f somatic -o t1.somatic.bcf -s samples.tsv geno.bcf', shell=True)
