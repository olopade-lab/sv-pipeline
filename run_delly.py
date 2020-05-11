#!/usr/bin/env python3
'''
This script executes the DELLY Somatic Variant Calling pipeline on multiple given pairs of samples in parallel.

Features:
1. Identifies the input json and extracts the reference genome along with the excl. files for genes to be excluded.
2. Runs the Parsl app for the DELLY pipeline in parallel for all given samples.
3. Final .bcf files are placed in the same directory as the application.
'''
import json
import subprocess
import csv
import os
import importlib
import logging
import argparse
logging.basicConfig(level=logging.INFO)
import parsl

''' Processing the Inputs Json file and extracting the reference genome and excl file '''

parser = argparse.ArgumentParser(description='DELLY Implementation in Parallel')
parser.add_argument('input_json', metavar='input_file', type= str, help='The input json file containing file paths to bam files.')
parser.add_argument('config', metavar='--config', type= str, help='The Parsl configuration for execution.')
args = parser.parse_args()

base_dir = '/'.join(os.path.abspath(__file__).split('/')[:-1])
config = os.path.join(base_dir, 'configs', '{}.py'.format(args.config))
spec = importlib.util.spec_from_file_location('', config)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
parsl.load(module.config)

inputs = json.load(open(args.input_json))
logging.info('Initialised Parsl and inputs.json file')

ref_fa = inputs['reference-fa']['path'] # type: str
excl = inputs['exclude']['path'] # type: str

logging.info('Reference file and excl file loaded.')

@python_app
def run_delly(tumor_bam: str, normal_bam: str , sample_num: int) -> None:
    '''
    This function takes the normal_bam and tumor_bam file strings and adds them to subprocess calls below.
    '''
    sv_discovery = 'delly call -x {} -o t1.bcf -g {} {} {}'.format(excl, ref_fa,tumor_bam, normal_bam)
    subprocess.call(sv_discovery, shell=True)
    subprocess.call('dos2unix samples.tsv', shell=True)
    subprocess.call('delly filter -f somatic -o t1.pre.bcf -s samples.tsv t1.bcf', shell=True)
    pre_filter_call =  'delly call -g {} -v t1.pre.bcf -o geno.bcf -x {} {} {}'.format(ref_fa, excl,tumor_bam, normal_bam)
    subprocess.call(pre_filter_call, shell=True)
    subprocess.call('delly filter -f somatic -o t1.somatic.bcf -s sample{}.tsv geno.bcf'.format(sample_num), shell=True)

@python_app
def make_tsv(tumor_bam:str, normal_bam: str, sample_num: int) -> None:
    '''
    This function takes the normal_bam and tumor_bam file strings and makes the .tsv file for that sample.
    '''
    tumor = tumor_bam.split('.bam')[0]
    control - normal_bam.split('.bam')[0]
    import csv
    out_file = open('sample{}.tsv'.format(sample_num), 'wt')
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow([tumor, 'tumor'])
    tsv_writer.writerow([control, 'control'])
    out_file.close()

num_samples = int((len(inputs) - 2)/2) # type: int
samples = [] # type: list
results = [] # type: list

logging.info('{} Samples detected'.format(num_samples))

for i in range(1,num_samples+1):
    tumor_bam = inputs['tumor-bam-{}'.format(i)]['path']
    normal_bam = inputs['normal-bam-{}'.format(i)]['path']

    samples.append(make_tsv(tumor_bam, normal_bam, sample_num))
    results.append(run_delly(tumor_bam, normal_bam, sample_num))
    logging.info('Compiling DELLY for Sample Number {}'.format(i))

logging.info('Executing DELLY using Parsl')
samples = [i.result() for i in samples]
results = [i.result() for i in results]
