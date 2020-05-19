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
from parsl import python_app, bash_app

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

try:
    inputs = json.load(open(args.input_json))
except:
    raise FileNotFoundError('The Input Json file was either not found or it is corrupted.')

logging.info('Initialised Parsl and inputs.json file')

ref_fa = inputs['reference-fa']['path'] # type: str
excl = inputs['exclude']['path'] # type: str

if not os.path.exists(ref_fa):
    raise FileNotFoundError('The Reference file was not found. Please denote the correct directory for the reference file.')

if not os.path.exists(excl):
    raise FileNotFoundError('The exclusion (excl) file was not found. Please denote the correct directory for the excl file.')

logging.info('Reference file and excl file loaded.')

@bash_app
def run_delly(tumor_bam: str, normal_bam: str , sample_num: int):
    '''
    This function takes the normal_bam and tumor_bam file strings and adds them to subprocess calls below.
    '''
    sv_discovery = 'delly call -x {} -o t{}.bcf -g {} {} {}'.format(excl, sample_num, ref_fa,tumor_bam, normal_bam)
    conversion = 'dos2unix samples.tsv'
    filter = 'delly filter -f somatic -o t{}.pre.bcf -s samples.tsv t{}.bcf'.format(sample_num, sample_num)
    pre_filter_call =  'delly call -g {} -v t{}.pre.bcf -o geno{}.bcf -x {} {} {}'.format(ref_fa, sample_num,
                                                                                          sample_num, excl,
                                                                                          tumor_bam, normal_bam)
    final_filter = 'delly filter -f somatic -o t{}.somatic.bcf -s sample{}.tsv geno{}.bcf'.format(sample_num, sample_num, sample_num)
    return '{} ; {} ; {} ; {} ; {}'.format(sv_discovery, conversion, filter, pre_filter_call, final_filter)

@python_app
def make_tsv(tumor_bam:str, normal_bam: str, sample_num: int) -> None:
    '''
    This function takes the normal_bam and tumor_bam file strings and makes the .tsv file for that sample.
    '''
    import csv

    tumor = tumor_bam.split('.bam')[0]
    control - normal_bam.split('.bam')[0]
    out_file = open('sample{}.tsv'.format(sample_num), 'wt')
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow([tumor, 'tumor'])
    tsv_writer.writerow([control, 'control'])
    out_file.close()

num_samples = int((len(inputs) - 2)/2) # type: int
samples = [] # type: list
results = [] # type: list

logging.info('{} Samples detected'.format(num_samples))

if num_samples < 1:
    raise ValueError('The Json file doesnt have any data for a given sample. It only has the reference file and excl file. \
                      Please add the normal bam and tumor bam file for each sample in the inputs json file.')

for i in range(1,num_samples+1):
    tumor_bam = inputs['tumor-bam-{}'.format(i)]['path']
    normal_bam = inputs['normal-bam-{}'.format(i)]['path']

    if not os.path.exists(tumor_bam) or not os.path.exists(normal_bam):
        raise FileNotFoundError('The Tumor bam or the Normal bam file couldnt be found. Please specify the right directory.')

    samples.append(make_tsv(tumor_bam, normal_bam, sample_num))
    results.append(run_delly(tumor_bam, normal_bam, sample_num))
    logging.info('Compiling DELLY for Sample Number {}'.format(i))

logging.info('Executing DELLY using Parsl')
parsl.wait_for_current_tasks()
