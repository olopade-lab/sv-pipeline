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
from caller_apps import run_delly, run_lumpy, make_tsv, merge_vcf_files

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

    results.append(run_lumpy(tumor_bam, normal_bam, sample_num))
    logging.info('Compiling DELLY for Sample Number {}'.format(i))

    else:
        raise ValueError('Please specify one of the available callers from the Readme document.')

logging.info('Executing DELLY using Parsl')
parsl.wait_for_current_tasks()

merging_results = []

logging.info('Proceeding to merging vcf files for each sample.')

for i in range(1,num_samples+1):
    merging_results.append(merge_vcf_files(i))

parsl.wait_for_current_tasks()
