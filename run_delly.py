#!/usr/bin/env python3
'''
This script executes the DELLY Somatic Variant Calling pipeline on multiple given pairs of samples in parallel.

Features:
1. Identifies the input json and extracts the reference genome along with the excl. files for genes to be excluded.
2. Runs the Parsl app for the DELLY pipeline in parallel for all given samples.
3. Final .bcf files are placed in the same directory as the application.
'''

import sys
import json
import subprocess
import csv
import logging
logging.basicConfig(level=logging.INFO)

import parsl
from parsl.data_provider.files import File
from parsl.app.app import python_app
from parsl.providers import LocalProvider
from parsl.channels import LocalChannel
from parsl.config import Config
from parsl.executors import HighThroughputExecutor

config = Config(
    executors=[
        HighThroughputExecutor(
            label="htex_local",
            cores_per_worker=1,
            provider=LocalProvider(
                channel=LocalChannel(),
                init_blocks=1,
                max_blocks=1,
            ),
        )
    ],
)

''' Loading the Parsl Configuration for a High Throughput Executor '''
parsl.load(config)

''' Processing the Inputs Json file and extracting the reference genome and excl file '''

inputs = json.load(open(sys.argv[1]))
logging.info('Initialised Parsl and inputs.json file')

ref_fa = inputs['reference-fa']['path'] # type: str
excl = inputs['exclude']['path'] # type: str

logging.info('Reference file and excl file loaded.')

@python_app
def run_delly(tumor_bam: str, normal_bam: str , samples_tsv: str) -> None:
    '''
    This function takes the normal_bam and tumor_bam file strings and adds them to subprocess calls below.
    '''

    sv_discovery = 'delly call -x {} -o t1.bcf -g {} {} {}'.format(excl, ref_fa,tumor_bam, normal_bam)
    subprocess.call(sv_discovery, shell=True)
    subprocess.call('dos2unix samples.tsv', shell=True)
    subprocess.call('delly filter -f somatic -o t1.pre.bcf -s samples.tsv t1.bcf', shell=True)
    pre_filter_call =  'delly call -g {} -v t1.pre.bcf -o geno.bcf -x {} {} {}'.format(ref_fa, excl,tumor_bam, normal_bam)
    subprocess.call(pre_filter_call, shell=True)
    subprocess.call('delly filter -f somatic -o t1.somatic.bcf -s {} geno.bcf'.format(samples_tsv), shell=True)


num_samples = int((len(inputs) - 2)/3) # type: int
results = [] # type: list

logging.info('{} Samples detected'.format(num_samples))

for i in range(1,num_samples+1):
    tumor_bam = inputs['tumor-bam-{}'.format(i)]['path']
    normal_bam = inputs['normal-bam-{}'.format(i)]['path']
    samples_tsv = inputs['tsv-{}'.format(i)]['path']
    results.append(run_delly(tumor_bam, normal_bam, samples_tsv))
    logging.info('Compiling DELLY for Sample Number {}'.format(i))

logging.info('Executing DELLY using Parsl')
results = [i.result() for i in results]
