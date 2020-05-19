# Implementation of Somatic Variant-Calling Pipelines in Parallel

### Steps to run the pipeline in Python:

<b>Step 1:</b> Install all dependencies

```zsh
git clone https://github.com/olopade-lab/sv-pipeline
bash ./dependencies.sh
```

Note: If you are using Lumpy, don't forget to export the path of the bin with the real as shown below:

```zsh
export PATH=$PATH:~/lumpy-sv/bin
```

<b>Step 2:</b> Build your own inputs.json file or run the get_files.sh script to download test files. If you make your own inputs.json file, replace the original one provided in this repository. Make sure to specify 2 inputs in the json file for every set of samples. Example: normal-bam-2, tumor-bam-2. The script will find the number of samples on its own and build the samples.tsv file automatically as well.

```zsh
bash ./get_files.sh
```

<b>Step 3:</b> Run the python file with inputs.json as secondary file, Parsl configuration as third file (Select the right configuration from the configs folder) and type of caller as the third parameter.

Type of callers current available in this script: delly, lumpy

```python
python3 run_caller.py test_inputs.json igsb_jupyter delly
```

<hr/>

### Steps to run the pipeline in Docker:

<b>Step 1:</b> CD into a folder that contains all the .bam files, the index files and the genome reference file. If you want to use the test inputs, you may run the following commmand to download them:

```zsh
bash ./get_files.sh
```

<b>Step 2:</b> Pull the Docker image.
```zsh
docker pull dellytools/delly
```

<b>Step 3:</b> Initialise and run the Docker Image. (Note, I have initialised the root directory in the present working directory of the user)
```zsh
docker run -it -v $(pwd):/root dellytools/delly  
```

<b>Step 4:</b> Run Delly inside the container:

```zsh
delly call -o /root/sv.bcf -g /root/genome.fa /root/HCC1143_ds/HCC1143.bam /root/HCC1143_ds/HCC1143_BL.bam
```

<b>Step 5:</b> Exit the container:

```zsh
exit
```
