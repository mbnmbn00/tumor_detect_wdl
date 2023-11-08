# Tumor detection pipeline in WDL
- Author: Byoungnam Min
- Last updated: 2023-10-27

It is an example pipeline for tumor detection from DNA sequencing.
The pipeline consists of:
- Dockerfile
- WDL (Workflow Description Language)
- Python scripts used for running programs
  - `run_fastp.py` (for read QC)
  - `run_bwa.py` (for read alignment)
  - `run_mutect2.py` (for finding tumor-related mutations)

# Procedure

This test was performed on Ubuntu 22.04

## 1. Install Docker and Cromwell

```bash
# Docker
sudo apt-get update
sudo apt-get remove docker docker-engine docker.io
sudo apt-get install -y docker.io
sudo systemctl start docker
sudo systemctl enable docker

# Cromwell
wget https://github.com/broadinstitute/cromwell/releases/download/86/cromwell-86.jar
```

## 2. Download Github repository and build the Docker image

```bash
git clone https://github.com/mbnmbn00/tumor_detect_wdl.git
cd tumor_detect_wdl/scripts
docker build --platform linux/amd64 --tag tumor_img .
```

## 3. Download a toy dataset

We will use mice tumor samples as reported in [McCreery et al, 2015](https://www.nature.com/articles/nm.3979). For the purpose of speed, we will focus on chromosome 7 where known mutations are reported.

```bash
# Mice genome GRCm39 chromosome 7
# https://www.ncbi.nlm.nih.gov/nuccore/CM001000.3/
export NCBI_SEQ_ID=CM001000.3
python3 tumor_detect_wdl/scripts/download_ncbi_seq.py \
  --ncbi_seq_id ${NCBI_SEQ_ID} \
  --output_file chromosome_7.fasta
```

We will use only one sample for testing, but you can download the full list of samples from [ENA](https://www.ebi.ac.uk/ena/browser/view/PRJEB4767).

```bash
curl --remote-name ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR100/005/ERR1008135/ERR1008135_1.fastq.gz
curl --remote-name ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR100/005/ERR1008135/ERR1008135_2.fastq.gz
```

## 4. Run the pipeline

```bash
java -jar cromwell.jar run mutect2_pipeline.wdl --i my_inputs.json
```
