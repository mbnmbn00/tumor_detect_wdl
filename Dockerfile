# Start from a base image with Ubuntu
FROM ubuntu:22.04

# Set noninteractive installation
ENV DEBIAN_FRONTEND=noninteractive

# Set the working directory
WORKDIR /app

# Install dependencies
RUN apt-get update && \
    apt-get install \
    -y samtools python3-pip wget unzip openjdk-18-jre

RUN pip3 install biopython pandas

# Download and install bwa-mem2
RUN wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 && \
    tar -xvjf bwa-mem2-2.2.1_x64-linux.tar.bz2
ENV PATH=$PATH:/app/bwa-mem2-2.2.1_x64-linux


# Download and install fastp
RUN wget http://opengene.org/fastp/fastp && \
    mv fastp /usr/local/bin/ && \
    chmod a+x /usr/local/bin/fastp

# Download GATK
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip && \
    unzip gatk-4.4.0.0.zip
ENV PATH=$PATH:/app/gatk-4.4.0.0
RUN ln -s /usr/bin/python3 /usr/bin/python

# Download Picard tools
RUN wget https://github.com/broadinstitute/picard/releases/download/3.1.0/picard.jar
ENV PATH=$PATH:/app

# Cleanup unnecessary files
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Add my script to the Docker image
COPY . /app
RUN chmod +x run_fastp.py
