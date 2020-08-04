FROM ubuntu:18.04

WORKDIR /opt

#### get necessary packages
RUN apt-get -y update && apt-get install -y wget build-essential zlib1g zlib1g-dev


#### download and compile cd-hit 4.8.1
RUN wget https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz && \
    tar xvf cd-hit-v4.8.1-2019-0228.tar.gz && \
    mv cd-hit-v4.8.1-2019-0228 cd-hit && \
    cd /opt/cd-hit && \
    make && \
    cd /opt/cd-hit/cd-hit-auxtools && \
    make 


#### get NCBI BLAST+ 2.8.1
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.8.1/ncbi-blast-2.8.1+-x64-linux.tar.gz && \
    tar xvf ncbi-blast-2.8.1+-x64-linux.tar.gz && \
    rm -f ncbi-blast-2.8.1+-x64-linux.tar.gz


#### set system PATH
ENV PATH="/opt/cd-hit:/opt/cd-hit/cd-hit-auxtools:/opt/cd-hit/psi-cd-hit:/opt/ncbi-blast-2.8.1+/bin:${PATH}"


