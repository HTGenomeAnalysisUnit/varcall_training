FROM debian:bullseye-slim AS binary
LABEL description="internal variant calling course at HT"
LABEL base_image="debian:bullseye-slim"
LABEL software="varcall"
LABEL about.home="https://github.com/HTGenomeAnalysisUnit/varcall_training"
LABEL about.license="GPLv3"

ARG DEBIAN_FRONTEND=noninteractive
#install basic libraries and python

WORKDIR /opt

RUN apt-get update

RUN apt-get -y install \
	bash \
	curl \
	python3-dev \
	python3-pip

RUN apt-get -y clean all \
	&& apt-get -y purge \
	&& rm -rf /var/cache \
	&& rm -rf /var/lib/apt/lists/*

RUN curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda
RUN conda create -y -n varcallenv -c conda-forge -c bioconda \
	gatk4=4.6.1.0 \
	minimap2=2.28 \
	fastp=0.24.0 \
	samblaster=0.1.26 \
	bcftools=1.21 \
	samtools=1.21 \
	igv=2.17.4 \
	vcfanno=0.3.5 \
	pandas=1.2.3 \
	numpy=1.22.3 \
	pyranges=0.0.95 \
	htslib \
	zlib \
	cyvcf2=0.30.4 \
	ncls=0.0.57 \
	hap.py=0.3.15 \
	truvari=4.3.1
RUN echo "source activate varcallenv" > ~/.bashrc
ENV PATH /miniconda/envs/varcallenv/bin:$PATH
RUN pip install git+https://github.com/fakedrtom/SVAFotate.git
