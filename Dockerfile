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
	build-essential \
	software-properties-common \
	bash \
	wget \
	curl \
	git \
	bzip2 \
	libbz2-dev \
	libgsl-dev \
	zlib1g \
	zlib1g-dev \
	liblzma-dev \
	libssl-dev \
	libncurses5-dev \
	libz-dev \
	python3-dev \
	python3-pip \ 
	libjemalloc-dev \
	cmake \
	make \
	g++ \
	libhts-dev \
	libzstd-dev \
	autoconf \
	libatomic-ops-dev \
	pkg-config \
	libomp5 \
	libomp-dev \
	libssl-dev \
	pkg-config \
	zip \
	unzip

##install fastp latest
RUN mkdir fastp \
	&& cd fastp \
	&& wget --no-check-certificate http://opengene.org/fastp/fastp \
	&& chmod a+x ./fastp \
	&& cd ..

ENV PATH /opt/fastp:$PATH

#install vcfanno 0.3.5
RUN mkdir -p vcfanno \
	&& cd vcfanno \
	&& wget https://github.com/brentp/vcfanno/releases/download/v0.3.5/vcfanno_linux64 \
	&& mv vcfanno_linux64 vcfanno \
	&& chmod a+x ./vcfanno \
	&& cd ..

ENV PATH /opt/vcfanno:$PATH

##install minimap2 2.28
RUN wget https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 \
	&& tar -jxvf minimap2-2.28_x64-linux.tar.bz2 \
	&& rm minimap2-2.28_x64-linux.tar.bz2

ENV PATH /opt/minimap2-2.28_x64-linux:$PATH

#install samblaster 0.1.26
RUN git clone git clone git://github.com/GregoryFaust/samblaster.git \
	&& cd samblaster \
	&& git checkout b642639117eafedc760d8b84c0d2c4872b0da084 \
	&& make

ENV PATH /opt/samblaster:$PATH
