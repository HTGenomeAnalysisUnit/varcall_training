FROM debian:bullseye-slim AS binary
LABEL description="varcall_training"
LABEL base_image="debian:bullseye-slim"
LABEL software="varcall_training"
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


##install minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 \
	&& tar -jxvf minimap2-2.28_x64-linux.tar.bz2 \
	&& rm minimap2-2.28_x64-linux.tar.bz2

ENV PATH /opt/minimap2-2.28_x64-linux:$PATH
