# Use the standard Microsoft DevContainer Ubuntu base image
FROM mcr.microsoft.com/devcontainers/base:ubuntu

LABEL maintainer="Edoardo Giacopuzzi <edoardo.giacopuzzi@fht.org>"

# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive
USER root

# 1. Install system dependencies required for Micromamba/Bioinformatics
RUN apt-get update && apt-get install -y --no-install-recommends \
	wget \
	bzip2 \
	ca-certificates \
	curl \
	git \
	&& apt-get clean && rm -rf /var/lib/apt/lists/*

# 2. Set up Mamba (Conda) environment variables
# We install to a location accessible by the 'vscode' user
ENV MAMBA_ROOT_PREFIX="/opt/conda"
ENV PATH="${MAMBA_ROOT_PREFIX}/bin:${PATH}"

# 3. Install Micromamba
# We verify the binary and install it to /usr/local/bin for global access
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba \
	&& mv bin/micromamba /usr/local/bin/micromamba

# 4. Create the directory
RUN mkdir -p ${MAMBA_ROOT_PREFIX}

# 6. Create the environment and install tools
# We create a 'bioinfo' env but also initialize the shell so commands work immediately
# - htslib includes: bgzip, tabix
# - bwa includes: bwa mem
# - gatk4 includes: gatk >= 4
RUN micromamba create -y -n base -c bioconda -c conda-forge \
	samtools \
	bcftools \
	htslib \
	vcfanno \
	"gatk4>=4" \
	bwa \
	samblaster \
	&& micromamba clean --all --yes

# 7. Init shell for bash and zsh (standard in Codespaces)
RUN micromamba shell init --shell bash --root-prefix=${MAMBA_ROOT_PREFIX} \
	&& micromamba shell init --shell zsh --root-prefix=${MAMBA_ROOT_PREFIX}