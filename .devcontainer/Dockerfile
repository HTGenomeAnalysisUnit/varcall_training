# Use the standard Microsoft DevContainer Ubuntu base image
FROM mcr.microsoft.com/devcontainers/base:ubuntu-22.04

LABEL maintainer="Edoardo Giacopuzzi <edoardo.giacopuzzi@fht.org>"

# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive

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
	&& mv bin/micromamba /usr/local/bin/micromamba \
	&& rm -rf bin

# 4. Create the directory and fix permissions for the 'vscode' user
# This allows the non-root user to install packages if needed later
RUN mkdir -p ${MAMBA_ROOT_PREFIX} \
	&& chown -R vscode:vscode ${MAMBA_ROOT_PREFIX} \
	&& chmod -R 775 ${MAMBA_ROOT_PREFIX}

# 5. Switch to 'vscode' user to install the bioinformatics tools
# This ensures the environment belongs to the user, not root
USER vscode

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

# Switch back to root to finish any standard devcontainer setup if needed
# (Codespaces usually handles the rest)
USER root