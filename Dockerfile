FROM continuumio/miniconda3:latest

# Install other dependencies (gcc)
RUN apt-get update
RUN apt-get -y install gcc

# Update conda
RUN conda update -n base conda --yes
RUN conda update --all --yes

# Download poreduck
RUN git clone -b 18.07.1-3 https://github.com/alexiswl/betaduck.git
WORKDIR ./betaduck

# Install matplotlib_venn through pip
RUN pip install matplotlib_venn

# Install required packages
RUN conda install --file requirements.txt --yes

# Re-update conda
RUN conda update --all --yes

# Upgrade pip
RUN pip install --upgrade pip

# Install poreduck using pip
RUN pip install -e . --ignore-installed

# Install deconcatenate fastqs
RUN pip install ont-fastq-deconcatenate

# Copy the entry point for the user
COPY ./docker-entrypoint.sh /

# Change to /data directory
WORKDIR /data

# Change user
RUN useradd -ms /bin/bash docker
USER docker

# Set the entrypoint to be 'betaduck'
ENTRYPOINT ["/docker-entrypoint.sh"]
