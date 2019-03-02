FROM continuumio/miniconda3:latest

# Install other dependencies (gcc)
RUN apt-get update
RUN apt-get -y install gcc zlib1g-dev liblzma-dev libcurl-dev libssl-dev

# Update conda
RUN conda update -n base conda --yes
RUN conda update --all --yes
# Install pip
RUN conda install pip --yes

# Download betaduck
RUN git clone -b 19.01.1-1-dev https://github.com/alexiswl/betaduck.git
WORKDIR ./betaduck

# Upgrade pip
RUN pip install --upgrade pip

# Install environment
RUN conda create environment.yaml --yes

# Source env
RUN source activate python_3.7

# Install betaduck
RUN python setup.py install

# Copy the entry point for the user
COPY ./docker-entrypoint.sh /

# Change to /data directory
WORKDIR /data

# Change user
RUN useradd -ms /bin/bash docker
USER docker

# Set the entrypoint to be 'betaduck'
ENTRYPOINT ["/docker-entrypoint.sh"]
