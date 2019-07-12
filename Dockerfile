FROM continuumio/miniconda3:4.5.12

# Copy the entry point for the user
COPY ./docker-entrypoint.sh /

# Install other dependencies (gcc)
RUN apt-get update
RUN apt-get install -y --no-install-recommends apt-utils
RUN apt-get -y install gcc zlib1g-dev liblzma-dev libcurl4-gnutls-dev libssl-dev

# Install pip
RUN conda install pip --yes

# Download betaduck
RUN git clone -b 19.01.1-1-dev https://github.com/alexiswl/betaduck.git
WORKDIR ./betaduck

# Upgrade pip
RUN pip install --upgrade pip

# Install environment
RUN conda env create --name betaduck --file environment.yaml

# Source env
# Pull the environment name out of the environment.yml
ENV PATH /opt/conda/envs/betaduck/bin:$PATH

# Install betaduck
RUN python setup.py install

# Install wub
RUN pip install git+https://github.com/nanoporetech/wub.git

# Change to /data directory
WORKDIR /data

# Change user
RUN useradd -ms /bin/bash docker
USER docker


# Set the entrypoint to be 'betaduck'
ENTRYPOINT ["/docker-entrypoint.sh"]
