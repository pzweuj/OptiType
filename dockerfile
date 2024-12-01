#################################################################
# Dockerfile
#
# Version:          1.1
# Software:         OptiType
# Software Version: 1.3.5.p
# Description:      Accurate NGS-based 4-digit HLA typing
# Website:          https://github.com/pzweuj/OptiType/
# Tags:             Genomics
# Provides:         OptiType 1.3.5.p
# Base Image:       alpine:3.20.3
# Build Cmd:        docker build --rm -t ghcr.io/pzweuj/optitype .
# Pull Cmd:         docker pull ghcr.io/pzweuj/optitype
# Run Cmd:          docker run -v /path/to/file/dir:/data ghcr.io/pzweuj/optitype
#################################################################

# Source Image
FROM alpine:3.20.3

################## BEGIN INSTALLATION ###########################
# install
RUN apk add --no-cache \
    gcc \
    g++ \
    make \
    cmake \
    zlib-dev \
    bzip2-dev \
    git \
    python3 \
    python3-dev \
    py3-pip \
    curl \
    musl-dev \
    libffi-dev \
    hdf5 \
    hdf5-dev \
    xz-dev \
    glpk \
    glpk-dev

# Install biosc from source code
RUN git clone https://github.com/Blosc/c-blosc2.git \
    && cd c-blosc2 \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && make install \
    && cd ../.. \
    && rm -rf c-blosc2

# HLA Typing
# OptiType dependencies
RUN python3 -m venv /opt/venv && /opt/venv/bin/pip install --no-cache-dir --upgrade pip && /opt/venv/bin/pip install --no-cache-dir \
    numpy \
    pyomo \
    pysam \
    matplotlib \
    cython \
    numexpr \
    tables \
    pandas 

# installing optitype from git repository (version Dec 09 2015) and writing config.ini
RUN git clone -b dev https://github.com/pzweuj/OptiType.git \
    && sed -i -e '1i#!/usr/bin/env python3' OptiType/OptiTypePipeline.py \
    && mv OptiType/ /usr/local/bin/ \
    && chmod 777 /usr/local/bin/OptiType/OptiTypePipeline.py

# installing razers3
RUN git clone https://github.com/seqan/seqan.git seqan-src \
    && cd seqan-src \
    && cmake -DCMAKE_BUILD_TYPE=Release \
    && make razers3 \
    && cp bin/razers3 /usr/local/bin \
    && cd .. \
    && rm -rf seqan-src

ENV PATH=/usr/local/bin/OptiType:$PATH
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Change workdir to /data/
WORKDIR /data/

# Define default command
CMD ["/bin/sh"]

##################### INSTALLATION END ##########################

# File Author / Maintainer
LABEL author="schubert" maintainer="pzweuj" version="v1.3.5.p1"
