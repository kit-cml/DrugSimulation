FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    build-essential \
    ca-certificates \
    cmake \
    git \
    libcurl4-openssl-dev \
    libjson-c-dev \
    libopenmpi-dev \
    make \
    openmpi-bin \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace/DrugSimulation

CMD ["/bin/bash"]
