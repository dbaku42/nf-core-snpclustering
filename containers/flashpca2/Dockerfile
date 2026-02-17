# containers/flashpca2/Dockerfile
FROM ubuntu:22.04 AS build

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential git wget ca-certificates \
    libboost-program-options-dev libeigen3-dev \
 && rm -rf /var/lib/apt/lists/*

# Spectra header-only
RUN wget -qO- https://github.com/yixuan/spectra/archive/v0.8.1.tar.gz | tar xz

# flashpca
RUN git clone --depth 1 https://github.com/gabraham/flashpca.git /src/flashpca
WORKDIR /src/flashpca
RUN make all \
    EIGEN_INC=/usr/include/eigen3 \
    BOOST_INC=/usr/include/boost \
    SPECTRA_INC=/spectra-0.8.1/include

FROM ubuntu:22.04
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 \
    libboost-program-options1.74.0 \
    python3 \
 && rm -rf /var/lib/apt/lists/*

COPY --from=build /src/flashpca/flashpca /usr/local/bin/flashpca
CMD ["bash"]
