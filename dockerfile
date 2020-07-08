FROM ubuntu

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
apt-get install -y \
mpich \
libhdf5-mpich-dev \
doxygen \
cmake \
clang-format \
python3-pip \
python3.6
RUN pip3 install numpy
RUN pip3 install h5py
