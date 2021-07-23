FROM ubuntu

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
apt-get install -y \
mpich \
libhdf5-mpich-dev \
doxygen \
cmake \
git \
clang-10 \
clang-format-10 \
python3-pip \
python3.6 \
gcovr

RUN pip3 install \
h5py==2.10.0 \
mpmath==1.1.0 \
numpy==1.19.0 \
pandas==1.0.5 \
python-dateutil==2.8.1 \
pytz==2020.1 \
scipy==1.5.1 \
six==1.15.0 \
sympy==1.6.1