# DockerFile for a Fludity development container

# Use a focal base image
FROM ubuntu:focal

# This DockerFile is looked after by
MAINTAINER Tim Greaves


# Installs shouldn't expect input
ENV DEBIAN_FRONTEND=noninteractive

# Package updates and installs
RUN apt-get update && \
     apt-get -y dist-upgrade && \
     apt-get -y install gnupg dirmngr && \
     echo "deb http://ppa.launchpad.net/fluidity-core/ppa/ubuntu focal main" > /etc/apt/sources.list.d/fluidity-core-ppa-focal.list && \
     gpg --keyserver keyserver.ubuntu.com --recv 0D45605A33BAC3BE && \
     gpg --export --armor 33BAC3BE | apt-key add - && \
     apt-get update && \
     echo "Europe/London" > /etc/timezone && \
     apt-get -y install fluidity-dev texlive-pstricks texlive texlive-latex-extra texlive-science python3-pip && \
     rm -rf /var/cache/apt && \
     rm -rf /var/lib/apt/lists

RUN python3 -m pip install --upgrade junit-xml

WORKDIR /usr/local
RUN curl -fsL https://gmsh.info/bin/Linux/gmsh-4.8.3-Linux64.tgz | tar --strip-components=1 -zxf -

ENV PETSC_DIR /usr/lib/petscdir/3.12
ENV LDFLAGS -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi
ENV CPPFLAGS -I/usr/include/hdf5/openmpi
ENV FCFLAGS -I/usr/include
ENV OMPI_MCA_btl_vader_single_copy_mechanism none
ENV OMPI_MCA_rmaps_base_oversubscribe 1

# Add an ICFERST user who will be the default user for this container
# Make sure the user has a userid matching the host system
# -- pass this as an argument at build time
ARG userid=1000
RUN adduser --disabled-password --gecos "" -u $userid icferst

# Add unzip
RUN apt-get -y install unzip

# Add the multifluids source
RUN mkdir -p /usr/local/src
ADD multifluids_icferst-master.zip /tmp
WORKDIR /usr/local/src
RUN unzip /tmp/multifluids_icferst-master.zip
WORKDIR /usr/local/src/multifluids_icferst-master
RUN ./configure --enable-2d-adaptivity --prefix=/usr/local
RUN make install

# Change to the icferst user and home directory
USER icferst
WORKDIR /home/icferst



