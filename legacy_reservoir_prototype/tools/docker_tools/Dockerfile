# DockerFile for a Fludity development container

# Use a Trusty base image
FROM ubuntu:trusty

# This DockerFile is looked after by
MAINTAINER Tim Greaves

# Add the Fluidity repository
RUN echo "deb http://ppa.launchpad.net/fluidity-core/ppa/ubuntu trusty main" > /etc/apt/sources.list.d/fluidity-core-ppa-trusty.list

# Import the Launchpad PPA public key
RUN gpg --keyserver keyserver.ubuntu.com --recv 0D45605A33BAC3BE
RUN gpg --export --armor 33BAC3BE | apt-key add -

# Update the system
RUN apt-get update
RUN apt-get -y dist-upgrade

# Install Fluidity development environment
RUN apt-get -y install fluidity-dev

# Add an icferst user who will be the default user for this container
RUN adduser --disabled-password --gecos "" icferst
ENV PETSC_DIR /usr/lib/petscdir/3.6.3

# Add unzip
RUN apt-get -y install unzip

# Add the multifluids source
RUN mkdir -p /usr/local/src
ADD multifluids_monthly-master.zip /tmp
WORKDIR /usr/local/src
RUN unzip /tmp/multifluids_monthly-master.zip
WORKDIR /usr/local/src/multifluids_monthly-master
RUN ./configure --enable-2d-adaptivity --prefix=/usr/local
RUN make install

# Change to the icferst user and home directory
USER icferst
WORKDIR /home/icferst



