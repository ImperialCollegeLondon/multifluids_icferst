# DockerFile for a Fludity development container

# Use a Xenial base image# FOR DIFFERENT UBUNTU DISTROS REMEMBER TO CHANGE ALL THE XENIAL INSTANCES!
FROM ubuntu:xenial

# This DockerFile is looked after by
MAINTAINER Tim Greaves

# Add the Fluidity repository
RUN echo "deb http://ppa.launchpad.net/fluidity-core/ppa/ubuntu xenial main" > /etc/apt/sources.list.d/fluidity-core-ppa-xenial.list

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
ADD multifluids_icferst-master.zip /tmp
WORKDIR /usr/local/src
RUN unzip /tmp/multifluids_icferst-master.zip
WORKDIR /usr/local/src/multifluids_icferst-master
RUN ./configure --enable-2d-adaptivity --prefix=/usr/local
RUN make install

# Change to the icferst user and home directory
USER icferst
WORKDIR /home/icferst



