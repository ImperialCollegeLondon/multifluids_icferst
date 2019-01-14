# DockerFile for a ICFerst runtime container
FROM doc-reg-ac.pds.nl/fluidity2rt:1.0

# This DockerFile is looked after by
MAINTAINER Sergey Koshelev, e-mail: sergey.koshelev@pds.nl

ENV PATH /usr/lib64/openmpi/bin:/usr/local/sbin:/usr/local/bin:/sbin:/bin:/usr/sbin:/usr/bin:/root/bin:/opt/icl/bin
ENV LD_LIBRARY_PATH /usr/lib64/openmpi/lib64/petsc/3.6.3/linux-gnu-c-opt/lib64:/usr/lib64/openmpi/lib
ENV PYTHONPATH=/glb/data/icl/lib/python2.7/site-packages

COPY ./icl/ /opt/icl/

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]