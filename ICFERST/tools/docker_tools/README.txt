Docker can be downloaded from https://www.docker.com/


The scripts Run_fldecomp.py and Run_icferst.py are set up so it makes easier for the final user to run simulations using docker.

When running, by default the user is icferst, make sure it has rights to read/write that folder or change the user.

DOCKER USEFUL COMMANDS:
Loading the container is required to have it in the system and to actually be able to run icferst from docker.
To load the container, run:

    docker load --input icferst.tar


The dockerfile is set up so new builds can be created so long the file multifluids_monthly-master.zip is present in the same folder. 

To rebuild the container, in the directory with the Dockerfile and multifluids tarball run:

    docker build -t icferst .

To see local docker builds

    docker images

To create another icferst.tar in the current location

    docker save icferst > icferst.tar

To remove local images

    docker rmi IMAGENAME

To make a new distribution archive:

    docker save icferst > icferst.tar

