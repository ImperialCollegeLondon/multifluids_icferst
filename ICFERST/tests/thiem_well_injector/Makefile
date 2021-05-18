SHELL = sh

SIM = thiem_well_model

input: clean

	gmsh -format msh2 -3 cylinder.geo

clean:
	rm -f *.stat *.vtu *.msh
