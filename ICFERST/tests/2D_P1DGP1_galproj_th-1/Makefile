SHELL = sh

input:  clean
	gmsh -2 2d_mesh.geo

clean:
	rm -f *.vtu *stat *.msh
	rm -f *.log-0 *.err-0
