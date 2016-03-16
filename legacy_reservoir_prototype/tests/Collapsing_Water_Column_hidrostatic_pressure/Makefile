SHELL = sh

input:  clean
	gmsh -2 2d_mesh.geo

clean:
	rm -f *.vtu *stat *.msh matrixdump*
	rm -f fluidity.log-0 fluidity.err-0
