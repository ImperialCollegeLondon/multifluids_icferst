SHELL = sh

SIM = Checkpointing

mesh:
	gmsh -2 -format msh2 -bin src/unstructured_quartered.geo

decomp: mesh
	../../../bin/fldecomp -n 2 -m gmsh src/unstructured_quartered

input: decomp

clean:
	rm -f *.d.* fluidity.*
	rm -f *.node *.bound *.ele *.face *.stat *.vtu *.pvtu 
	rm -rf BL_radial_0
	rm -f matrixdump*
	rm -f src/*.halo
