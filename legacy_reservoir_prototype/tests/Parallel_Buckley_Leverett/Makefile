SHELL = sh

SIM = Checkpointing

input:  clean

	gmsh -2 -bin src/unstructured_quartered.geo
	../../../bin/fldecomp -n 2 -m gmsh src/unstructured_quartered


clean:
	rm -f *.d.* fluidity.*
	rm -f *.node *.bound *.ele *.face *.stat *.vtu *.pvtu 
	rm -rf BL_radial_0
	rm -f matrixdump*
	rm -f src/*.msh src/*.halo
