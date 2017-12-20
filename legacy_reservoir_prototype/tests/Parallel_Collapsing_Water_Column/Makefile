SHELL = sh

SIM = Checkpointing

input:  clean

	gmsh -2 -bin src/2d_mesh.geo
	../../../bin/fldecomp -n 2 -m gmsh src/2d_mesh


clean:
	rm -f *.d.* fluidity.*
	rm -f *.node *.bound *.ele *.face *.stat *.vtu *.pvtu 
	rm -rf cwc_P1DGP?_[0-9]
	rm -f matrixdump*
	rm -f src/*.msh src/*.halo
