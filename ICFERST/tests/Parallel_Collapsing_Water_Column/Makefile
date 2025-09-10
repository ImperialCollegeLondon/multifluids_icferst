SHELL = sh

SIM = Checkpointing

mesh:
	gmsh -2 -format msh2 -bin src/2d_mesh.geo

decomp: mesh
	../../../bin/fldecomp -n 2 -m gmsh src/2d_mesh

input: decomp

clean:
	rm -f *.d.* fluidity.*
	rm -f *.node *.bound *.ele *.face *.stat *.vtu *.pvtu 
	rm -rf cwc_P1DGP?_[0-9]
	rm -f matrixdump*
	rm -f src/*.halo
