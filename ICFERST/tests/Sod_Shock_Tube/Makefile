SHELL = sh

SIM = Sod_Shock_Tube_2D

input: clean

	gmsh -2 src/2d_mesh_20layers.geo
	gmsh -2 src/2d_mesh_80layers.geo

clean:
	rm -f *.d.* fluidity.* ocvn
	rm -f flog.dat time_count tempfile CheckSparsityMatrix.dat filedim fvalues mirror_int_data.dat fort.*
	rm -f *.node *.bound *.ele *.face *.stat *.vtu
	rm -f mirror_new.dat matrixdump matrixdump.info
	rm -f *.ele *.bound *.face *.node *.edge
	rm -f src/*.msh
#	rm -f *_80layers.mpml
