SHELL = sh

SIM = sod_shock3d.mpml
input: clean

	gmsh -3 src/3d_mesh.geo
	gmsh -3 src/3d_mesh_f.geo
clean:
	rm -f *.d.* fluidity.* ocvn
	rm -f flog.dat time_count tempfile CheckSparsityMatrix.dat filedim fvalues mirror_int_data.dat fort.*
	rm -f *.node *.bound *.ele *.face *.stat *.vtu
	rm -f mirror_new.dat matrixdump matrixdump.info
	rm -f *.ele *.bound *.face *.node *.edge
	rm -f src/*.msh
