all: input
	fluidity  -v5 balzano1_2plus1.flml &> output.txt

input: clean

clean:
	rm -f  *.ele *.edge *.node *.poly *vtu *.s *.d.1 *.stat *.msh *.detectors matrixdump* *.pdf *.png *.dat
	rm -rf balzano1_2plus1_? balzano1_2plus1_?? balzano1_2plus1_???
	rm -rf *flredecomp*
	cd mesh; rm -f  *.ele *.edge *.node *.poly *.vtu *.s *.d.1 *.stat *.detectors matrixdump* *.halo 
config:
	diamond -s ../../tools/fluidity_options.rng  balzano1_2plus1.flml 

