@ECHO OFF
rem OutputFile: %2\%1.log
del %2\%1.log
cd ..
del 1.tmp
echo "Generating output file ......" >>%2\%1.log
%3\GID_B3D.exe %2\%1.DAT %3\B3D.PAR
if exist 1.tmp goto end
echo "Running job......" >>%2\%1.log
"C:\Program Files (x86)\IC-QMUL\VGW\Y3D.exe"  %1.Y3D
echo "Converting into VTU files......" >>%2\%1.log
"C:\Program Files (x86)\IC-QMUL\VGW\m2vtu\m2vtu3D.exe" %1 %1.Y3D
echo "Visulizing the results......" >>%2\%1.log
"C:\Program Files (x86)\MayaVi\mayavi" -d %10.vtu -m SurfaceMap -f ExtractTensorComponents

:end 
del 1.tmp
del %2\%1.log

