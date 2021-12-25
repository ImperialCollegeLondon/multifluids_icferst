#!/usr/bin/env python3

# arguments:: input file for icferst
# This script runs the input test case (.mpml file) using the docker package icferst.
# Example: python Run_icferst.py test.mpml
import sys
import os

try:
	testname = sys.argv[1]
except:
	print('An input file must be introduced')
	print('For an example run with -h')
	exit()

if (testname[-2:] == '-h'):
	print('This script runs the input test case (.mpml file) using the docker package icferst.')
	print('The terminal HAS to be where the input file is located.')
	print('Example: python Run_icferst.py test.mpml')
	print('By default docker runs as the user icferst. Make sure it has the rights to write and read in the working folder.')
	exit()

print('Running the test case: ' + testname)
pwdpath = os.getcwd()
run_command  = "docker run -v " + os.getcwd()+":/rundir" # run and tell docker which folder it can see
run_command += " -w /rundir" #sets up the working directory
run_command += " -a stdout -t icferst" #Connects the output of the program to this terminal and specify the docker enviroment
run_command += " icferst" #name of the executable
run_command += " /rundir/"+testname
run_command = run_command.replace("./","")#Remove any potential "./" to access a file
print(run_command)
#Run the command on the terminal
os.system(run_command)
