These two scripts need to be modified by hand and then copied using sudo into /usr/bin. 
Note that the text between <> NEEDS to be modified!

RunICFERST.sh makes a copy of the ICFERST executable in the current directory, this is intended to ensure that if in the future the same test needs to be run a copy of the last executable is available.

Then You can easily run ICFERST and Diamond using the following commands:

Diamond.sh <INPUT>.mpml

RunICFERST.sh <INPUT>.mpml
