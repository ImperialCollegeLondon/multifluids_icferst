/*  Copyright (C) 2006 Imperial College London and others.
    
    Please see the AUTHORS file in the main source directory for a full list
    of copyright holders.

    Prof. C Pain
    Applied Modelling and Computation Group
    Department of Earth Science and Engineering
    Imperial College London

    amcgsoftware@imperial.ac.uk
    
    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation,
    version 2.1 of the License.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
    USA
*/

#include <cstring>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "confdefs.h"

#include "c++debug.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace std;

extern "C"{
  void checkmesh(const char *, size_t);
}

void Usage(){
  cout << "Usage: checkmesh MESH_BASENAME\n"
       << "\n"
       << "Checks the validity of the supplied mesh" << endl;
}

int main(int argc, char** argv){
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  // Undo some MPI init shenanigans
  int ierr = chdir(getenv("PWD"));
  if (ierr == -1) {
        cerr << "Unable to switch to directory " << getenv("PWD");
        abort();
  }
#endif
    
  if(argc < 2){
    Usage();
    return -1;
  }
  
  // Logging
  int nprocs = 1;
#ifdef HAVE_MPI
  int init_flag;
  MPI_Initialized(&init_flag);
  if(init_flag){
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  }
#endif
if(nprocs > 1){
  int rank = 0;
#ifdef HAVE_MPI
  int init_flag;
  MPI_Initialized(&init_flag);
  if (init_flag){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  }
#endif
  ostringstream buffer;
  buffer << "checkmesh.log-" << rank;
  if(freopen(buffer.str().c_str(), "w", stdout) == NULL){
    cerr << "Failed to redirect stdout" << endl;
    exit(-1);
  }
  buffer.str("");
  buffer << "checkmesh.err-" << rank;
  if(freopen(buffer.str().c_str(), "w", stderr) == NULL){
    cerr << "Failed to redirect stderr" << endl;
    exit(-1);
  }
  buffer.str("");
  }

  // Verbosity
  int verbosity = 2;
  set_global_debug_level_fc(&verbosity);

  string basename;
  basename.append(argv[1]);
  size_t basename_len = basename.size();

  checkmesh(basename.c_str(), basename_len);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
