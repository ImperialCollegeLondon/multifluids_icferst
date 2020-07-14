all: $(SRCS)

# Object Files
OBJS=m2vtu.o frame.o  

# Executable Program Name
PROG=m2vtu

# Compiler
CC = gcc
CXX = g++

# Include Directories
INCLUDES = -I. -I../include/ -I/usr/local/include -I/usr/include/vtk-5.8 -I/usr/include/vtk-5.10/

# Compiler Flags
CFLAGS = -Wall -g $(INCLUDES)
CXXFLAGS = -Wall -g $(INCLUDES) -Wno-deprecated

# Library Paths
LIBS = -L../ -L/usr/lib/i -lm -lvtkCommon -lvtkGraphics \
-lvtkIO -lvtkFiltering -lvtkRendering -lvtkImaging
#  -lpolyobj
#       -L/usr/X11R6/lib -lGL -lGLU

VPATH = ./:./src

.SUFFIXES: .c .cpp .o

.c.o:
	$(CC) $(CFLAGS) -c $<

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

all: $(OBJS) $(PROG)

# Compilation command
$(PROG) : $(OBJS)
	$(CXX) $(INCLUDES) $(CXXFLAGS) $(OBJS) $(LIBS) -o $(PROG)

# Depend command
depend : $(SRCS)
	makedepend -- $(CXXFLAGS) -- $(INCLUDES) $(SRCS)

# Clean command
clean:
	rm -f $(PROG) $(OBJS) core

# DO NOT DELETE

