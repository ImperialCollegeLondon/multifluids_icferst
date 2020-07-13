all: $(SRCS)


# Source Files
SRCS = frame.c  GID_B3D.c 

# Object Files
OBJS=$(SRCS:.c=.o)

# Executable Program Name
PROG=GID_B3DMC_Linux

# Compiler
CXX = gcc

# Include Directories
INCLUDES = -I. -I../include/ -I/usr/local/include

# Compiler Flags
CXXFLAGS = -O3 -ffast-math -funroll-loops -ftree-vectorize -fomit-frame-pointer


# Objects
%.o : %.c
	$(CXX) $(CXXFLAGS) -c -o $@ $^

# Library Paths
LIBS = -L../ -L/usr/local/lib64 -lm\
       -L/usr/X11/lib64 -lX11 -lXext -lXt\
       -L/usr/X11R6/lib64 -lGL -lGLU


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
