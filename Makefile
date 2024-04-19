#
# Laboratory for Scientific Computing
# http://www.lam-mpi.org/tutorials/
# University of Notre Dame
#
# Sample Makefile
#
#
CXX       = mpic++
LIBS      = -lX11 -lm -fopenmp
CFLAGS    =
LDFLAGS   = -L/usr/X11R6/lib $(LIBS)
PROGRAM   = a.out                       #name of the binary
SRCS      = Life.cpp
default: all
all: $(PROGRAM) 
$(PROGRAM): $(SRCS)
        $(CXX) $(SRCS) $(CFLAGS) $(LDFLAGS)
clean:
        /bin/rm -f $(PROGRAM)
