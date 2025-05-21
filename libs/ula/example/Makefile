CC = g++
CCFLAGS  =  -O -Wall -ftemplate-depth-150
 
#LDFLAGS     = -L/home/itp/karen/scr/ITP/Lib/FFT/dfftpack
IDFLAGS     = -I. -I/home/karen/Work/Programming/Lib/ula/include -I/home/karen/Work/Programming/Lib/ula/include/ula
LIBS_LAPACK = -lgfortran -lm -llapack -lblas #-larpack   #-ldfftpack 
.cpp:
	$(CC) $(CCFLAGS) $@.cpp -o $@ $(LDFLAGS) $(IDFLAGS) $(LIBS_LAPACK) 

clean:
	\rm *~ 
