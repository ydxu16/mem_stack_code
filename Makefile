CXX = icc
#CXX = gcc
LIBS = -lfftw3 -lm # -liomp5 -lpthread
#CXXFLAGS
CXXFLAGS = -lstdc++ -g -Wall -O3 -openmp -I/opt/intel/composer_xe_2013_sp1.4.211/compiler/include/ #-L/opt/intel/composer_xe_2013_sp1.4.211/compiler/lib/intel64/

PROG_SOLVE = solver 
SOLVE_objects = solver.o OneLayer.o FFTW_Mat_C_2D.o FFTW_Mat_R_2D.o FFTW_Mat_R_3D.o functions.o


all: $(PROG_SOLVE) 

$(PROG_SOLVE) : $(SOLVE_objects)
	$(CXX) ${CXXFLAGS} -o $@ $^ -L/opt/fftw3-gcc/lib $(LIBS) 

$(SOLVE_objects) : 
	$(CXX) ${CXXFLAGS} -c $(SOLVE_objects:.o=.cc) -I/opt/fftw3-gcc/include 

clean:
	$(RM) *.o *~
	$(RM) .depend

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend
