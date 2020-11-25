EXE =	chromo

SRC = \
	initConfig.cpp \
	dynamics.cpp \
	chromoCell_main.cpp

INC = \
	initConfig.hpp \
	dynamics.hpp

OBJ = \
	initConfig.o \
	dynamics.o \
	chromoCell_main.o	

# if compiling for serial or openmp code
CXX = g++
# if compiling MPI c++ code for single machine (mpiexec from mpich) 
#CXX = mpic++

DEBUG = -g
CFLAGS = -fopenmp -O3 -std=c++11 $(DEBUG)
LFLAGS = -fopenmp

all : $(EXE)

$(EXE) : $(OBJ)
	$(CXX) -o $@  $(LFLAGS) $(OBJ)

#$(OBJ) : $(INC)
#	$(CXX) $(CFLAGS) -c $(SRC)

initConfig.o : initConfig.hpp initConfig.cpp
	$(CXX) $(CFLAGS) -c initConfig.cpp	
dynamics.o : initConfig.hpp dynamics.hpp dynamics.cpp
	$(CXX) $(CFLAGS) -c dynamics.cpp	
chromoCell_main.o : initConfig.hpp dynamics.hpp chromoCell_main.cpp
	$(CXX) $(CFLAGS) -c chromoCell_main.cpp

clean:
	\rm $(OBJ) $(EXE)
