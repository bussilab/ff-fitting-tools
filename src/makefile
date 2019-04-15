CC 				:=
	
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux)
		CC += g++
	endif
	ifeq ($(UNAME_S),Darwin)
		CC += g++-6
	endif

CFLAGS=-I.
OPT=-O2
MAIN=main.cpp
LEPTON=include/lepton
default: main 

System.o: System.cpp System.h
	$(CC) -O2 -c System.cpp -I./include -std=c++11

#Print.o: Print.cpp Print.h
#	$(CC) $(OPT) -c Print.cpp -I./include
ReadFromFile.o: System.o ReadFromFile.cpp ReadFromFile.h Print.h
	$(CC) -O2 -c ReadFromFile.cpp ./include/lepton/*.o -I./include -std=c++11

ObjFunction.o: ObjFunction.cpp ObjFunction.h 
	$(CC) -O2 -c ObjFunction.cpp -I./include -std=c++11

main: $(MAIN) ReadFromFile.o ObjFunction.o
	$(CC) -O2 -o fit.x $(MAIN) System.cpp ReadFromFile.cpp ObjFunction.cpp ./include/lepton/*.o -I./include -std=c++11 -fopenmp 

.PHONY: clean
clean: 
	rm -f out.x *.o $(LEPTON)/*.o


