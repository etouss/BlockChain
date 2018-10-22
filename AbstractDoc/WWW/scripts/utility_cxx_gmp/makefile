CC      = g++
CFLAGS  = -Wall -O3 -std=c++0x -lgmp -lgmpxx
RM      = rm -f

default: af wf

af: af.cpp tools.cpp
	$(CC) af.cpp tools.cpp -o af.o $(CFLAGS)

wf: window_fork.cpp tools.cpp
	$(CC) window_fork.cpp tools.cpp -o wf.o $(CFLAGS)

clean: 
	rm *.o 
	rm *.temp
