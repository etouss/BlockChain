CC      = g++
CFLAGS  = -Wall -O3 -std=c++0x -lgmp -lgmpxx
RM      = rm -f

default: btc

btc: btc.cpp
	$(CC) btc.cpp -o btc.o $(CFLAGS)

clean: 
	rm *.o
