CPP=g++ 
CC=gcc
OPTIMISER = -O3
CPPFLAGS =  -Wall -Wmissing-prototypes -Wshadow -fmessage-length=0 -std=c++11 -msse2 -mfpmath=sse
CFLAGS = 

INC = -I/usr/local/include
PROGRAM = msascorer

# Headers
HDR = Sequence.h 

# Source
CPPS = MSAscorer.cpp Sequence.cpp 
CPPO = MSAscorer.o Sequence.o 

all : $(PROGRAM)

$(CPPO) : $(CPPS) $(HDR)
	$(CPP) $(OPTIMISER) $(CPPFLAGS) $(INC) -c $(CPPS)

$(PROGRAM) : $(CPPO) 
	$(CPP) $(CPPFLAGS) $(OPTIMISER) $(INC) $(LIB) $(CPPO) -o $(PROGRAM)


clean:
	rm -f $(CPPO)
	rm -f $(PROGRAM)
	rm -f core

