CC=g++
ARMS=arms/arms.o
TAJD=tajd.o
HOME=/home/ss545
CFLAGS= -c -g -I /opt/boost-1.57.0/include/ -I ${HOME}/bin/alglib/cpp/src/  -I ${HOME}/bin/gzstream/ -I /opt/nag/libC/include/ -Wno-write-strings -Wpointer-arith -pg -std=c++0x
#CFLAGS= -c -g -I ./include -I ${HOME}/bin/alglib/cpp/src/  -I ${HOME}/bin/gzstream/ -I /opt/nag/libC/include/  -Wno-write-strings -Wpointer-arith -std=c++0x -pg 
#LIB= -L ${HOME}/bin/alglib/cpp/src -lalglib  -lgsl -lgslcblas -lm -lz -L ${HOME}/bin/gzstream -lgzstream -L /opt/nag/libC/lib -lnagc_nag
LIB= -L ${HOME}/bin/alglib/cpp/src -lalglib -L /opt/boost-1.57.0/lib -lgsl -lgslcblas -lm -lboost_filesystem -lboost_system -lz -L ${HOME}/bin/gzstream -lgzstream -L /opt/nag/libC/lib -lnagc_nag
SRC1 =  functions.cpp io.cpp data.cpp snp.cpp ind.cpp packedgtype.cpp interval.cpp intervalmap.cpp genotype.cpp fileio.cpp computed.cpp 
SRC12 =  functions.cpp io.cpp data.cpp snp.cpp ind.cpp packedgtype.cpp interval.cpp intervalmap.cpp genotype.cpp snpmap.cpp fileio.cpp getposterior.cpp
HEADER = samplealpha.h


OBJ1=$(SRC1:.cpp=.o)
OBJ12=$(SRC12:.cpp=.o)

EXEC=computed getposterior

#$(warning OBJ=$(OBJ))


.SUFFIXES:
.SUFFIXES:  .cpp .o

.cpp.o:
	$(CC) $(CFLAGS) -o $@ -c $<

all:$(EXEC) 

clean: 
	rm -f $(EXEC) *.o

computed: $(OBJ1)
	$(CC) -o $@ $(OBJ1)  $(LIB) 	

getposterior: $(OBJ12)
	$(CC) -o $@ $(OBJ12) $(ARMS)  $(LIB) 	

install:
	cp computed getposterior ../bin
