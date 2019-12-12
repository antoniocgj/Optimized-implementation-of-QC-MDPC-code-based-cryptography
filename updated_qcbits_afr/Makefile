CC=gcc
FLAGS= -O3 -march=native -funroll-all-loops -mpopcnt

ifdef UNIFORM
FLAGS += -DUNIFORM
endif

ifdef OPENSSL
FLAGS += -DOPENSSL
endif

all: qcbits.c mainTest.c 
	$(CC) qcbits.c mainTest.c $(FLAGS) -o qcbits

clean:
	rm qcbits

