.PHONY: all
.PHONY: clean

# Option
OP=0
SIZE=4

# Flags
CC=mpicc
CFLAGS=-O${OP} 
LDFLAGS=

# Files
BRC=$(shell grep -il "main" ./src/*.c)
TARGET=$(BRC:.c=)

all: $(TARGET) ex_4

hello_world:
	mpirun -hostfile hosts -n 10 ./src/$@

ex_4:
	mpirun -n 2 ./src/$@