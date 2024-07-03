# ---------------------	macros --------------------------

# compiler and flags
CC := g++
CCFLAGS := -std=c++11 -Wall -Werror
INCFLAGS := -I.

COMPILE := $(CC) $(CCFLAGS) $(INCFLAGS)

# --------------------- utils ----------------------

TARGET := main

all: $(TARGET)

$(TARGET): main.o qgate.o matrix.o qcircuit.o
	$(COMPILE) -o $@ $^

%.o: %.cpp %.h
	$(COMPILE) -c $< -o $@

clean:
	del *.o *.exe

.PHONY: all clean
