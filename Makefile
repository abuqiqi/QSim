# ---------------------	macros --------------------------

# compiler and flags
CC := g++
CCFLAGS := -std=c++11 -Wall -Werror
INCFLAGS := -I.

COMPILE := $(CC) $(CCFLAGS) $(INCFLAGS) #-pthread 

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

# UTIL_CPPS := matrix.cpp
# UTIL_CPPS := $(wildcard $(UTIL_CPPS))

# UTIL_OBJS := $(addprefix $(OBJ_DIR)/, $(patsubst %.cpp,%.o,$(UTIL_CPPS)))

# # --------------------- main -----------------------

# MAIN_CPPS := qgate.cpp
# MAIN_CPPS := $(wildcard $(MAIN_CPPS))

# MAIN_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%,$(MAIN_CPPS))

# # -------------------- targets ---------------------

# TARGETS := $(MAIN_OBJS)

# .PHONY: all clean
# .PRECIOUS: %.o

# all: $(TARGETS)

# # utility functions: util/*.cpp -> obj/util/*.o
# # simulator functions: simulator/*.cpp -> obj/simulator/*.o
# $(OBJ_DIR)/%.o: %.cpp
# 	@echo "[INFO] Compiling" $< ...
# 	@$(COMPILE) -c $< -o $@

# # executable files: main/*.cpp -> obj/*
# $(OBJ_DIR)/%: $(OBJ_DIR) $(MAIN_CPPS) $(UTIL_OBJS)
# 	@echo "[INFO] Linking" $@ ...
# 	@$(COMPILE) $(patsubst $(OBJ_DIR)/%,$(SRC_DIR)/%.cpp,$@) $(UTIL_OBJS) -o $@
# 	@echo "[INFO]" $@ "has been built. "
