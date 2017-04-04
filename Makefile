CONFIG_FILE := Makefile.config
ifeq ($(wildcard $(CONFIG_FILE)),)
$(error $(CONFIG_FILE) not found.)
endif
include $(CONFIG_FILE)

ifeq ($(wildcard $(BOOST_DIR)),)
$(error Boost $(BOOST_DIR) not found.)
endif

CC := mpiCC
CFLAGS := -std=c++11 -g -Wall -O3 -lboost_mpi # -I $(BOOST_DIR)
SRC_DIR := src
OBJ_DIR := build
SRCS := $(shell find $(SRC_DIR) ! -name "main.cc" -name "*.cc")
HEADERS := $(shell find $(SRC_DIR) -name "*.h")
OBJS := $(SRCS:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)
MAIN := src/main.cc
EXE := hci

ifeq ($(wildcard $(OBJ_DIR)),)
$(shell mkdir $(OBJ_DIR))
endif

.PHONY: all clean

all: $(EXE)

clean:
	rm -f $(OBJ_DIR)/*

$(EXE): $(OBJS) $(MAIN)
	$(CC) $(CFLAGS) $(MAIN) $(OBJS) -o hci
	
$(OBJS): $(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@
