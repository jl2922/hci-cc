CONFIG_FILE := Makefile.config
ifeq ($(wildcard $(CONFIG_FILE)),)
$(error $(CONFIG_FILE) not found.)
endif
include $(CONFIG_FILE)

ifeq ($(wildcard $(BOOST_DIR)),)
$(error Boost folder $(BOOST_DIR) not found.)
endif

CC := g++
CFLAGS := -std=c++11 -O3 -I $(BOOST_DIR)
SRC_DIR := src
OBJ_DIR := build
SRCS := $(shell find $(SRC_DIR)/ ! -name "main.cc" -name "*.cc")
OBJS := $(SRCS:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)
MAIN := src/main.cc
EXE := hci

.PHONY: all clean

all: $(EXE)

clean:
	rm -f $(OBJ_DIR)/*

$(EXE): $(OBJS) $(MAIN)
	$(CC) $(CFLAGS) $(MAIN) $(OBJS) -o hci
	
$(OBJS): $(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc $(SRC_DIR)/%.h
	$(CC) $(CFLAGS) -c $< -o $@
