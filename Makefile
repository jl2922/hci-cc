CONFIG_FILE := Makefile.config
ifeq ($(wildcard $(CONFIG_FILE)),)
$(error $(CONFIG_FILE) not found.)
endif
include $(CONFIG_FILE)

ifeq ($(wildcard $(BOOST_DIR)),)
$(error Boost folder $(BOOST_DIR) not found.)
endif

CC := g++
CFLAGS := -std=c++11 -O3 -Wall -I $(BOOST_DIR)
SRC_DIR := src
OBJ_DIR := build
SRCS := $(shell find $(SRC_DIR)/ -name "*.cc")
OBJS := $(SRCS:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)

.PHONY: all clean

all: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o hci

clean:
	rm -f $(OBJ_DIR)/*

$(OBJS): $(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	$(CC) $(CFLAGS) -c $< -o $@
