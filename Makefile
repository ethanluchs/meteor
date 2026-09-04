CC := gcc
CFLAGS := -Wall -Wextra -std=c11 -Iinclude
LDFLAGS := $(shell pkg-config --libs librtlsdr sqlite3 libcurl) -lm
CFLAGS += $(shell pkg-config --cflags librtlsdr sqlite3 libcurl)

SRC := $(wildcard src/*.c)
OBJ := $(patsubst src/%.c,build/%.o,$(SRC))
TARGET := build/meteor-detector

$(TARGET): $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDFLAGS)

build/%.o: src/%.c
	@mkdir -p build
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf build/*

.PHONY: clean
