CC=gcc
DEPS=linalg.h gs.h parse.h

all: linalg.o gs.o parse.o main.c
	$(CC) -o main main.c linalg.o gs.o parse.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $<

clean:
	rm -rf main.exe linalg.o gs.o parse.o