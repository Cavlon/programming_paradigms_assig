CC=gcc
DEPS=gs.h parse.h

all: gs.o parse.o main.c
	$(CC) -o main main.c gs.o parse.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $<

clean:
	rm -rf main.exe gs.o parse.o