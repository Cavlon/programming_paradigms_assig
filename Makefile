CC=g++
DEPS=linalg.h parse.h

all: linalg.o parse.o main.cpp
	$(CC) -o runme main.cpp linalg.o parse.o

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $<

clean:
	rm -rf runme.exe linalg.o parse.o