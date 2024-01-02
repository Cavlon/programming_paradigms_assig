CC=g++
DEPS=linalg.h parse.h svp.h

all: linalg.o parse.o svp.o main.cpp
	$(CC) -o runme.exe main.cpp linalg.o parse.o svp.o

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $<

clean:
	rm -rf runme.exe test.exe linalg.o parse.o svp.o result.txt

test: linalg.o parse.o svp.o test.cpp
	$(CC) -o test.exe test.cpp linalg.o parse.o svp.o