CC=g++
DEPS=linalg.h gs.h parse.h

all: linalg.o gs.o parse.o main.cpp
	$(CC) -o main main.cpp linalg.o gs.o parse.o

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $<

clean:
	rm -rf main.exe linalg.o gs.o parse.o