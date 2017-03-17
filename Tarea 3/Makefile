all:test

test: obj/main.o 
	g++ obj/main.o -o test

main.o: main.cpp
	g++ -c src/main.cpp

clean:
	rm -rf *.o test
