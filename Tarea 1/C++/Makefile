IDIR = ./include
SDIR = ./src
ODIR = ./obj
LDIR = ./lib

CC = gcc
CXX = g++

CFLAGS = -I$(IDIR)
CPPFLAGS = -I$(IDIR)

EXEC = tarea_1

_DEPS = tarea_1.h 
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = tarea_1.o $(EXEC).o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

LIBS =

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CPPFLAGS)

$(EXEC): $(OBJ)
	g++ -o $@ $^ $(CPPFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o $(EXEC)