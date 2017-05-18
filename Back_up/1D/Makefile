SRC = src
INCLUDE = include
BIN = bin

CXX = g++
CXXFLAGS = -g -lm -std=c++11 -Wall

all: $(BIN)/test1D

# Edition des liens et génération de l'exécutable:
$(BIN)/test1D : $(BIN)/test1D.o $(BIN)/Interpolation.o $(BIN)/BinaryTree.o
		$(CXX) $(CXXFLAGS) -o $@ $^

# Construction des objets (file.o)
$(BIN)/test1D.o : $(SRC)/test1D.cpp $(INCLUDE)/Interpolation.hpp
		$(CXX) $(CXXFLAGS) -c -o $@ $<
$(BIN)/Interpolation.o : $(SRC)/Interpolation.cpp $(INCLUDE)/Interpolation.hpp
		$(CXX) $(CXXFLAGS) -c -o $@ $<
$(BIN)/BinaryTree.o : $(SRC)/BinaryTree.cpp $(INCLUDE)/BinaryTree.hpp
		$(CXX) $(CXXFLAGS) -c -o $@ $<


# Nettoyage du projet
clean:
	@rm -rf $(BIN)/*
