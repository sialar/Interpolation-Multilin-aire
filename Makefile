SRC = src
INCLUDE = include
BIN = bin

CXX = g++
CXXFLAGS = -O3 -lm -std=c++14 -Wall -g

OBJ_FILES = $(BIN)/LagrangeInterpolation.o $(BIN)/PiecewiseInterpolation.o
OBJ_FILES += $(BIN)/Utils.o $(BIN)/BinaryTree.o
INCLUDE_FILES = $(INCLUDE)/MultiVariatePoint.hpp $(INCLUDE)/Utils.hpp
TEST_FILES = $(BIN)/TestLagrangeInterpolation $(BIN)/TestPiecewiseInterpolation

all: $(TEST_FILES)

# Edition des liens et génération de l'exécutable:
$(BIN)/Test% : $(BIN)/Test%.o $(OBJ_FILES) $(INCLUDE_FILES)
		$(CXX) $(CXXFLAGS) -o $@ $^

# Construction des objets (file.o)
$(BIN)/Test%.o : $(SRC)/Test%.cpp $(INCLUDE_FILES)
		$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BIN)/%.o : $(SRC)/%.cpp $(INCLUDE_FILES)
		$(CXX) $(CXXFLAGS) -c -o $@ $<

# Nettoyage du projet
clean:
	@rm -rf $(BIN)/*
