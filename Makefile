SRC = src
INCLUDE = include
BIN = bin

CXX = g++
CXXFLAGS = -O3 -lm -std=c++14 -Wall -g

OBJ_FILES = $(BIN)/Interpolation.o $(BIN)/Utils.o $(BIN)/BinaryTree.o
OBJ_FILES += $(BIN)/LagrangeInterpolation.o $(BIN)/PiecewiseInterpolation.o
INCLUDE_FILES = $(INCLUDE)/MultiVariatePoint.hpp $(INCLUDE)/Utils.hpp
TEST_FILES = $(BIN)/TestLejaSequence $(BIN)/TestLagrangeInterpolation
TEST_FILES += $(BIN)/TestPiecewiseInterpolation

all: $(TEST_FILES)

# Edition des liens et génération de l'exécutable:
$(BIN)/Test% : $(BIN)/Test%.o $(OBJ_FILES)
		$(CXX) $(CXXFLAGS) -o $@ $^

# Construction des objets (file.o)
$(BIN)/Test%.o : $(SRC)/Test%.cpp $(INCLUDE_FILES)
		$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BIN)/%.o : $(SRC)/%.cpp $(INCLUDE_FILES)
		$(CXX) $(CXXFLAGS) -c -o $@ $<

# Nettoyage du projet
clean:
	@rm -rf $(BIN)/*
