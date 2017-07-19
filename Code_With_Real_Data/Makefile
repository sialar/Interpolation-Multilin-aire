SRC = AI/src
TEST = AI/tests
INCLUDE = AI/include
BIN = AI/bin

CXX = g++
CXXFLAGS = -O3 -lm -std=c++14 -g -Wall

OUTPUT_FILES = AI/data/MOX/* AI/data/UOX/* AI/data/UOX-Gd/* AI/images/MOX/* AI/images/UOX/* AI/images/UOX-Gd/*
OBJ_FILES = $(BIN)/LagrangeInterpolation.o $(BIN)/PiecewiseInterpolation.o $(BIN)/Functions.o
OBJ_FILES += $(BIN)/MixedInterpolation.o $(BIN)/Utils.o $(BIN)/BinaryTree.o
INCLUDE_FILES = $(INCLUDE)/MultiVariatePoint.hpp $(INCLUDE)/Utils.hpp $(INCLUDE)/Functions.hpp
TEST_FILES += $(BIN)/TestAdaptativeInterpolation $(BIN)/Tests

all: $(TEST_FILES)

# Edition des liens et génération de l'exécutable:
$(BIN)/Test% : $(BIN)/Test%.o $(OBJ_FILES) $(INCLUDE_FILES)
		$(CXX) $(CXXFLAGS) -o $@ $^

# Construction des objets (file.o)
$(BIN)/Test%.o : $(TEST)/Test%.cpp $(INCLUDE_FILES)
		$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BIN)/%.o : $(SRC)/%.cpp $(INCLUDE_FILES)
		$(CXX) $(CXXFLAGS) -c -o $@ $<

# Nettoyage du projet
clean:
	@rm -rf $(BIN)/*.o $(BIN)/Tucker/* $(BIN)/Test* #$(OUTPUT_FILES)
