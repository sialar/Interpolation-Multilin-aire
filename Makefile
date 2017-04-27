SRC = src
INCLUDE = include
BIN = bin

CXX = g++
CXXFLAGS = -O3 -g -lm -std=gnu++11 -Wall
OPENMP_FLAGS = -fopenmp

OBJ_FILES = $(BIN)/Interpolation.o $(BIN)/Utils.o
INCLUDE_FILES = $(INCLUDE)/MultiVariatePoint.hpp $(INCLUDE)/Utils.hpp
TEST_FILES = $(BIN)/TestLejaSequence $(BIN)/TestAlgoAI

all: $(TEST_FILES)

# Edition des liens et génération de l'exécutable:
$(BIN)/Test% : $(BIN)/Test%.o $(OBJ_FILES)
		$(CXX) $(OPENMP_FLAGS) $(CXXFLAGS) -o $@ $^

# Construction des objets (file.o)
$(BIN)/Test%.o : $(SRC)/Test%.cpp $(INCLUDE_FILES)
		$(CXX) $(OPENMP_FLAGS) $(CXXFLAGS) -c -o $@ $<

$(BIN)/Interpolation.o : $(SRC)/Interpolation.cpp $(INCLUDE_FILES)
		$(CXX) $(OPENMP_FLAGS) $(CXXFLAGS) -c -o $@ $<
$(BIN)/Utils.o : $(SRC)/Utils.cpp $(INCLUDE_FILES)
		$(CXX) $(CXXFLAGS) -c -o $@ $<

# Nettoyage du projet
clean:
	@rm -rf $(BIN)/*
