SRC = src
INCLUDE = include
BIN = bin

CXX = g++
CXXFLAGS = -g -lm -std=c++11 -Wall

OBJ_FILES = $(BIN)/Interpolation.o $(BIN)/Utils.o
INCLUDE_FILES = $(INCLUDE)/MultiVariatePoint.hpp $(INCLUDE)/Utils.hpp
TEST_FILES = $(BIN)/TestLejaSequence $(BIN)/TestAlgoAI

all: $(TEST_FILES)

# Edition des liens et génération de l'exécutable:
$(BIN)/Test% : $(BIN)/Test%.o $(OBJ_FILES)
		$(CXX) $(CXXFLAGS) -o $@ $^

# Construction des objets (file.o)
$(BIN)/Test%.o : $(SRC)/Test%.cpp $(INCLUDE_FILES)
		$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BIN)/Interpolation.o : $(SRC)/Interpolation.cpp $(INCLUDE_FILES)
		$(CXX) $(CXXFLAGS) -c -o $@ $<
$(BIN)/Utils.o : $(SRC)/Utils.cpp $(INCLUDE_FILES)
		$(CXX) $(CXXFLAGS) -c -o $@ $<

# Nettoyage du projet
clean:
	@rm -rf $(BIN)/*
