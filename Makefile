SRC = src
INCLUDE = include
BIN = bin

CXX = g++
CXXFLAGS = -g -lm -std=c++11 -Wall

all: $(BIN)/TestLejaSequence $(BIN)/TestAlgoAI

# Edition des liens et génération de l'exécutable:
$(BIN)/TestAlgoAI : $(BIN)/TestAlgoAI.o $(BIN)/Interpolation.o $(BIN)/Utils.o $(BIN)/MultiVariatePoint.o
		$(CXX) $(CXXFLAGS) -o $@ $^
$(BIN)/TestLejaSequence : $(BIN)/TestLejaSequence.o $(BIN)/Utils.o $(BIN)/MultiVariatePoint.o
		$(CXX) $(CXXFLAGS) -o $@ $^

# Construction des objets (file.o)
$(BIN)/TestLejaSequence.o : $(SRC)/TestLejaSequence.cpp $(INCLUDE)/Utils.hpp
		$(CXX) $(CXXFLAGS) -c -o $@ $<
$(BIN)/TestAlgoAI.o : $(SRC)/TestAlgoAI.cpp $(INCLUDE)/Interpolation.hpp
		$(CXX) $(CXXFLAGS) -c -o $@ $<
$(BIN)/Interpolation.o : $(SRC)/Interpolation.cpp $(INCLUDE)/Interpolation.hpp
		$(CXX) $(CXXFLAGS) -c -o $@ $<
$(BIN)/MultiVariatePoint.o : $(SRC)/MultiVariatePoint.cpp $(INCLUDE)/MultiVariatePoint.hpp
		$(CXX) $(CXXFLAGS) -c -o $@ $<
$(BIN)/Utils.o : $(SRC)/Utils.cpp $(INCLUDE)/Utils.hpp
		$(CXX) $(CXXFLAGS) -c -o $@ $<

# Nettoyage du projet
clean:
	@rm -rf $(BIN)/*
