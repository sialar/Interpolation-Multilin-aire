SRC = src
INCLUDE = include
BIN = bin

CXX = g++
CXXFLAGS = -g -lm -std=c++11 -Wall

all: $(BIN)/test1D $(BIN)/test2D $(BIN)/testLejaSequence 

# Edition des liens et génération de l'exécutable:
$(BIN)/test1D : $(BIN)/test1D.o $(BIN)/LagrangeInterpolation1D.o $(BIN)/Utils.o
		$(CXX) $(CXXFLAGS) -o $@ $^
$(BIN)/test2D : $(BIN)/test2D.o $(BIN)/LagrangeInterpolation2D.o $(BIN)/Utils.o
		$(CXX) $(CXXFLAGS) -o $@ $^
$(BIN)/testND : $(BIN)/testND.o $(BIN)/LagrangeInterpolationND.o $(BIN)/Utils.o
		$(CXX) $(CXXFLAGS) -o $@ $^
$(BIN)/testLejaSequence : $(BIN)/testLejaSequence.o $(BIN)/Utils.o
		$(CXX) $(CXXFLAGS) -o $@ $^
# Construction des objets (file.o)
$(BIN)/test1D.o : $(SRC)/test1D.cpp $(INCLUDE)/LagrangeInterpolation1D.hpp
		$(CXX) $(CXXFLAGS) -c -o $@ $<
$(BIN)/test2D.o : $(SRC)/test2D.cpp $(INCLUDE)/LagrangeInterpolation2D.hpp
		$(CXX) $(CXXFLAGS) -c -o $@ $<
$(BIN)/testND.o : $(SRC)/testND.cpp $(INCLUDE)/LagrangeInterpolationND.hpp
		$(CXX) $(CXXFLAGS) -c -o $@ $<
$(BIN)/testLejaSequence.o : $(SRC)/testLejaSequence.cpp $(INCLUDE)/Utils.hpp
		$(CXX) $(CXXFLAGS) -c -o $@ $<
$(BIN)/LagrangeInterpolation%.o : $(SRC)/LagrangeInterpolation%.cpp $(INCLUDE)/LagrangeInterpolation%.hpp
		$(CXX) $(CXXFLAGS) -c -o $@ $<
$(BIN)/Utils.o : $(SRC)/Utils.cpp $(INCLUDE)/Utils.hpp
		$(CXX) $(CXXFLAGS) -c -o $@ $<

# Nettoyage du projet
clean:
	@rm -rf $(BIN)/*
