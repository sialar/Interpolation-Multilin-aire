SRC = src
INCLUDE = include
BIN = bin

CXX = g++
CXXFLAGS = -g -lm -std=c++11 -Wall

all: test1D test2D

test% : $(BIN)/test%.o $(BIN)/LagrangeInterpolation%.o $(BIN)/Utils.o
		$(CXX) $(CXXFLAGS) -o $@ $^

$(BIN)/test%.o : $(SRC)/test%.cpp $(INCLUDE)/LagrangeInterpolation%.hpp
		$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BIN)/LagrangeInterpolation%.o : $(SRC)/LagrangeInterpolation%.cpp $(INCLUDE)/LagrangeInterpolation%.hpp
		$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BIN)/Utils.o : $(SRC)/Utils.cpp $(INCLUDE)/Utils.hpp
		$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	@rm -rf test1D test2D $(BIN)/*.o
