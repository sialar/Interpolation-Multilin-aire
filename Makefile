SRC = src
TEST = tests
INCLUDE = include
BIN = bin

CXX = g++
CXXFLAGS = -O3 -lm -std=c++14 -g -Wall

OBJ_FILES = $(BIN)/LagrangeInterpolation.o $(BIN)/PiecewiseInterpolation.o
OBJ_FILES += $(BIN)/MixedInterpolation.o $(BIN)/Utils.o $(BIN)/BinaryTree.o $(BIN)/Functions.o
INCLUDE_FILES = $(INCLUDE)/MultiVariatePoint.hpp $(INCLUDE)/Utils.hpp
TEST_FILES = $(BIN)/TestLagrangeInterpolation $(BIN)/TestPiecewiseInterpolation
TEST_FILES += $(BIN)/TestMixedInterpolation $(BIN)/TestAutoMixedInterpolation
TEST_FILES += $(BIN)/TestErrors $(BIN)/TestSameFunctionWithDifferentPaths
TEST_FILES += $(BIN)/TestX

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
	@rm -rf $(BIN)/* build-window-Desktop_Qt_5_7_0_GCC_64bit-Debug/
