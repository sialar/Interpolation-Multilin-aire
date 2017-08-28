SRC = AI/src
TEST = AI/tests
INCLUDE = AI/include
BIN = AI/bin

CXX = g++
CXXFLAGS = -O3 -lm -std=c++14 -g -Wall


OBJ_FILES = $(BIN)/LagrangeInterpolation.o $(BIN)/PiecewiseInterpolation.o $(BIN)/Functions.o
OBJ_FILES += $(BIN)/RealDataFunctions.o $(BIN)/AnalyticalFunctions.o
OBJ_FILES += $(BIN)/MixedInterpolation.o $(BIN)/Utils.o $(BIN)/BinaryTree.o
OBJ_FILES += $(BIN)/Tucker/LagrangePolynomial.o $(BIN)/Tucker/TuckerApproximation.o

INCLUDE_FILES = $(INCLUDE)/MultiVariatePoint.hpp $(INCLUDE)/Utils.hpp $(INCLUDE)/Functions.hpp
INCLUDE_FILES += $(INCLUDE)/Tucker/LagrangePolynomial.hpp $(INCLUDE)/Tucker/TuckerApproximation.hpp

TEST_FILES += $(BIN)/TestWithRealFunction $(BIN)/TestWithAnalyticalFunction

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
	@rm -rf $(BIN)/*.o $(BIN)/Test*
