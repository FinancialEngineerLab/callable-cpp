CXX            := -c++
CXXFLAGS       := -O2 # -pedantic-errors -Wall -Wextra -Werror
LDFLAGS        := -lgtest -lpthread
INCLUDE        := -Iinclude/ -IEigen/ -I/usr/src/gtest/include

BUILD          := ./build
OBJ_DIR        := $(BUILD)/objects
MAIN_DIR       := $(BUILD)/main
TEST_DIR       := $(BUILD)/tests
TARGET         := beagle
TEST_TARGET    := run_tests

MAIN_CPP_NAME  := beagle.cpp
SRC_DIR        := src
SRC            := $(wildcard $(SRC_DIR)/*.cpp)
NON_MAIN_SRC   := $(filter-out $(SRC_DIR)/$(MAIN_CPP_NAME), $(SRC))
OBJS           := $(SRC:%.cpp=$(OBJ_DIR)/%.o)

TEST_SRC_DIR   := test
TEST_SRC       := $(wildcard $(TEST_SRC_DIR)/*.cpp $(NON_MAIN_SRC))
TEST_OBJS      := $(TEST_SRC:%.cpp=$(OBJ_DIR)/%.o)

beagle: build $(MAIN_DIR)/$(TARGET)
test:   build $(TEST_DIR)/$(TEST_TARGET)
all: beagle test

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@ $(LDFLAGS)

$(MAIN_DIR)/$(TARGET): $(OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $(MAIN_DIR)/$(TARGET) $^ $(LDFLAGS)

$(TEST_DIR)/$(TEST_TARGET): $(TEST_OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $(TEST_DIR)/$(TEST_TARGET) $^ $(LDFLAGS)

.PHONY: all build clean debug release

build:
	@mkdir -p $(MAIN_DIR)
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(TEST_DIR)

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(MAIN_DIR)/*
	-@rm -rvf $(TEST_DIR)/*
