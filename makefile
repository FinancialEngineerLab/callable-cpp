CXX            := -c++
CXXFLAGS       := -O2 # -pedantic-errors -Wall -Wextra -Werror
LDFLAGS        := -lgtest -lpthread
INCLUDE        := -Iinclude/ -IEigen/ -I/usr/src/gtest/include

BUILD          := ./build
OBJ_DIR        := $(BUILD)/objects
APP_DIR        := $(BUILD)/apps
TEST_DIR       := $(BUILD)/tests
TARGET         := program
TEST_TARGET    := run_tests

MAIN_CPP_NAME  := test.cpp
SRC_DIR        := src
SRC            := $(wildcard $(SRC_DIR)/*.cpp)
NON_MAIN_SRC   := $(filter-out $(SRC_DIR)/$(MAIN_CPP_NAME), $(SRC))
OBJS           := $(SRC:%.cpp=$(OBJ_DIR)/%.o)

TEST_SRC_DIR   := test
TEST_DIR_SRC   := $(wildcard $(TEST_SRC_DIR)/*.cpp)
TEST_SRC       := $(wildcard $(NON_MAIN_SRC) $(TEST_DIR_SRC)/*.cpp)
TEST_OBJS      := $(TEST_SRC:%.cpp=$(OBJ_DIR)/%.o)


all: build $(APP_DIR)/$(TARGET)       \
		 build $(TEST_DIR)/$(TEST_TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@ $(LDFLAGS)

$(APP_DIR)/$(TARGET): $(OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $(APP_DIR)/$(TARGET) $^ $(LDFLAGS)

$(TEST_DIR)/$(TEST_TARGET): $(TEST_OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $(TEST_DIR)/$(TEST_TARGET) $^ $(LDFLAGS)

.PHONY: all build clean debug release

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(TEST_DIR)

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O2
release: all

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*
	-@rm -rvf $(TEST_DIR)/*
