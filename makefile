GOOGLE_TEST_LIB = gtest
GOOGLE_TEST_INCLUDE = /usr/src/gtest/include

CXX       := -c++
CXXFLAGS := -O2 -I $(GOOGLE_TEST_INCLUDE) # -pedantic-errors -Wall -Wextra -Werror
LDFLAGS  := -L/usr/lib/ -l$(GOOGLE_TEST_LIB) -lpthread
BUILD    := ./build
OBJ_DIR  := $(BUILD)/objects
APP_DIR  := $(BUILD)/apps
TEST_DIR := $(BUILD)/tests
TARGET   := program
TESTS    := run_tests
INCLUDE  := -Iinclude/ -IEigen/
SRC      :=                      \
   $(wildcard src/*.cpp test/*.cpp)         \

OBJECTS  := $(SRC:%.cpp=$(OBJ_DIR)/%.o)

test: build $(TEST_DIR)/$(TESTS)

all: build $(APP_DIR)/$(TARGET) \
	   build $(TEST_DIR)/$(TESTS)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@ $(LDFLAGS)

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $(APP_DIR)/$(TARGET) $^ $(LDFLAGS)

$(TEST_DIR)/$(TESTS): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $(TEST_DIR)/$(TESTS) $^ $(LDFLAGS)

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
