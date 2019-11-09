GOOGLE_TEST_LIB = gtest
GOOGLE_TEST_INCLUDE = /usr/local/include

G++ = g++
G++_FLAGS = -c -Wall -I $(GOOGLE_TEST_INCLUDE) -I ../ -I . -std=c++11
LD_FLAGS = -L /usr/local/lib -l $(GOOGLE_TEST_LIB) -l pthread

OBJECTS = main.o arithmetic.o
TARGET = arithmetic

all: $(TARGET)

$(TARGET): $(OBJECTS)
			g++ -o $(TARGET) $(OBJECTS) $(LD_FLAGS)

%.o : %.cpp
			$(G++) $(G++_FLAGS) $<

clean:
			rm -f $(TARGET) $(OBJECTS)
                    
.PHONY: all clean