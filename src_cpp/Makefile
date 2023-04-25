# Makefile
CXX = g++
CXXFLAGS = -O2 -std=c++11 -fopenmp

TARGET = delta_stepping
OBJFILES = delta_stepping.o random_graph.o

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJFILES)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(TARGET) $(OBJFILES)
