CXX = g++
CXXFLAGS = -fopenmp
TARGETS = delta_stepping_parallel check dijkstra delta_stepping_serial

all: $(TARGETS)

delta_stepping_parallel: delta_stepping.cpp load_graph.o
	$(CXX) $(CXXFLAGS) -o $@ $< load_graph.o

check: check.cpp load_graph.o
	$(CXX) $(CXXFLAGS) -o $@ $< load_graph.o

dijkstra: dijkstra.cpp load_graph.o
	$(CXX) $(CXXFLAGS) -o $@ $< load_graph.o

delta_stepping_serial: delta_stepping_serial.cpp load_graph.o
	$(CXX) $(CXXFLAGS) -o $@ $< load_graph.o

load_graph.o: load_graph.cpp load_graph.h
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f delta_stepping dijkstra check load_graph delta_stepping_serial *.o
