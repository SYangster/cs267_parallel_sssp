CXX = g++
CXXFLAGS = -std=c++11 -O3 -fopenmp
LDFLAGS = -L/opt/cray/pe/gcc/11.2.0/snos/lib64 -L/opt/cray/pe/libsci/22.11.1.2/GNU/9.1/x86_64/lib
LDLIBS = -lsci_gnu_82_mpi_mp -lsci_gnu_82_mp -lgomp
TARGETS = delta_stepping_parallel check dijkstra delta_stepping_serial radius_stepping_serial

all: $(TARGETS)

delta_stepping_parallel: delta_stepping.cpp load_graph.o
	$(CXX) $(CXXFLAGS) -o $@ $< load_graph.o

check: check.cpp load_graph.o
	$(CXX) $(CXXFLAGS) -o $@ $< load_graph.o

dijkstra: dijkstra.cpp load_graph.o
	$(CXX) $(CXXFLAGS) -o $@ $< load_graph.o

delta_stepping_serial: delta_stepping_serial.cpp load_graph.o
	$(CXX) $(CXXFLAGS) -o $@ $< load_graph.o

radius_stepping_serial: radius_stepping_serial.cpp load_graph.o
	$(CXX) $(CXXFLAGS) -o $@ $< load_graph.o

load_graph.o: load_graph.cpp load_graph.h
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f delta_stepping dijkstra check load_graph delta_stepping_serial radius_stepping_serial *.o
