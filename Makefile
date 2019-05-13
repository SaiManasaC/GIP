cpp_source=compress.cpp verifier.cpp aligner.cpp utils.cpp buildindex.cpp refcom.cpp
src_dir=src
objs_dir=objs
objs+= $(patsubst %.cpp, $(objs_dir)/%.o, $(cpp_source))

cxx=g++
cxxflags=-std=c++11 -Wall -O3 -funroll-all-loops -fopenmp -march=native -lpthread # -Wextra -pg
#cxxflags=-std=c++11 -Wall -O3 -funroll-all-loops -fopenmp -march=native -lpthread -DNDEBUG

ldflags= 
exec=REFCOM

all: dir $(exec)
	
dir:
	-mkdir -p $(objs_dir)

$(exec): $(objs)
	$(cxx) $(cxxflags) $(objs) -o $(exec)
	
clean:
	-rm -rf $(objs)

$(objs_dir)/%.o: $(src_dir)/%.cpp
	$(cxx) $(cxxflags) -o $@ -c $<

