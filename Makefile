CXXFLAGS = -O3 -march=native -DNDEBUG -std=c++14
CXX = g++
EIGEN_INCLUDE = ./eigen3

svd: main.cpp ImplicitQRSVD.h Tools.h
	$(CXX) $(CXXFLAGS) -o svd main.cpp -I$(EIGEN_INCLUDE)

clean:
	rm -f svd
