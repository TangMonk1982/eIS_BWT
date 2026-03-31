all:build_bwt

build_bwt:
	/opt/intel/oneapi/compiler/2021.2.0/linux/bin/intel64/icpc src/build.cpp -std=c++14 -DNDEBUG -march=native -O3 -o build_bwt -pthread -DCILK
clean:
	rm build_bwt
