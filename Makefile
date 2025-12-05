BUILD_DIR=./build

all: debug

debug: clean_install cmake_debug install

release: clean_install cmake_release install

pristine: clean_all cmake_debug install

pristine_release: clean_all cmake_release install

run: run-parallel

run-serial: clean_out
	cd $(BUILD_DIR) && mpiexec -n 1 ./numsim ../parameters/parameters.txt

run-serial-3: clean_out
	cd $(BUILD_DIR) && mpiexec -n 2 ./numsim ../parameters/parameters-3.txt

run-parallel: clean_out
	cd $(BUILD_DIR) && mpiexec -n 2 ./numsim ../parameters/parameters.txt

run-parallel-3: clean_out
	cd $(BUILD_DIR) && mpiexec -n 2 ./numsim ../parameters/parameters-3.txt

cmake_debug:
	mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR) && cmake -DCMAKE_BUILD_TYPE=Debug ..

cmake_release:
	mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR) && cmake -DCMAKE_BUILD_TYPE=Release ..

install:
	cd $(BUILD_DIR) && make install

clean_install:
	@if [ -f $(BUILD_DIR)/Makefile ]; then cd $(BUILD_DIR) && make clean; fi

clean_out:
	rm -rf $(BUILD_DIR)/out

clean_all:
	rm -rf $(BUILD_DIR)
