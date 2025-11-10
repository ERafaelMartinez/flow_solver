BUILD_DIR=./build

all: debug

debug: clean_install cmake_debug install

release: clean_install cmake_release install

prestine_debug: clean_all cmake_debug install

prestine_release: clean_all cmake_release install

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

clean_all:
	rm -rf $(BUILD_DIR)
