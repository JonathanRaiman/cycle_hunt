all:
	clang++ --std=c++11 -stdlib=libc++ cycle_removal.cpp -o cycle_removal
clean:
	rm -f cycle_removal
