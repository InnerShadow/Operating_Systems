#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <fstream>

const int PAGE_SIZE = 4096;

//#define TIME_COUNT

#define MYALLOC

//#define SOMETIMES

#define MY_TEST

#include "pageingPreInit.h"

template <typename T>
struct MyAllocator {
    using value_type = T;

    MemoryManager manager;

    MyAllocator(void* initialMemory, size_t initialSize) noexcept {
        initializeMemoryManager(&manager, initialMemory, initialSize);
    }

    template <typename U>
    MyAllocator(const MyAllocator<U>&) noexcept {}

    T* allocate(std::size_t n) {
        return static_cast<T*>(myalloc(&manager, n * sizeof(T)));
    }

    void deallocate(T* p, std::size_t n) {
        myfree(&manager, static_cast<void*>(p));
    }

    ~MyAllocator() {
        freeMemoryManager(&manager);
    }
};

int main(int argc, char* argv[]) {

    #ifdef SOMETIMES
        for(std::size_t j = 0; j < 100; ++j){
    #endif

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <number_of_elements>" << std::endl;
        return 1;
    }

    int numElements = std::atoi(argv[1]);

    #ifdef TIME_COUNT
        auto start = std::chrono::high_resolution_clock::now();
    #endif

    #ifdef MYALLOC
        size_t blockSize = PAGE_SIZE * 10;
        void* memoryBlock = std::malloc(blockSize);

        MyAllocator<int> allocator(memoryBlock, blockSize);

        std::vector<int, MyAllocator<int>> vec(allocator);
    #endif

    #ifndef MYALLOC
        std::vector<int> vec;
    #endif

    for (int i = 1; i <= numElements; ++i) {
        vec.push_back(i * 10);
    }

    #ifdef TIME_COUNT
        auto end = std::chrono::high_resolution_clock::now();
    #endif

    std::cout << "vec: ";
    for (const auto& element : vec) {
        std::cout << element << " ";
    }

    std::cout << "\n\n";

    #ifdef TIME_COUNT
        auto elapsedTime = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        std::cout << "Time taken: " << elapsedTime << " microseconds" << std::endl;

        #ifdef MYALLOC
            std::ofstream outFile("time_" + std::to_string(numElements) + ".txt", std::ios::app);
        #endif

        #ifndef MYALLOC
            std::ofstream outFile("CPP_time_" + std::to_string(numElements) + ".txt", std::ios::app);
        #endif

        outFile << elapsedTime << "\n";

        std::free(memoryBlock);
    #endif

    #ifdef SOMETIMES
        }
    #endif

    return 0;
}

