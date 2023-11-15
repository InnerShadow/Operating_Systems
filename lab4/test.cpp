#include <iostream>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <fstream>

#define MY_TEST
const int PAGE_SIZE = 4096 * 1024;

#define TIME_COUNT

//#define SOMETIMES

#define MYALLOC

#include "pageing.h"

template <typename T>
struct MyAllocator {
    typedef T value_type;

    MemoryManager manager;

    MyAllocator() noexcept {
        initializeMemoryManager(&manager); 
    }

    template <class U>
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
        std::vector<int, MyAllocator<int>> vec;
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

    #endif

    #ifdef SOMETIMES
        }
    #endif

    return 0;
}

