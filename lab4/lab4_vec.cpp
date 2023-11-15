#include <iostream>
#include <vector>

#define MY_TEST
#include "myaloc.h"

template <typename T>
struct MyAllocator {
    using value_type = T;

    MyAllocator() noexcept {}

    template <typename U>
    MyAllocator(const MyAllocator<U>&) noexcept {}

    T* allocate(std::size_t n) {
        return static_cast<T*>(myalloc(n * sizeof(T)));
    }

    void deallocate(T* p, std::size_t n) {
        myfree(p);
    }
};

int main() {
    std::size_t num_elements = 10;

    const std::size_t memSize = 1024;

    void* memoryBlock = std::malloc(memSize);

    mysetup(memoryBlock, memSize);

    std::vector<int, MyAllocator<int>> vec;

    for(std::size_t i = 0; i < num_elements; ++i){
        vec.push_back((i + 20) * (10 + i));
    }    

    std::cout << "\nVec: ";
    for (const int& value : vec) {
        std::cout << value << " ";
    }
    std::cout << "\n\n";

    std::free(memoryBlock);
    return 0;
}

