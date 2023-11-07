#include <iostream>

#define MY_TEST
#include "myaloc.h"

int main() {
    const std::size_t memSize = 1024;

    void* memoryBlock = std::malloc(memSize);

    mysetup(memoryBlock, memSize);

    std::size_t* intPtr = static_cast<std::size_t*>(myalloc(sizeof(std::size_t)));

    *intPtr = 42;
    std::cout << "\nintPtr = " << *intPtr << "\n\n";

    myfree(intPtr);

    std::free(memoryBlock);

    return 0;
}

