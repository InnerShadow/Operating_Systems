#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>

#ifdef MY_TEST
    #include <cstdlib>
    #include <iomanip>
    #include <iostream>

    #define MY_DEBUG
#endif

extern const int PAGE_SIZE;

typedef struct Page {
    size_t size;
    struct Page* next;
} Page;

typedef struct MemoryManager {
    Page* head;
} MemoryManager;

void initializeMemoryManager(MemoryManager* manager, void* initialMemory, size_t initialSize) {
    manager->head = (Page*)initialMemory;
    manager->head->size = initialSize - sizeof(Page);
    manager->head->next = NULL;
}

void* myalloc(MemoryManager* manager, size_t size) {
    #ifdef MY_DEBUG
        std::cout << "\n--------------------START MEMORY ALLOC--------------------\n\n";
    #endif

    if (manager->head == NULL || manager->head->size < size) {
        printf("Cannot allocate memory. Not enough space.\n");
        return NULL;
    }

    void* memory = (void*)((char*)manager->head + sizeof(Page));
    manager->head->size -= size;

    #ifdef MY_DEBUG
        std::cout << "\nAllocated " << size << " bytes of memory at address = " << memory << "\n";
    #endif

    #ifdef MY_DEBUG
        std::cout << "\n--------------------MEMORY ALLOC--------------------\n\n";
    #endif

    return memory;
}

void myfree(MemoryManager* manager, void* ptr) {
    #ifdef MY_DEBUG
        std::cout << "\n--------------------START MEMORY FREE--------------------\n\n";
    #endif

    if (ptr == NULL) {
        return;
    }

    #ifdef MY_DEBUG
        std::cout << "myfree(" << ptr << ')' << std::endl;
    #endif

    Page* currentPage = manager->head;
    while (currentPage != NULL) {
        char* pageStart = (char*)(currentPage + 1);
        char* pageEnd = pageStart + currentPage->size;

        if (ptr >= pageStart && ptr < pageEnd) {
            // Increase the size by the originally allocated size
            currentPage->size += PAGE_SIZE;

            #ifdef MY_DEBUG
                std::cout << "current page == " << currentPage << "\n";
            #endif

            #ifdef MY_DEBUG
                std::cout << "\n--------------------MEMORY FREE--------------------\n\n";
            #endif

            return;
        }
        currentPage = currentPage->next;
    }

    #ifdef MY_DEBUG
        std::cout << "\n--------------------CANNOT HANDLE MEMORY FREE--------------------\n\n";
    #endif
}

void freeMemoryManager(MemoryManager* manager) {
    try{
        while (manager->head != NULL) {
            Page* temp = manager->head;
            manager->head = manager->head->next;
            free(temp);
        }

        #ifdef MY_DEBUG
            std::cout << "\n--------------------FREE ALL ALLOCATED MEMORY--------------------\n\n";
        #endif

    } catch (const std::exception& e){
        return;
    }
}

