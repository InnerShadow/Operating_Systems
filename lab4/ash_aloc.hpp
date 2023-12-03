
// int test_1(void) {
//     const std::size_t memSize = 24;

//     void *memoryBlock = std::malloc(memSize);

//     mysetup(memoryBlock, memSize);

//     std::size_t *intPtr = static_cast<std::size_t *>(myalloc(sizeof(std::size_t)));

//     if (intPtr != nullptr) {
//         *intPtr = 42;
//         std::cout << "intPtr = " << *intPtr << "\n";
//         myfree(intPtr);
//     }

//     std::free(memoryBlock);

//     return 0;
// }