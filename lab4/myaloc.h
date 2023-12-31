#include <iostream>
#include <climits>

#ifdef MY_TEST
  #include <cstdlib>
  #include <iomanip>
  
  #define MY_DEBUG
#endif

#define MIN_BLOCK_SIZE 16

typedef struct header* (*HeaderIterator)(struct header*);

struct header {
    unsigned char free;
    std::size_t actualSize;
};

struct tail {
    unsigned char free; 
    std::size_t actualSize;
};

struct header* _startHeader;
struct tail* _endTail;
std::size_t _size;
std::size_t minSize;

struct header* initBlock(void *buf, std::size_t size) {
    #ifdef MY_DEBUG
        std::cout << "initBlock(" << buf << ", " << size << ')' << "\n";
    #endif

    struct header* ph = (struct header*)buf;
    ph->free = 1;
    ph->actualSize = size - sizeof(struct header) - sizeof(struct tail);

    #ifdef MY_DEBUG
        std::cout << "   ph addr = " << ph << "; free = " << (int)(ph->free) << "; actualSize = " << ph->actualSize << "\n";
    #endif

    struct tail* pt = (struct tail*)((unsigned char*)buf + size - sizeof(struct tail));
    pt->free = 1;
    pt->actualSize = ph->actualSize;

    #ifdef MY_DEBUG
        std::cout << "   pt addr = " << pt << "; free = " << (int)(pt->free) << "; actualSize = " << pt->actualSize << "\n";
    #endif

    return ph;
}

struct tail* getTail(struct header* phead) {
    unsigned char* base = (unsigned char*)phead;
    base += sizeof(struct header) + phead->actualSize;
    return (struct tail*)base;
}

struct header* getHeader(struct tail* ptail) {
    unsigned char* base = (unsigned char*)ptail;
    base -= sizeof(struct header) + ptail->actualSize;
    return (struct header*)base;
}

struct header* getNext(struct header* phead) {
    #ifdef MY_DEBUG
        std::cout << "getNext(" << phead << ')' << "\n";
    #endif

    struct tail* ptail = getTail(phead);
    if(ptail == _endTail){
        return NULL;
    }

    unsigned char* base = (unsigned char*)ptail;
    base += sizeof(struct tail);
    return (struct header*)base;
}

struct header* getPrevious(struct header* phead) {
    #ifdef MY_DEBUG
        std::cout << "getPrevious(" << phead << ')' << "\n";
    #endif

    if(phead == _startHeader){
        return NULL;
    }

    unsigned char* base = (unsigned char*)phead;
    base -= sizeof(struct tail);
    return getHeader((struct tail*)base);
}

std::size_t getAllSize(std::size_t size) {
    return size + sizeof(struct header) + sizeof(struct tail);
}

std::size_t getActualSize(struct header* ph, struct tail* pt) {
    if((unsigned char*)ph >= (unsigned char*)pt){
        return 0;
    }

    return (unsigned char*)pt - (unsigned char*)ph - sizeof(struct header);
}

void joinBlocks(struct header* phStart, struct header* phEnd) {
    if(phStart == phEnd){
        return;
    }
  
    #ifdef MY_DEBUG
        std::cout << "joinBlocks(" << phStart << ", " << phEnd << ')' << "\n";
    #endif
  
    if(phStart > phEnd) {
        struct header* tmp = phStart;
        phStart = phEnd;
        phEnd = tmp;
        #ifdef MY_DEBUG
            std::cout << "  swap block. phtStart " << phStart << "; phEnd " << phEnd << "\n";
        #endif
    }
  
    struct tail* ptEnd = getTail(phEnd);
    #ifdef MY_DEBUG
        std::cout << "   ptEnd = " << ptEnd << "\n";
    #endif
  
    phStart->actualSize = getActualSize(phStart, ptEnd);
    ptEnd->actualSize = phStart->actualSize;

    #ifdef MY_DEBUG
        std::cout << "   actualSize = " << phStart->actualSize << "\n";
    #endif
}

struct header* utilizeBlock(struct header* ph, std::size_t size) {
    if(!ph->free || ph->actualSize < size){
        return NULL;
    }
  
    #ifdef MY_DEBUG
        std::cout << "utilizeBlock(" << ph << ", " << size << ")" << "\n";
    #endif
  
    std::size_t allSize = getAllSize(size);
  
    struct tail* ptEnd = getTail(ph);
  
    if(ph->actualSize <= allSize) {
        ph->free = 0;
        ptEnd->free = 0;
        return ph;
    }
  
    std::size_t newSize = ph->actualSize - allSize;
  
    if(newSize < MIN_BLOCK_SIZE) {
        ph->free = 0;
        ptEnd->free = 0;
        return ph;
    }
  
    unsigned char* base = (unsigned char*)ptEnd;
    base -= allSize;
  
    struct tail* ptStart = (struct tail*)base;
    ptStart->free = 1;
    ptStart->actualSize = getActualSize(ph, ptStart);
    ph->actualSize = ptStart->actualSize;
  
    struct header* phEnd = (struct header*)(base + sizeof(struct tail));
    phEnd->free = 0;
    phEnd->actualSize = getActualSize(phEnd, ptEnd);
    
    ptEnd->free = 0;
    ptEnd->actualSize = phEnd->actualSize;
    
    return phEnd;
}

struct header* joinNearestFreeBlocks(struct header* ph, HeaderIterator iterator) {
    if(ph == NULL || !ph->free)
        return ph;
  
    #ifdef MY_DEBUG
        std::cout << "joinNearestFreeBlocks(" << ph << ')' << "\n";
    #endif
  
    struct header* next = iterator(ph);
    #ifdef MY_DEBUG
        std::cout << "   next header = " << next << " " << (next == NULL ? '=': '!') << "= NULL";
    #endif
  
    while(next != NULL && next->free) {
        #ifdef MY_DEBUG
            std::cout << "; free = " << (unsigned)(next->free) << "\n";
        #endif
        joinBlocks(ph, next);

        ph = next;
        next = iterator(ph);

        #ifdef MY_DEBUG
            std::cout << "   next header = " << next << " " << (next == NULL ? '=': '!') << "= NULL";
        #endif
    }
    #ifdef MY_DEBUG
        std::cout << "\n";
    #endif
  
    return ph;
}

void mysetup(void *buf, std::size_t size) {
    _size = size;
    _startHeader = initBlock(buf, size);
    _endTail = getTail(_startHeader);
    minSize = getAllSize(MIN_BLOCK_SIZE);
}

void* myalloc(std::size_t size) {
    #ifdef MY_DEBUG
        std::cout << "\n--------------------START MEMORY ALLOC--------------------\n\n";
    #endif

    std::size_t allSize = getAllSize(size);
    struct header* ph = _startHeader;
    while(ph != NULL) {
        if(ph->free && (ph->actualSize >= allSize || ph->actualSize >= size)) {
            unsigned char *base = (unsigned char*)utilizeBlock(ph, size);
            base += sizeof(struct header);

            #ifdef MY_DEBUG
                std::cout << "\n--------------------MEMORY ALLOC--------------------\n\n";
            #endif

            return (void*)base;
        }
      ph = getNext(ph);
    }
    return NULL;
}

void myfree(void *p) {
    #ifdef MY_DEBUG
        std::cout << "\n--------------------START MEMORY FREE--------------------\n\n";
    #endif

    #ifdef MY_DEBUG
        std::cout << "myfree(" << p << ')' << "\n";
    #endif
    
    if(p == NULL)
        return;
    
    struct header* ph = (struct header*)((unsigned char*)p - sizeof(struct header));
    ph->free = 1;
    
    #ifdef MY_DEBUG
        std::cout << "   ph addr = " << ph << "; free = " << (unsigned int)(ph->free) << "; actualSize = " << ph->actualSize << "\n";
    #endif
  
    struct tail* pt = getTail(ph);
    pt->free = 1;
    
    #ifdef MY_DEBUG
        std::cout << "   pt addr = " << pt << "; free = " << (unsigned int)(pt->free) << "; actualSize = " << pt->actualSize << "\n";
    #endif
  
    ph = joinNearestFreeBlocks(ph, getPrevious);
    joinNearestFreeBlocks(ph, getNext);

    #ifdef MY_DEBUG
        std::cout << "\n--------------------MEMORY FREE--------------------\n\n";
    #endif
}

