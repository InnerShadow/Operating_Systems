#include <iostream>
#include <climits>
#include <vector>

#define MIN_SEGMENT_SIZE 4

typedef struct header* (*HeaderIterator)(struct header*);

struct header {
    unsigned char free;
    std::size_t actualSize;
};

struct tail {
    unsigned char free; 
    std::size_t actualSize;
};

struct segment {
    struct header* startHeader;
    struct tail* endTail;
    std::size_t size;
    std::size_t minSize;
};

struct header* StartHeader;
struct tail* EndTail;
std::size_t Size;
std::size_t minSize;

struct header* InitSegment(void* buf, std::size_t size) {
    struct header* ph = (struct header*)buf;
    ph->free = 1;
    ph->actualSize = size - sizeof(struct header) - sizeof(struct tail);

    struct tail* pt = (struct tail*)((unsigned char*)buf + size - sizeof(struct tail));
    pt->free = 1;
    pt->actualSize = ph->actualSize;

    return ph;
}

struct tail* GetTail(struct header* phead) {
    unsigned char* base = (unsigned char*)phead;
    base += sizeof(struct header) + phead->actualSize;
    return (struct tail*)base;
}

struct header* GetHeader(struct tail* ptail) {
    unsigned char* base = (unsigned char*)ptail;
    base -= sizeof(struct header) + ptail->actualSize;
    return (struct header*)base;
}

struct header* GetNext(struct header* phead) {
    std::cout << "getNext(" << phead << ')' << "\n";

    struct tail* ptail = GetTail(phead);
    if (ptail == EndTail) {
        return nullptr;
    }

    unsigned char* base = (unsigned char*)ptail;
    base += sizeof(struct tail);
    return (struct header*)base;
}

struct header* GetPrev(struct header* phead) {
    if (phead == StartHeader) {
        return nullptr;
    }

    unsigned char* base = (unsigned char*)phead;
    base -= sizeof(struct tail);
    return GetHeader((struct tail*)base);
}

std::size_t AllSize(std::size_t size) {
    return size + sizeof(struct header) + sizeof(struct tail);
}

std::size_t ActualSize(struct header* ph, struct tail* pt) {
    if ((unsigned char*)ph >= (unsigned char*)pt) {
        return 0;
    }

    return (unsigned char*)pt - (unsigned char*)ph - sizeof(struct header);
}

struct header* InitSegmentData(struct header* ph, std::size_t size) {
    if(!ph->free || ph->actualSize < size){
        return NULL;
    }
  
    std::size_t allSize = AllSize(size);
  
    struct tail* ptEnd = GetTail(ph);
  
    if(ph->actualSize <= allSize) {
        ph->free = 0;
        ptEnd->free = 0;
        return ph;
    }
  
    std::size_t newSize = ph->actualSize - allSize;
  
    if(newSize < MIN_SEGMENT_SIZE) {
        ph->free = 0;
        ptEnd->free = 0;
        return ph;
    }
  
    unsigned char* base = (unsigned char*)ptEnd;
    base -= allSize;
  
    struct tail* ptStart = (struct tail*)base;
    ptStart->free = 1;
    ptStart->actualSize = ActualSize(ph, ptStart);
    ph->actualSize = ptStart->actualSize;
  
    struct header* phEnd = (struct header*)(base + sizeof(struct tail));
    phEnd->free = 0;
    phEnd->actualSize = ActualSize(phEnd, ptEnd);
    
    ptEnd->free = 0;
    ptEnd->actualSize = phEnd->actualSize;
    
    return phEnd;
}

void Setup(void* buf, std::size_t size) {
    Size = size;
    StartHeader = InitSegment(buf, size);
    EndTail = GetTail(StartHeader);
    minSize = AllSize(MIN_SEGMENT_SIZE);
}

void* Alloc(std::size_t size) {
    std::size_t allSize = AllSize(size);
    struct header* ph = StartHeader;
    while (ph != nullptr && (ph->actualSize >= allSize || ph->actualSize >= size)) {
        if (ph->free) {
            unsigned char* base = (unsigned char*)InitSegmentData(ph, size);
            base += sizeof(struct header);

            std::cout << "Allocated " << size << " bytes at address: " << (void*)base << std::endl;

            return (void*)base;
        }
        ph = GetNext(ph);

        std::cout << "\nNEXT: " << ph << "\n\n";
    }

    std::cout << "Allocation failed for " << size << " bytes." << std::endl;

    return nullptr;
}

void JoinHelper(struct header* phStart, struct header* phEnd) {
    if (phStart == phEnd) {
        return;
    }

    if (phStart > phEnd) {
        struct header* tmp = phStart;
        phStart = phEnd;
        phEnd = tmp;
    }

    struct tail* ptEnd = GetTail(phEnd);

    phStart->actualSize = ActualSize(phStart, ptEnd);
    ptEnd->actualSize = phStart->actualSize;
}

struct header* JoinFreeBlocks(struct header* ph, HeaderIterator iterator) {
    if(ph == NULL || !ph->free)
        return ph;
  
    struct header* next = iterator(ph);
 
    while(next != NULL && next->free) {
        JoinHelper(ph, next);

        ph = next;
        next = iterator(ph);
    }
  
    return ph;
}

void Free(void* p) {
    if (p == nullptr){
        return;
    }

    struct header* ph = (struct header*)((unsigned char*)p - sizeof(struct header));
    ph->free = 1;

    struct tail* pt = GetTail(ph);
    pt->free = 1;

    std::cout << "Freed memory at address: " << p << std::endl;

    ph = JoinFreeBlocks(ph, GetNext);
    JoinFreeBlocks(ph, GetPrev);

    StartHeader = ph;
}

template <typename T>
struct Allocator {
    using value_type = T;

    Allocator(struct segment* seg) noexcept : segment(seg) {}

    template <typename U>
    Allocator(const Allocator<U>& other) noexcept : segment(other.segment) {}

    T* allocate(std::size_t n) {
        return static_cast<T*>(Alloc(n * sizeof(T)));
    }

    void deallocate(T* p, std::size_t n) {
        Free(p);
    }

    struct segment* segment;
};

int main(void) {
    const std::size_t memSize = 1024 * 4;

    void* memoryBlock = std::malloc(memSize);

    struct segment mySegment;
    Setup(memoryBlock, memSize);
    mySegment.startHeader = StartHeader;
    mySegment.endTail = EndTail;
    mySegment.size = Size;
    mySegment.minSize = minSize;

    Allocator<int> allocator(&mySegment);

    auto p1 = allocator.allocate(1);
    auto p2 = allocator.allocate(1);
    auto p3 = allocator.allocate(1);

    *p1 = -1;
    *p2 = -2;
    *p3 = -3;

    std::cout << "\n\nAfter allegate)\n";
    std::cout << "P1 = " << p1 << " p2 = " << p2 << " p3 = " << p3 << "\n";
    std::cout << "value p1 = " << *p1 << " p2 = " << *p2 << " p3 = " << *p3 << "\n";

    allocator.deallocate(p2, 1);

    std::cout << "\n\nAfter deallegate)\n";
    std::cout << "P1 = " << p1 << " p2 = " << p2 << " p3 = " << p3 << "\n";
    std::cout << "value p1 = " << *p1 << " p2 = " << *p2 << " p3 = " << *p3 << "\n";

    auto p4 = allocator.allocate(1);

    *p4 = -4;

    std::cout << "\n\nAfter additional allegate)\n";
    std::cout << "P1 = " << p1 << " p2 = " << p2 << " p3 = " << p3 << " p4 = " << p4 << "\n";
    std::cout << "value p1 = " << *p1 << " p2 = " << *p2 << " p3 = " << *p3 << " p4 = " << *p4 <<"\n";

    std::free(memoryBlock);
    return 0;
}

