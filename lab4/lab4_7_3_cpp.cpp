#include <iostream>
#include <memory>

constexpr std::size_t MIN_SEGMENT_SIZE = 4;

class Header {
public:
    Header() = default;

    unsigned char& GetFree() { return free; }  // Изменено здесь
    std::size_t& GetActualSize() { return actualSize; }  // Изменено здесь

private:
    unsigned char free;
    std::size_t actualSize;
};

class Tail {
public:
    Tail() = default;

    unsigned char& GetFree() { return free; }  // Изменено здесь
    std::size_t& GetActualSize() { return actualSize; }  // Изменено здесь

private:
    unsigned char free;
    std::size_t actualSize;
};

class Segment {
public:
    Header* GetStartHeader() const { return startHeader; }
    Tail* GetEndTail() const { return endTail; }
    std::size_t GetSize() const { return size; }
    std::size_t GetMinSize() const { return minSize; }

    void SetStartHeader(Header* header) { startHeader = header; }
    void SetEndTail(Tail* tail) { endTail = tail; }
    void SetSize(std::size_t s) { size = s; }
    void SetMinSize(std::size_t min) { minSize = min; }

private:
    Header* startHeader;
    Tail* endTail;
    std::size_t size;
    std::size_t minSize;
};


Header* StartHeader;
Tail* EndTail;
std::size_t Size;
std::size_t minSize;

Header* InitSegment(void* buf, std::size_t size) {
    auto ph = new (buf) Header();
    ph->GetFree() = 1;
    ph->GetActualSize() = size - sizeof(Header) - sizeof(Tail);

    auto pt = new (reinterpret_cast<unsigned char*>(buf) + size - sizeof(Tail)) Tail();
    pt->GetFree() = 1;
    pt->GetActualSize() = ph->GetActualSize();

    return ph;
}

Tail* GetTail(Header* phead) {
    unsigned char* base = reinterpret_cast<unsigned char*>(phead);
    base += sizeof(Header) + phead->GetActualSize();
    return reinterpret_cast<Tail*>(base);
}

Header* GetHeader(Tail* ptail) {
    unsigned char* base = reinterpret_cast<unsigned char*>(ptail);
    base -= sizeof(Header) + ptail->GetActualSize();
    return reinterpret_cast<Header*>(base);
}

Header* GetNext(Header* phead) {
    std::cout << "getNext(" << phead << ')' << "\n";

    Tail* ptail = GetTail(phead);
    if (ptail == EndTail) {
        return nullptr;
    }

    unsigned char* base = reinterpret_cast<unsigned char*>(ptail);
    base += sizeof(Tail);
    return reinterpret_cast<Header*>(base);
}

Header* GetPrev(Header* phead) {
    if (phead == StartHeader) {
        return nullptr;
    }

    unsigned char* base = reinterpret_cast<unsigned char*>(phead);
    base -= sizeof(Tail);
    return GetHeader(reinterpret_cast<Tail*>(base));
}


std::size_t AllSize(std::size_t size) {
    return size + sizeof(Header) + sizeof(Tail);
}

std::size_t ActualSize(Header* ph, Tail* pt) {
    if (reinterpret_cast<unsigned char*>(ph) >= reinterpret_cast<unsigned char*>(pt)) {
        return 0;
    }

    return reinterpret_cast<unsigned char*>(pt) - reinterpret_cast<unsigned char*>(ph) - sizeof(Header);
}

Header* InitSegmentData(Header* ph, std::size_t size) {
    if (!ph->GetFree() || ph->GetActualSize() < size) {
        return nullptr;
    }

    std::size_t allSize = AllSize(size);

    Tail* ptEnd = GetTail(ph);

    if (ph->GetActualSize() <= allSize) {
        ph->GetFree() = 0;
        ptEnd->GetFree() = 0;
        return ph;
    }

    std::size_t newSize = ph->GetActualSize() - allSize;

    if (newSize < MIN_SEGMENT_SIZE) {
        ph->GetFree() = 0;
        ptEnd->GetFree() = 0;
        return ph;
    }

    unsigned char* base = reinterpret_cast<unsigned char*>(ptEnd);
    base -= allSize;

    Tail* ptStart = new (base) Tail();
    ptStart->GetFree() = 1;
    ptStart->GetActualSize() = ActualSize(ph, ptStart);
    ph->GetActualSize() = ptStart->GetActualSize();

    Header* phEnd = new (base + sizeof(Tail)) Header();
    phEnd->GetFree() = 0;
    phEnd->GetActualSize() = ActualSize(phEnd, ptEnd);

    ptEnd->GetFree() = 0;
    ptEnd->GetActualSize() = phEnd->GetActualSize();

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
    Header* ph = StartHeader;
    while (ph != nullptr && (ph->GetActualSize() >= allSize || ph->GetActualSize() >= size)) {
        if (ph->GetFree()) {
            unsigned char* base = reinterpret_cast<unsigned char*>(InitSegmentData(ph, size));
            base += sizeof(Header);

            std::cout << "Allocated " << size << " bytes at address: " << static_cast<void*>(base) << std::endl;

            return static_cast<void*>(base);
        }
        ph = GetNext(ph);

        std::cout << "\nNEXT: " << ph << "\n\n";
    }

    std::cout << "Allocation failed for " << size << " bytes." << std::endl;

    return nullptr;
}

void JoinHelper(Header* phStart, Header* phEnd) {
    if (phStart == phEnd) {
        return;
    }

    if (phStart > phEnd) {
        Header* tmp = phStart;
        phStart = phEnd;
        phEnd = tmp;
    }

    Tail* ptEnd = GetTail(phEnd);

    phStart->GetActualSize() = ActualSize(phStart, ptEnd);
    ptEnd->GetActualSize() = phStart->GetActualSize();
}

Header* JoinFreeBlocks(Header* ph, Header* (*iterator)(Header*)) {
    if (ph == nullptr || !ph->GetFree())
        return ph;

    Header* next = iterator(ph);

    while (next != nullptr && next->GetFree()) {
        JoinHelper(ph, next);

        ph = next;
        next = iterator(ph);
    }

    return ph;
}

void Free(void* p) {
    if (p == nullptr) {
        return;
    }

    Header* ph = reinterpret_cast<Header*>(static_cast<unsigned char*>(p) - sizeof(Header));
    ph->GetFree() = 1;

    Tail* pt = GetTail(ph);
    pt->GetFree() = 1;

    std::cout << "Freed memory at address: " << p << std::endl;

    ph = JoinFreeBlocks(ph, GetNext);
    JoinFreeBlocks(ph, GetPrev);

    StartHeader = ph;
}

template <typename T>
class Allocator {
public:
    using value_type = T;

    Allocator(Segment* seg) noexcept : segment(seg) {}

    template <typename U>
    Allocator(const Allocator<U>& other) noexcept : segment(other.segment) {}

    T* allocate(std::size_t n) {
        return static_cast<T*>(Alloc(n * sizeof(T)));
    }

    void deallocate(T* p, std::size_t n) {
        Free(p);
    }

private:
    Segment* segment;
};

int main(void) {
    const std::size_t memSize = 1024 * 4;

    void* memoryBlock = std::malloc(memSize);

    Segment mySegment;
    Setup(memoryBlock, memSize);
    mySegment.SetStartHeader(StartHeader);
    mySegment.SetEndTail(EndTail);
    mySegment.SetSize(Size);
    mySegment.SetMinSize(minSize);

    Allocator<int> allocator(&mySegment);

    auto p1 = allocator.allocate(1);
    auto p2 = allocator.allocate(1);
    auto p3 = allocator.allocate(1);

    *p1 = -1;
    *p2 = -2;
    *p3 = -3;

    std::cout << "\n\nAfter allocate)\n";
    std::cout << "P1 = " << p1 << " p2 = " << p2 << " p3 = " << p3 << "\n";
    std::cout << "value p1 = " << *p1 << " p2 = " << *p2 << " p3 = " << *p3 << "\n";

    allocator.deallocate(p2, 1);

    std::cout << "\n\nAfter deallocate)\n";
    std::cout << "P1 = " << p1 << " p2 = " << p2 << " p3 = " << p3 << "\n";
    std::cout << "value p1 = " << *p1 << " p2 = " << *p2 << " p3 = " << *p3 << "\n";

    auto p4 = allocator.allocate(1);

    *p4 = -4;

    std::cout << "\n\nAfter additional allocate)\n";
    std::cout << "P1 = " << p1 << " p2 = " << p2 << " p3 = " << p3 << " p4 = " << p4 << "\n";
    std::cout << "value p1 = " << *p1 << " p2 = " << *p2 << " p3 = " << *p3 << " p4 = " << *p4 << "\n";

    std::free(memoryBlock);
    return 0;
}
