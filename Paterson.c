#include <stdio.h>
#include <stdatomic.h>

void critical_section(){
    printf("AAAAAAAAAAAAAAAAAAAA\n");
}

struct lock {
    atomic_int flag[2]; 
};

void lock_init(struct lock *lock){
    atomic_store(&lock->flag[0], 0);
    atomic_store(&lock->flag[1], 0);
}

void lock(struct lock *lock, int threadId){
    const int me = threadId;
    const int other = 1 - me;

    atomic_store(&lock->flag[me], 1);
    while (atomic_load(&lock->flag[other])){}
}

void unlock(struct lock *lock, int threadId){
    const int me = threadId;
    atomic_store(&lock->flag[me], 0);
}

int main(void){
    struct lock my_lock;
    lock_init(&my_lock);

    int threadId = 0;

    lock(&my_lock, threadId);
    critical_section();
    unlock(&my_lock, threadId);
    
    return 0;
}


