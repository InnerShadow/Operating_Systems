#include <stdio.h>
#include <stdatomic.h>

void critical_section(){
    printf("AAAAAAAAAAAAAAAAAAAA\n");
}

struct lock {
    atomic_int flag[2];
    atomic_int turn;
};

void acquire_lock(struct lock *my_lock, int thread_id){
    atomic_store(&my_lock->flag[thread_id], 1);
    atomic_store(&my_lock->turn, thread_id);

    while (atomic_load(&my_lock->flag[1 - thread_id]) && atomic_load(&my_lock->turn) == thread_id) {}
}

void release_lock(struct lock *my_lock, int thread_id){
    atomic_store(&my_lock->flag[thread_id], 0);
}

int main(void){
    struct lock my_lock = { {0, 0}, -1 };

    int thread_id = 0;

    acquire_lock(&my_lock, thread_id);
    critical_section();
    release_lock(&my_lock, thread_id);

    return 0;
}

