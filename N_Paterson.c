#include <stdio.h>
#include <stdlib.h>
#include <stdatomic.h>

const int N = 10;

void critical_section() {
    printf("AAAAAAAAAAAAAAAAAAAA\n");
}

struct lock_one {
    atomic_int last;
    atomic_int flag;
};

struct lock {
    struct lock_one *locks;
};

int flags_clear(const struct lock_one *lock, int me) {
    for (size_t i = 0; i < N; ++i) {
        if (i != me && atomic_load(&lock[i].flag)) {
            return 0;
        }
    }
    return 1;
}

void lock_one(struct lock_one *lock, int me) {
    atomic_store(&lock->flag, 1);
    atomic_store(&lock->last, me);

    while ((!flags_clear(lock, me)) && atomic_load(&lock->last) == me) {}
}

void unlock_one(struct lock_one *lock) {
    atomic_store(&lock->flag, 0);
}

void lock(struct lock *locks, int me) {
    for (size_t i = 0; i < N; ++i) {
        lock_one(&locks->locks[i], me);
    }
}

void unlock(struct lock *locks) {
    for (int i = N - 1; i >= 0; --i) {
        unlock_one(&locks->locks[i]);
    }
}

int main(void) {
    struct lock locks;
    locks.locks = (struct lock_one *)malloc(N * sizeof(struct lock_one));
    
    for (int thread_id = 0; thread_id < N; ++thread_id) {
        lock(&locks, thread_id);
        critical_section();
        unlock(&locks);
    }

    free(locks.locks); 

    return 0;
}

