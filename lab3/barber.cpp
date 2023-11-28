#include <iostream>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <semaphore.h>
#include <sys/mman.h>

const char *fifoPath = "/tmp/barber_fifoP";
const char *semaphoreName = "/barber_semaphore";

const int MAX_CHAIRS = 5;

int main() {
    mkfifo(fifoPath, 0666);

    int fifo = open(fifoPath, O_RDONLY);

    int shm_fd = shm_open(semaphoreName, O_CREAT | O_RDWR, 0666);
    if (shm_fd == -1) {
        std::cout << "shm_open error" << "\n";
        return 1;
    }

    if (ftruncate(shm_fd, sizeof(sem_t)) == -1) {
        std::cout << "ftruncate error" << "\n";
        return 2;
    }

    sem_t *semaphore = (sem_t *)mmap(NULL, sizeof(sem_t), PROT_READ | PROT_WRITE, MAP_SHARED, shm_fd, 0);
    if (semaphore == MAP_FAILED) {
        std::cout << "mmap error" << "\n";
        return 3;
    }

    sem_init(semaphore, 1, MAX_CHAIRS);

    int customerId;
    while (true) {
        int value = 0;
        sem_getvalue(semaphore, &value);
        //std::cout << " sem: " << value << "\n";

        if (value == MAX_CHAIRS){
            std::cout << "Barber is going to sleep." << "\n";
            sleep(2);
            continue;
        }

        read(fifo, &customerId, sizeof(customerId));
        std::cout << "Barber is cutting hair for Customer " << customerId << "." << "\n";
        sleep(4);

        sem_post(semaphore);
    }

    close(fifo);
    sem_destroy(semaphore);
    munmap(semaphore, sizeof(sem_t));
    shm_unlink(semaphoreName);

    return 0;
}

