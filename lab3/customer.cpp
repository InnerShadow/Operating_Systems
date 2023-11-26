#include <iostream>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <semaphore.h>
#include <sys/mman.h>

const char *fifoPath = "/tmp/barber_fifoP";
const char *semaphoreName = "/barber_semaphore";

int main() {
    int fifo = open(fifoPath, O_WRONLY);

    int shm_fd = shm_open(semaphoreName, O_CREAT | O_RDWR, 0666);
    if (shm_fd == -1) {
        perror("shm_open");
        exit(EXIT_FAILURE);
    }

    sem_t *semaphore = (sem_t *)mmap(NULL, sizeof(sem_t), PROT_READ | PROT_WRITE, MAP_SHARED, shm_fd, 0);
    if (semaphore == MAP_FAILED) {
        perror("mmap");
        exit(EXIT_FAILURE);
    }

    for (std::size_t i = 0; i <= 20; ++i) {
        if (sem_trywait(semaphore) == 0) {
            int customerId = i;
            write(fifo, &customerId, sizeof(customerId));
            std::cout << "Customer " << customerId << " takes a seat in the reception area." << std::endl;
        } else {
            std::cout << "------ No available chairs. Customer " << i << " leaves ------" << std::endl;
        }

        sleep(1);
        if (i > 7) {
            sleep(1);
        }
        if (i == 15) {
            sleep(20);
        }        
    }

    close(fifo);
    sem_destroy(semaphore);
    munmap(semaphore, sizeof(sem_t));
    return 0;
}
