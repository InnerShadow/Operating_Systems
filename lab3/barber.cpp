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
        perror("shm_open");
        exit(EXIT_FAILURE);
    }

    if (ftruncate(shm_fd, sizeof(sem_t)) == -1) {
        perror("ftruncate");
        exit(EXIT_FAILURE);
    }

    sem_t *semaphore = (sem_t *)mmap(NULL, sizeof(sem_t), PROT_READ | PROT_WRITE, MAP_SHARED, shm_fd, 0);
    if (semaphore == MAP_FAILED) {
        perror("mmap");
        exit(EXIT_FAILURE);
    }

    sem_init(semaphore, 1, MAX_CHAIRS);

    int customerId;
    while (true) {

        int value = 0;
        sem_getvalue(semaphore, &value);
        //std::cout << value << "\n";

        if (value == MAX_CHAIRS){
            std::cout << "Barber is going to sleep." << "\n"; 
        }

        read(fifo, &customerId, sizeof(customerId));
        std::cout << "Barber is cutting hair for Customer " << customerId << "." << std::endl;
        sleep(4);

        sem_post(semaphore);
    }

    close(fifo);
    sem_destroy(semaphore);
    munmap(semaphore, sizeof(sem_t));
    shm_unlink(semaphoreName);

    return 0;
}
