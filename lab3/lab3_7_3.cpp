#include <iostream>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <semaphore.h>
#include <sys/mman.h>
#include <sys/wait.h>

const char *fifoPath = "/tmp/barber_fifoP2";
const char *semaphoreName = "/barber_semaphoreP2";

const int MAX_CHAIRS = 5;

void barberShopSimulation() {
    int fifo = open(fifoPath, O_RDONLY);

    int shm_fd = shm_open(semaphoreName, O_CREAT | O_RDWR, 0666);
    if (shm_fd == -1) {
        std::cerr << "shm_open error" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (ftruncate(shm_fd, sizeof(sem_t)) == -1) {
        std::cerr << "ftruncate error" << std::endl;
        exit(EXIT_FAILURE);
    }

    sem_t *semaphore = (sem_t *)mmap(NULL, sizeof(sem_t), PROT_READ | PROT_WRITE, MAP_SHARED, shm_fd, 0);
    if (semaphore == MAP_FAILED) {
        std::cerr << "mmap error" << std::endl;
        exit(EXIT_FAILURE);
    }

    sem_init(semaphore, 1, MAX_CHAIRS);

    int customerId;
    while (true) {
        int value = 0;
        sem_getvalue(semaphore, &value);

        if (value == MAX_CHAIRS) {
            std::cout << "Barber is going to sleep." << std::endl;
            sleep(2);
            continue;
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
}

void simulateCustomerArrivals() {
    int fifo = open(fifoPath, O_WRONLY);

    int shm_fd = shm_open(semaphoreName, O_CREAT | O_RDWR, 0666);
    if (shm_fd == -1) {
        std::cerr << "shm_open error" << std::endl;
        exit(EXIT_FAILURE);
    }

    sem_t *semaphore = (sem_t *)mmap(NULL, sizeof(sem_t), PROT_READ | PROT_WRITE, MAP_SHARED, shm_fd, 0);
    if (semaphore == MAP_FAILED) {
        std::cerr << "mmap error" << std::endl;
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
}

int main() {
    mkfifo(fifoPath, 0666);

    pid_t barber_pid = fork();

    if (barber_pid == 0) {
        barberShopSimulation();
    } else {
        simulateCustomerArrivals();
    }

    return 0;
}
