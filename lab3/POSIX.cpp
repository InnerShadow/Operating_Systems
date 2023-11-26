#include <iostream>
#include <queue>
#include <pthread.h>
#include <semaphore.h>
#include <unistd.h>

const int MAX_CHAIRS = 5;

pthread_mutex_t barberMutex = PTHREAD_MUTEX_INITIALIZER;
std::queue<int> waitingCustomers;
sem_t availableChairs;
bool barber = false;

void* barberCutHair(void* arg) {
    int customerId = *((int*)arg);
    pthread_mutex_lock(&barberMutex);
    std::cout << "Barber is cutting hair for Customer " << customerId << "." << std::endl;
    sleep(4);
    pthread_mutex_unlock(&barberMutex);
    return nullptr;
}

void customerArrives(int customerId) {
    if (sem_trywait(&availableChairs) == 0) {
        waitingCustomers.push(customerId);
        std::cout << "Customer " << customerId << " takes a seat in the reception area." << std::endl;
    } else {
        std::cout << "------ No available chairs. Customer " << customerId << " leaves ------" << std::endl;
    }
}

void* simulateCustomerArrivals(void* arg) {
    for (std::size_t i = 0; i <= 20; ++i) {
        customerArrives(i);

        sleep(1);
        if (i > 7) {
            sleep(1);
        }
        if (i == 15) {
            sleep(20);
        }
    }
    return nullptr;
}

void* barberShopSimulation(void* arg) {
    while (true) {
        while (waitingCustomers.empty()) {
            std::cout << "Barber is going to sleep." << std::endl;
            sleep(2);
        }

        int id = waitingCustomers.front();
        waitingCustomers.pop();

        std::cout << "Customer " << id << " is invited by the barber." << std::endl;

        pthread_t barberThread;
        pthread_create(&barberThread, nullptr, barberCutHair, &id);
        pthread_join(barberThread, nullptr);

        sem_post(&availableChairs);
    }
    return nullptr;
}

int main(void) {
    sem_init(&availableChairs, 0, MAX_CHAIRS);

    if(barber){
        pthread_t simulationThread;
        pthread_create(&simulationThread, nullptr, barberShopSimulation, nullptr);

        pthread_join(simulationThread, nullptr);
    } else {
        pthread_t arrivalThread;
        pthread_create(&arrivalThread, nullptr, simulateCustomerArrivals, nullptr);
        pthread_join(arrivalThread, nullptr);
    }    

    sem_destroy(&availableChairs);

    return 0;
}