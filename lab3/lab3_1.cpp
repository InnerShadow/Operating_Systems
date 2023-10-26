#include <iostream>
#include <pthread.h>
#include <unistd.h>

static const int NUM_PHILOSOPHERS = 5;
pthread_mutex_t forks[NUM_PHILOSOPHERS];

void philosopher(int id) {
    int left_fork = (id + 1) % NUM_PHILOSOPHERS;
    int right_fork = id;

    while (true){
        std::cout << "The Philosopher " << id << " is thinking..." << std::endl;
        usleep(1000000 * (id + 1));

        pthread_mutex_lock(&forks[left_fork]);
        std::cout << "The Philosopher " << id << " take left fork." << std::endl;

        usleep(1000000 * (id + 1));

        pthread_mutex_lock(&forks[right_fork]);
        std::cout << "The Philosopher " << id << " take right fork and start eating." << std::endl;

        usleep(1000000 * (id + 1));

        std::cout << "The Philosopher " << id << " leave forks and start thinking..." << std::endl;

        pthread_mutex_unlock(&forks[left_fork]);
        pthread_mutex_unlock(&forks[right_fork]);
    }
}

int main(void) {
    pthread_t philosophers[NUM_PHILOSOPHERS];

    for (size_t i = 0; i < NUM_PHILOSOPHERS; ++i){
        pthread_create(&philosophers[i], NULL, reinterpret_cast<void*(*)(void*)>(philosopher), reinterpret_cast<void*>(i));
    }

    for (size_t i = 0; i < NUM_PHILOSOPHERS; ++i){
        pthread_join(philosophers[i], NULL);
    }

    return 0;
}

