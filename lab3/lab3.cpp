#include <iostream>
#include <thread>
#include <mutex>

static const int NUM_PHILOSOPHERS = 5;
std::mutex forks[NUM_PHILOSOPHERS];

void philosopher(int id) {
    int left_fork = (id + 1) % NUM_PHILOSOPHERS;
    int right_fork = id;

    while (true) {
        std::cout << "The Philosopher " << id << " is thinking..." << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(1000 * (id + 1)));

        std::unique_lock<std::mutex> left(forks[left_fork]);
        std::cout << "The Philosopher " << id << " take left fork." << std::endl;

        std::this_thread::sleep_for(std::chrono::milliseconds(1000 * (id + 1)));

        std::unique_lock<std::mutex> right(forks[right_fork]);
        std::cout << "The Philosopher " << id << " take right fork and start eating." << std::endl;

        std::this_thread::sleep_for(std::chrono::milliseconds(1000 * (id + 1)));

        std::cout << "The Philosopher " << id << " leave forks and start thinking..." << std::endl;
    }
}

int main(void) {
    std::thread philosophers[NUM_PHILOSOPHERS];
    
    for (size_t i = 0; i < NUM_PHILOSOPHERS; ++i) {
        philosophers[i] = std::thread(philosopher, i);
    }
    
    for (size_t i = 0; i < NUM_PHILOSOPHERS; ++i) {
        philosophers[i].join();
    }

    return 0;
}
