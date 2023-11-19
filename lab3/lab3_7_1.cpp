#include <iostream>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>
#include <boost/interprocess/sync/interprocess_semaphore.hpp>

const int MAX_CHAIRS = 5;

std::mutex barberMutex;
std::mutex waitingRoomMutex;
std::queue<int> waitingCustomers;
std::condition_variable cv;
boost::interprocess::interprocess_semaphore availableChairs(MAX_CHAIRS);

void barberCutHair(int customerId) {
    std::lock_guard<std::mutex> lock(barberMutex);
    std::cout << "Barber is cutting hair for Customer " << customerId << "." << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(4));
}

void barberSleep() {
    std::unique_lock<std::mutex> lock(barberMutex);
    std::cout << "Barber is sleeping." << std::endl;
    cv.wait(lock);
}

void customerArrives(int customerId) {
    if (availableChairs.try_wait()) {
        std::lock_guard<std::mutex> lock(waitingRoomMutex);
        waitingCustomers.push(customerId);
        std::cout << "Customer " << customerId << " takes a seat in the reception area." << std::endl;
        cv.notify_one();
    } else {
        std::cout << "------ No available chairs. Customer " << customerId << " leaves ------" << std::endl;
    }
}

void simulateCustomerArrivals() {
    for (std::size_t i = 0; i <= 20; ++i) {
        customerArrives(i);

        std::this_thread::sleep_for(std::chrono::seconds(1));
        if (i > 7) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
        }
        if(i == 15){
            std::this_thread::sleep_for(std::chrono::seconds(25));
        }
    }
}

void barberShopSimulation() {
    while (true) {
        std::unique_lock<std::mutex> lock(waitingRoomMutex);

        if (!waitingCustomers.empty()) {
            int id = waitingCustomers.front();
            waitingCustomers.pop();
            lock.unlock();

            std::cout << "Customer " << id << " is invited by the barber." << std::endl;
            barberCutHair(id);
            availableChairs.post();
        } else {
            lock.unlock();
            barberSleep();
        }
    }
}

int main(void) {
    std::thread arrivalThread(simulateCustomerArrivals);
    std::thread simulationThread(barberShopSimulation);

    arrivalThread.join();
    simulationThread.join();

    return 0;
}
