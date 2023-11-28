#include <iostream>
#include <unistd.h>
#include <fcntl.h>
#include <fstream>
#include <semaphore.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <string>
#include <vector>

const char* semaphoreName = "/barber_semaphoreP28";
const char* mutexName = "/queue_mutexP28";
const char* file_path = "queue.txt";

const int MAX_CHAIRS = 5;

void barberShopSimulation() {
    sem_t* semaphore = sem_open(semaphoreName, O_CREAT, 0666, MAX_CHAIRS);
    sem_t* mutex = sem_open(mutexName, O_CREAT, 0666, 1);

    if (semaphore == SEM_FAILED || mutex == SEM_FAILED) {
        std::cout << "sem_open error" << std::endl;
        exit(EXIT_FAILURE);
    }

    int sleep_num = 0;

    while (true) {
        int value = 0;

        sem_getvalue(semaphore, &value);

        if (value == MAX_CHAIRS) {
            std::cout << "Barber is going to sleep. Till end " << 2 * (10 - sleep_num) << " seconds" << std::endl;
            sleep(2);
            sleep_num++;
            if (sleep_num > 10) {
                std::cout << "FINISHED\n";
                break;
            }
            continue;
        }

        sleep_num = 0;

        sem_wait(mutex);

        std::ifstream ifF(file_path);
        if (!ifF.is_open()){
            std::cout << "Cannot open file. Finish!\n";
            break; 
        }

        std::string id = "";

        std::vector<int> queue;

        while (std::getline(ifF, id)) {
            queue.push_back(std::stoi(id));
        }

        ifF.close();

        std::ofstream ofF(file_path);
        if (!ofF.is_open()) {
            std::cout << "Cannot open file. Finish!\n";
            break; 
        }

        for (std::size_t i = 1; i < queue.size(); ++i) {
            ofF << std::to_string(queue.at(i)) << "\n";
        }

        ofF.close();

        sem_post(mutex);

        std::cout << "Barber is cutting hair for Customer " << queue.at(0) << "." << std::endl;
        sleep(4);

        sem_post(semaphore);
    }

    sem_close(semaphore);
    sem_unlink(semaphoreName);
    sem_close(mutex);
    sem_unlink(mutexName);
}

void simulateCustomerArrivals() {
    sem_t* semaphore = sem_open(semaphoreName, O_CREAT, 0666, MAX_CHAIRS);
    sem_t* mutex = sem_open(mutexName, O_CREAT, 0666, 1);

    if (semaphore == SEM_FAILED || mutex == SEM_FAILED) {
        std::cout << "sem_open error" << std::endl;
        exit(EXIT_FAILURE);
    }

    pid_t main_pid = getpid();

    for (std::size_t i = 0; i <= 20; ++i) {
        if (getpid() == main_pid) {

            pid_t child_pid = fork();

            if (child_pid == 0) {
                if (sem_trywait(semaphore) == 0) {
                    sem_wait(mutex);

                    int customerId = i;
                    std::cout << "P_id = " << getpid() << " Customer " << customerId << " takes a seat in the reception area." << std::endl;

                    std::ofstream f(file_path, std::ios::app);
                    if (!f.is_open()) {
                        std::cout << "P_id = " << getpid() << " ------ No available chairs. Customer " << i << " leaves ------" << std::endl;
                        continue;
                    }

                    f << std::to_string(i) << "\n";

                    f.close();

                    sem_post(mutex);
                } else {
                    std::cout << "P_id = " << getpid() << " ------ No available chairs. Customer " << i << " leaves ------" << std::endl;
                }
            }
        }
        sleep(1);
        if (i > 7) {
            sleep(1);
        }
        if (i == 15) {
            sleep(20);
        }
    }

    sem_close(semaphore);
    sem_close(mutex);
}

int main() {
    pid_t barber_pid = fork();

    if (barber_pid == 0) {
        sleep(1);
        simulateCustomerArrivals();
    } else {
        barberShopSimulation();
    }

    return 0;
}


