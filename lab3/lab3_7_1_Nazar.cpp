#include <iostream>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>
#include <boost/interprocess/sync/interprocess_semaphore.hpp>

using namespace std;

const int STULIJA = 5;

mutex parichmacherMutex;
mutex zalaOzhidaniyaMutex;
queue<int> ozhidayushchieKlienty;
condition_variable cv;
boost::interprocess::interprocess_semaphore dostupnyeKresla(STULIJA);

void strizhkaVolosUParichmachera(int idKlienta) {
    lock_guard<mutex> lock(parichmacherMutex);
    cout << "Парикмахер стрижет волосы клиенту " << idKlienta << "." << endl;
    this_thread::sleep_for(chrono::seconds(4));
}

void parichmacherSpat() {
    unique_lock<mutex> lock(parichmacherMutex);
    cout << "Парикмахер спит." << endl;
    cv.wait(lock);
}

void pribytieKlienta(int idKlienta) {
    if (dostupnyeKresla.try_wait()) {
        lock_guard<mutex> lock(zalaOzhidaniyaMutex);
        ozhidayushchieKlienty.push(idKlienta);
        cout << "Клиент " << idKlienta << " занимает место в приемной." << endl;
        cv.notify_one();
    } else {
        cout << "Нет свободных мест. Клиент " << idKlienta << " облысел." << endl;
    }
}

void simulirovatPribytieKlientov() {
    for (size_t i = 0; i <= 20; ++i) {
        pribytieKlienta(i);

        this_thread::sleep_for(chrono::seconds(1));
        if (i > 7) {
            this_thread::sleep_for(chrono::seconds(1));
        }
        if(i == 15){
            this_thread::sleep_for(chrono::seconds(25));
        }
    }
}

void simulirovatParikmacherskuyu() {
    while (true) {
        unique_lock<mutex> lock(zalaOzhidaniyaMutex);

        if (!ozhidayushchieKlienty.empty()) {
            int id = ozhidayushchieKlienty.front();
            ozhidayushchieKlienty.pop();
            lock.unlock();

            cout << "Клиент " << id << " приглашен парикмахером." << endl;
            strizhkaVolosUParichmachera(id);
            dostupnyeKresla.post();
        } else {
            lock.unlock();
            parichmacherSpat();
        }
    }
}

int main(void) {
    thread pribytieThread(simulirovatPribytieKlientov);
    thread simulirovanieThread(simulirovatParikmacherskuyu);

    pribytieThread.join();
    simulirovanieThread.join();

    return 0;
}
