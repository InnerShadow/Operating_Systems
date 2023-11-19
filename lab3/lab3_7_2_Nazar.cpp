#include <iostream>
#include <queue>
#include <pthread.h>
#include <semaphore.h>
#include <unistd.h>

const int STULIJA = 5;

using namespace std;

pthread_mutex_t parichmacherMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t zalaOzhidaniyaMutex = PTHREAD_MUTEX_INITIALIZER;
queue<int> ozhidayushchieKlienty;
sem_t dostupnyeKresla;
pthread_cond_t konditsiyaDostupnogoKresla = PTHREAD_COND_INITIALIZER;

void* strizhkaVolosUParichmachera(void* arg) {
    int idKlienta = *((int*)arg);
    pthread_mutex_lock(&parichmacherMutex);
    cout << "Парикмахер стрижет волосы клиенту " << idKlienta << "." << endl;
    pthread_mutex_unlock(&parichmacherMutex);
    sleep(4);
    return nullptr;
}

void pribytieKlienta(int idKlienta) {
    if (sem_trywait(&dostupnyeKresla) == 0) {
        pthread_mutex_lock(&zalaOzhidaniyaMutex);
        ozhidayushchieKlienty.push(idKlienta);
        cout << "Клиент " << idKlienta << " занимает место в приемной." << endl;
        pthread_mutex_unlock(&zalaOzhidaniyaMutex);
    } else {
        cout << "Нет свободных мест. Клиент " << idKlienta << " облысел." << endl;
    }
}

void* simulirovatPribytieKlientov(void* arg) {
    for (size_t i = 0; i <= 20; ++i) {
        pribytieKlienta(i);

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

void* simulirovatParikmacherskuyu(void* arg) {
    while (true) {
        pthread_mutex_lock(&zalaOzhidaniyaMutex);

        while (ozhidayushchieKlienty.empty()) {
            cout << "Парикмахер идет спать." << endl;
            pthread_mutex_unlock(&zalaOzhidaniyaMutex);
            sleep(2);

            pthread_mutex_lock(&zalaOzhidaniyaMutex);
        }

        int id = ozhidayushchieKlienty.front();
        ozhidayushchieKlienty.pop();
        pthread_mutex_unlock(&zalaOzhidaniyaMutex);

        cout << "Клиент " << id << " приглашен парикмахером." << endl;

        pthread_t parichmacherThread;
        pthread_create(&parichmacherThread, nullptr, strizhkaVolosUParichmachera, &id);
        pthread_join(parichmacherThread, nullptr);

        sem_post(&dostupnyeKresla);
    }
    return nullptr;
}

int main(void) {
    sem_init(&dostupnyeKresla, 0, STULIJA);

    pthread_t pribytieThread;
    pthread_create(&pribytieThread, nullptr, simulirovatPribytieKlientov, nullptr);

    pthread_t simulirovanieThread;
    pthread_create(&simulirovanieThread, nullptr, simulirovatParikmacherskuyu, nullptr);

    pthread_join(pribytieThread, nullptr);
    pthread_join(simulirovanieThread, nullptr);

    sem_destroy(&dostupnyeKresla);
    pthread_cond_destroy(&konditsiyaDostupnogoKresla);

    return 0;
}
