#include <iostream>
#include <pthread.h>
#include <math.h>
#include <memory>
#include <fstream>
#include <string>
#include <chrono>
#include <ctime>

pthread_mutex_t count_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t write_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t log_mutex = PTHREAD_MUTEX_INITIALIZER;

const char* data_file_path = "data.txt";
const char* log_file_path = "log.txt";

std::string count_time = "";
std::string write_time = "";

const int lim = 10;
const double step = 0.1;
const int n = 100000;

class Point {
private:
    double x;
    double y;

public:
    Point(double x, double y) : x(x), y(y) {}

    void set_x(double x) {
        this->x = x;
    }

    void set_y(double y) {
        this->y = y;
    }

    double get_x() const {
        return x;
    }

    double get_y() const {
        return y;
    }
};

std::shared_ptr<Point> point = std::make_shared<Point>(0, 0);

double f(double x) {
    double h = (x - step) / n;
    double sum = 0.0f;

    for (std::size_t i = 0; i < n; ++i) {
        double x_i = step + i * h;
        sum += std::exp(x_i * std::sin(x_i) * std::cos(x_i)) * std::cos(x_i) * std::pow(std::sin(x_i), 2) * std::cos(std::pow(x_i, 2)) * std::log(x_i) * x_i;;
    }

    return sum;
}

void* do_cout(void* arg) {
    for (double i = step; i < lim; i += step) {
        pthread_mutex_lock(&count_mutex);
        point->set_x(i);
        point->set_y(f(i));

        auto currentTime = std::chrono::system_clock::now();
        auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime.time_since_epoch());
        long long milliseconds_count = milliseconds.count();
        count_time = std::to_string(milliseconds_count);

        std::cout << "x = " << i << " y = " << f(i) << "\n";

        pthread_mutex_unlock(&write_mutex);
    }

    return nullptr;
}

void* do_write(void* arg) {
    for (std::size_t i = 0; i < lim * (1 / step); ++i) {
        pthread_mutex_lock(&write_mutex);

        std::ofstream f_data(data_file_path, std::ios::app);
        if (!f_data.is_open()) {
            std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAA\n";
            pthread_mutex_unlock(&log_mutex);
            return nullptr;
        }

        f_data << "x = " + std::to_string(point->get_x()) + " y = " + std::to_string(point->get_y()) + "\n";
        std::cout << "Write data in data file\n";

        f_data.close();

        auto currentTime = std::chrono::system_clock::now();
        auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime.time_since_epoch());
        long long milliseconds_count = milliseconds.count();
        write_time = std::to_string(milliseconds_count);

        pthread_mutex_unlock(&log_mutex);
    }

    return nullptr;
}

void* do_log(void* arg) {
    for (std::size_t i = 0; i < lim * (1 / step); ++i) {
        pthread_mutex_lock(&log_mutex);

        std::ofstream f_log(log_file_path, std::ios::app);
        if (!f_log.is_open()) {
            std::cout << "aaaaaaaaaaaaaaaaaaaaaaaaaaa\n";
            pthread_mutex_unlock(&count_mutex);
            return nullptr;
        }

        f_log << "count time = " + count_time + " write time = " + write_time + "\n";
        std::cout << "Write data in log file\n";

        f_log.close();

        pthread_mutex_unlock(&count_mutex);
    }

    return nullptr;
}

int main() {
    pthread_t thread1, thread2, thread3;

    pthread_mutex_lock(&log_mutex);
    pthread_mutex_lock(&write_mutex);

    pthread_create(&thread1, nullptr, do_cout, nullptr);
    pthread_create(&thread2, nullptr, do_write, nullptr);
    pthread_create(&thread3, nullptr, do_log, nullptr);

    pthread_join(thread1, nullptr);
    pthread_join(thread2, nullptr);
    pthread_join(thread3, nullptr);

    pthread_mutex_destroy(&count_mutex);
    pthread_mutex_destroy(&write_mutex);
    pthread_mutex_destroy(&log_mutex);

    return 0;
}
