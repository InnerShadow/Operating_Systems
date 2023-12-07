#include <iostream>
#include <thread>
#include <mutex>
#include <cmath>
#include <memory>
#include <fstream>
#include <string>
#include <chrono>

std::mutex count_mutex;
std::mutex write_mutex;
std::mutex log_mutex;

const std::string data_file_path = "data.txt";
const std::string log_file_path = "log.txt";

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
    double sum = 0.0;

    for (std::size_t i = 0; i < n; ++i) {
        double x_i = step + i * h;
        sum += std::exp(x_i * std::sin(x_i) * std::cos(x_i)) * std::cos(x_i) * std::pow(std::sin(x_i), 2) * std::cos(std::pow(x_i, 2)) * std::log(x_i) * x_i;
    }

    return sum;
}

void do_cout() {
    for (double i = step; i < lim; i += step) {
        double x = i;
        double y = f(i);

        auto currentTime = std::chrono::system_clock::now();
        auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime.time_since_epoch());
        long long milliseconds_count = milliseconds.count();
        count_time = std::to_string(milliseconds_count);

        std::lock_guard<std::mutex> lock(count_mutex);

        point->set_x(x);
        point->set_y(y);

        write_mutex.unlock();

        std::cout << "x = " << x << " y = " << y << "\n";
    }
}

void do_write() {
    for (std::size_t i = 0; i < lim * (1 / step); ++i) {
        std::ofstream f_data(data_file_path, std::ios::app);
        if (!f_data.is_open()) {
            std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAA\n";
            log_mutex.unlock();
            return;
        }

        write_mutex.lock();

        f_data << "x = " + std::to_string(point->get_x()) + " y = " + std::to_string(point->get_y()) + "\n";

        auto currentTime = std::chrono::system_clock::now();
        auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime.time_since_epoch());
        long long milliseconds_count = milliseconds.count();
        write_time = std::to_string(milliseconds_count);

        log_mutex.unlock();

        std::cout << "Write data in data file\n";

        f_data.close();
    }
}

void do_log() {
    for (std::size_t i = 0; i < lim * (1 / step); ++i) {
        std::ofstream f_log(log_file_path, std::ios::app);
        if (!f_log.is_open()) {
            std::cout << "aaaaaaaaaaaaaaaaaaaaaaaaaaa\n";
            count_mutex.unlock();
            return;
        }

        log_mutex.lock();

        f_log << "count time = " + count_time + " write time = " + write_time + "\n";

        count_mutex.unlock();

        std::cout << "Write data in log file\n";

        f_log.close();
    }
}

int main() {
    std::thread thread1(do_cout);
    std::thread thread2(do_write);
    std::thread thread3(do_log);

    log_mutex.lock();
    write_mutex.lock();

    thread1.join();
    thread2.join();
    thread3.join();

    return 0;
}
