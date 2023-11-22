#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {
    int exec = 4;

    int length = argc - 1;
    int *arr = new int[length];

    for (int i = 0; i < length; i++) {
        arr[i] = atoi(argv[i + 1]);
    }

    int *pid_arr = new int[length];

    pid_arr[0] = getpid();
    cout << "1'st process (ID" << getpid() << ") - master process (ID " << getppid() << ")\n";

    for (int i = 1; i < length; i++) {
        pid_arr[i] = 0; 
            if (getpid() == pid_arr[arr[i] - 1]) {
            pid_t child_pid = fork();

            if (child_pid == 0) {
                cout << "process (ID" << getpid() << ") - master process (ID " << getppid() << ")\n";
                pid_arr[i] = getpid();
            }
        }
    }

    wait(NULL);

     if(getpid() == pid_arr[exec - 1]){
         execlp("whoami", "whoami", "--version", NULL); 
    }

    wait(NULL);

        cout << "process (ID" << getpid() << ") - master process (ID " << getppid() << ") FINISHED\n";

    delete[] arr;
    delete[] pid_arr;

    system("sleep 60");

    return 0;
}