#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

int main() {
    int arr[] = {0, 1, 1, 1, 3, 3, 5};
    int pid_arr[7] = { 0 };

    pid_arr[0] = getpid();
    printf("1'st process (PID %d) - dangen master process (PPID %d)\n", 1, getppid());

    int length = sizeof(arr) / sizeof(arr[0]);

    for (size_t i = 1; i < length; ++i) {
        if (getpid() == pid_arr[arr[i] - 1]) {
            pid_t child_pid = fork();

            if (child_pid == 0) {
                printf("Slave process (PID %d) - dangen master process (PPID %d)\n", getpid() - pid_arr[0] + 1, 
                    getppid() - pid_arr[0] + 1);
                pid_arr[i] = getpid();
            }
        }
    }

    wait(NULL);

    return 0;
}
