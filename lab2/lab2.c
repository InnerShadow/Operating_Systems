#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {

    int exec = 1;

    if (argc < 2) {
        fprintf(stderr, "Usage: %s <arr[0]> <arr[1]> ... <arr[n-1]>\n", argv[0]);
        return 1;
    }

    int length = argc - 1;
    int *arr = (int *)malloc(length * sizeof(int));

    for (size_t i = 0; i < length; ++i) {
        arr[i] = atoi(argv[i + 1]);
    }

    int *pid_arr = (int *)malloc(length * sizeof(int));

    pid_arr[0] = getpid();
    printf("\n1'st process (ID %d) - dangen master process (ID %d)\n", 1, getppid());

    for (int i = 1; i < length; ++i) {
        if (getpid() == pid_arr[arr[i] - 1]) {
            pid_t child_pid = fork();

            if (child_pid == 0) {
                printf("Slave process (ID %d) - dangen master process (Master ID %d)\n", 
                    getpid() - pid_arr[0] + 1, getppid() - pid_arr[0] + 1);
                pid_arr[i] = getpid();
            }
        }
    }

    wait(NULL);

    if(getpid() == pid_arr[exec]){
        execlp("ls", "ls", "-h", NULL);
        printf("\n");   
    }

    wait(NULL);

    if(getppid() - pid_arr[0] >= 0){
        printf("Process %d with master PID %d finished.\n", 
            getpid() - pid_arr[0] + 1, getppid() - pid_arr[0] + 1);
    } else {
        printf("Process %d with master PID %d finished.\n", 
            getpid() - pid_arr[0] + 1, getppid());
    }

    free(arr);
    free(pid_arr);

    return 0;
}
