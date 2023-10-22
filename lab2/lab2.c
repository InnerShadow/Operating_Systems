#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// Функция для создания процессов и сообщений о генеалогическом дереве
void create_process(int generations[], int num_generations, int current_generation) {
    int pid = getpid();
    int parent_pid = getppid();

    printf("Процесс с ID %d породил процесс с ID %d\n", parent_pid, pid);

    if (current_generation < num_generations) {
        for (int i = 0; i < generations[current_generation]; i++) {
            if (fork() == 0) {
                create_process(generations, num_generations, current_generation + 1);
                exit(0);
            }
        }
    }

    printf("Процесс с ID %d и ID родителя %d завершает работу\n", pid, parent_pid);
}

int main() {
    int generations[] = {0, 1, 1, 1, 3, 3, 5};
    int num_generations = sizeof(generations) / sizeof(generations[0]);

    int pid = getpid();
    int parent_pid = getppid();

    printf("Процесс с ID %d и ID родителя %d\n", pid, parent_pid);

    // Запускаем ls -h в процессе с ID 1
    if (pid == 1) {
        execlp("ls", "ls", "-h", NULL);
        perror("execlp");
    }

    create_process(generations, num_generations, 1);

    return 0;
}
