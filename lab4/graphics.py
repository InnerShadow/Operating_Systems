import matplotlib.pyplot as plt
import numpy as np

def main():
    N = 100

    with open(f"time_{N}.txt", 'r')  as f:
        data1 = [int(line.strip()) for line in f.readlines()]

    with open(f"CPP_time_{N}.txt", 'r')  as f:
        data2 = [int(line.strip()) for line in f.readlines()]

    fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize = (10, 5))

    axes[0].boxplot(data1)
    axes[0].set_title('MY ALOC')
    axes[0].set_ylim(min(data1 + data2) * 0.99, max(data1 + data2) * 1.01)

    axes[1].boxplot(data2)
    axes[1].set_title('CPP ALOC')
    axes[1].set_ylim(min(data1 + data2) * 0.99, max(data1 + data2) * 1.01)

    fig.suptitle('MY ALOC vs CPP ALOC')

    plt.show()


if __name__ == '__main__':
    main()

