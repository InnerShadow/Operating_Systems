#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define HEADER_SIZE 128
#define FILE_DATA_SIZE 384
#define NUM_BLOCKS 25

struct Block {
    int size;
    int free;
    char header[HEADER_SIZE];
    char data[FILE_DATA_SIZE];
};

struct Filesystem {
    struct Block* blocks;
    int numBlocks;
};

void initBlocks(struct Filesystem* fs) {
    fs->blocks = (struct Block*)malloc(NUM_BLOCKS * sizeof(struct Block));
    fs->numBlocks = NUM_BLOCKS;

    for (size_t i = 0; i < fs->numBlocks; ++i) {
        fs->blocks[i].free = 1;
        fs->blocks[i].size = sizeof(struct Block);
        memset(fs->blocks[i].header, 0, sizeof(fs->blocks[i].header));
        memset(fs->blocks[i].data, 0, sizeof(fs->blocks[i].data));
    }
}

void expandFileSystem(struct Filesystem* fs, int newSize) {
    fs->blocks = (struct Block*)realloc(fs->blocks, newSize * sizeof(struct Block));

    for (size_t i = fs->numBlocks; i < newSize; ++i) {
        fs->blocks[i].free = 1;
        fs->blocks[i].size = sizeof(struct Block);
        memset(fs->blocks[i].header, 0, sizeof(fs->blocks[i].header));
        memset(fs->blocks[i].data, 0, sizeof(fs->blocks[i].data));
    }

    fs->numBlocks = newSize;
}

void writeToBinaryFile(const char* fileName, struct Filesystem* fs) {
    FILE* file = fopen(fileName, "wb");
    if (!file) {
        perror("Cannot open file for writing");
        return;
    }

    fwrite(fs->blocks, sizeof(struct Block), fs->numBlocks, file);

    fclose(file);
}

void readFromBinaryFile(const char* fileName, struct Filesystem* fs) {
    FILE* file = fopen(fileName, "rb");
    if (!file) {
        perror("Cannot open file for reading");
        return;
    }

    fread(fs->blocks, sizeof(struct Block), fs->numBlocks, file);

    for (size_t i = 0; i < fs->numBlocks; ++i) {
        printf("Block %zu - Header: %s,\n Data: \n%s\n", i, fs->blocks[i].header, fs->blocks[i].data);
    }

    fclose(file);
}

struct Filesystem initFileSystem() {
    const char* fileName = "fs.bin";

    struct Filesystem fs;
    initBlocks(&fs);

    FILE* file = fopen(fileName, "wb");
    if (!file) {
        perror("Cannot open file");
        return fs;
    }

    fwrite(fs.blocks, sizeof(struct Block), fs.numBlocks, file);

    fclose(file);

    return fs;
}

void createFile(struct Filesystem* fs, const char* filePath, const char* fileData) {
    for (size_t i = 0; i < fs->numBlocks; ++i) {
        if (fs->blocks[i].free) {
            int requiredBlocks = (strlen(fileData) + sizeof(struct Block) - 1) / sizeof(struct Block);

            if (requiredBlocks > fs->numBlocks) {
                int newSize = requiredBlocks * 2;
                expandFileSystem(fs, newSize);
            }

            strncpy(fs->blocks[i].header, filePath, HEADER_SIZE - 1);
            fs->blocks[i].header[HEADER_SIZE - 1] = '\0';

            strncpy(fs->blocks[i].data, fileData, FILE_DATA_SIZE - 1);
            fs->blocks[i].data[FILE_DATA_SIZE - 1] = '\0';

            fs->blocks[i].free = 0;

            writeToBinaryFile("fs.bin", fs);

            return;
        }
    }

    fprintf(stderr, "No free blocks available for file creation\n");
}

void deleteFile(struct Filesystem* fs, const char* filePath) {
    for (size_t i = 0; i < fs->numBlocks; ++i) {
        if (!fs->blocks[i].free && strcmp(fs->blocks[i].header, filePath) == 0) {
            memset(fs->blocks[i].header, 0, sizeof(fs->blocks[i].header));
            memset(fs->blocks[i].data, 0, sizeof(fs->blocks[i].data));

            writeToBinaryFile("fs.bin", fs);

            fs->blocks[i].free = 1;

            return;
        }
    }

    fprintf(stderr, "File not found for deletion\n");
}

void copyFile(struct Filesystem* fs, const char* sourcePath, const char* destinationPath) {
    int sourceIndex = -1;
    for (size_t i = 0; i < fs->numBlocks; ++i) {
        if (!fs->blocks[i].free && strcmp(fs->blocks[i].header, sourcePath) == 0) {
            sourceIndex = i;
            break;
        }
    }

    if (sourceIndex == -1) {
        fprintf(stderr, "Source file not found\n");
        return;
    }

    int destinationIndex = -1;
    for (size_t i = 0; i < fs->numBlocks; ++i) {
        if (fs->blocks[i].free) {
            destinationIndex = i;
            break;
        }
    }

    if (destinationIndex == -1) {
        fprintf(stderr, "No free blocks available for copying\n");
        return;
    }

    strncpy(fs->blocks[destinationIndex].header, destinationPath, HEADER_SIZE - 1);
    fs->blocks[destinationIndex].header[HEADER_SIZE - 1] = '\0';

    strncpy(fs->blocks[destinationIndex].data, fs->blocks[sourceIndex].data, FILE_DATA_SIZE - 1);
    fs->blocks[destinationIndex].data[FILE_DATA_SIZE - 1] = '\0';

    fs->blocks[destinationIndex].free = 0;

    writeToBinaryFile("fs.bin", fs);
}

void moveFile(struct Filesystem* fs, const char* sourcePath, const char* destinationPath) {
    copyFile(fs, sourcePath, destinationPath);

    deleteFile(fs, sourcePath);
}

int main() {
    struct Filesystem fs = initFileSystem();

    createFile(&fs, "/root/folder1/file.txt",
               "Borscht Recipe\n"
               "Ingredients:\n"
               "- 500g meat (beef/pork)\n"
               "- 2 potatoes\n"
               "- 1 onion\n"
               "- 2 beets\n"
               "- 1 carrot\n"
               "- 1/2 cabbage\n"
               "- 2 tomatoes\n"
               "- 3 tbsp tomato paste\n"
               "- greens, salt, pepper to taste\n");

    createFile(&fs, "/root/folder1/file1.txt",
               "Beer Recipe\n"
               "Ingredients:\n"
               "- 4 lbs malted barley\n"
               "- 1 oz hops\n"
               "- 1 packet yeast\n"
               "- 5 gallons water\n"
               "- Priming sugar for bottling\n");

    createFile(&fs, "/root/folder2/file2.txt",
               "Healthy Relationship Recipe\n"
               "Ingredients:\n"
               "- Communication\n"
               "- Trust\n"
               "- Mutual Respect\n"
               "- Empathy\n"
               "- Quality Time\n"
               "- Shared Goals\n"
               "- Individual Space\n");

    createFile(&fs, "/root/folder1/folder3/file3.txt",
               "Krabby Patty Recipe\n"
               "Ingredients:\n"
               "- 1 lb Krabby Patty secret formula (shhh!)\n"
               "- Krabby Patty bun\n"
               "- Fresh lettuce\n"
               "- Tomato slices\n"
               "- Pickles\n"
               "- Ketchup and Mustard\n"
               "- Cheese (optional)\n");

    readFromBinaryFile("fs.bin", &fs);

    copyFile(&fs, "/root/folder1/file.txt", "/root/folder2/copied_file.txt");

    printf("\n\nAfter copying:\n");
    readFromBinaryFile("fs.bin", &fs);

    moveFile(&fs, "/root/folder1/file.txt", "/root/folder2/moved_file.txt");

    printf("\n\nAfter moving:\n");
    readFromBinaryFile("fs.bin", &fs);

    deleteFile(&fs, "/root/folder2/moved_file.txt");

    printf("\n\nAfter deletion:\n");
    readFromBinaryFile("fs.bin", &fs);

    free(fs.blocks); // Освобождаем выделенную память

    return 0;
}
