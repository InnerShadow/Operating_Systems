#include <stdio.h>
#include <string.h>

#define HEADER_SIZE 128
#define FILE_DATA_SIZE 384
#define NUM_BLOCKS 25

struct Block {
    char header[HEADER_SIZE];
    char data[FILE_DATA_SIZE];
    int size;
    int free;
};

struct Filesystem {
    struct Block blocks[NUM_BLOCKS];
};

void writeToBinaryFile(const char* fileName, struct Filesystem* fs) {
    FILE* file = fopen(fileName, "wb");
    if (!file) {
        fprintf(stderr, "Cannot open file\n");
        return;
    }

    fwrite(fs->blocks, sizeof(struct Block), NUM_BLOCKS, file);

    fclose(file);
}

void readFromBinaryFile(const char* fileName, struct Filesystem* fs) {
    FILE* file = fopen(fileName, "rb");
    if (!file) {
        fprintf(stderr, "Cannot open file\n");
        return;
    }

    fread(fs->blocks, sizeof(struct Block), NUM_BLOCKS, file);

    for (size_t i = 0; i < NUM_BLOCKS; ++i) {
        printf("Block %zu - Header: %s,\n Data: \n%s\n", i, fs->blocks[i].header, fs->blocks[i].data);
    }

    fclose(file);
}

            struct Filesystem initFileSystem() {
                const char* fileName = "fs.bin";

    struct Filesystem fs;

    FILE* file = fopen(fileName, "wb");
    if (!file) {
        fprintf(stderr, "Cannot open file\n");
        return fs;
    }

    for (size_t i = 0; i < NUM_BLOCKS; ++i) {
        fs.blocks[i].free = 1;
        fs.blocks[i].size = sizeof(struct Block);
        memset(fs.blocks[i].header, 0, sizeof(fs.blocks[i].header));
        memset(fs.blocks[i].data, 0, sizeof(fs.blocks[i].data));
        fwrite(&fs.blocks[i], sizeof(struct Block), 1, file);
    }

    fclose(file);

    return fs;
}

void createFile(struct Filesystem* fs, const char* filePath, const char* fileData) {
    for (size_t i = 0; i < NUM_BLOCKS; ++i) {
        if (fs->blocks[i].free) {
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
    for (size_t i = 0; i < NUM_BLOCKS; ++i) {
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
    for (size_t i = 0; i < NUM_BLOCKS; ++i) {
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
    for (size_t i = 0; i < NUM_BLOCKS; ++i) {
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

void dumpFilesystem(const char* dumpFileName, struct Filesystem* fs) {
    FILE* dumpFile = fopen(dumpFileName, "w");
    if (!dumpFile) {
        fprintf(stderr, "Cannot open file\n");
        return;
    }

    for (size_t i = 0; i < NUM_BLOCKS; ++i) {
        int freeBytes = fs->blocks[i].size - (strlen(fs->blocks[i].header) + strlen(fs->blocks[i].data));

        fprintf(dumpFile, "Block %zu - Size: %d\n", i, fs->blocks[i].size);
        fprintf(dumpFile, "  Header: %s\n", fs->blocks[i].header);
        fprintf(dumpFile, "  Data:\n%s\n", fs->blocks[i].data);
        fprintf(dumpFile, "  Free Bytes: %d\n", freeBytes);
    }

    fclose(dumpFile);
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
               "- 1 lb Krabby Patty secret formula\n"
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

    dumpFilesystem("filesystem_dump.txt", &fs);

    return 0;
}

