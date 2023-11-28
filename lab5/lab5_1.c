#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define HEADER_SIZE 128
#define FILE_DATA_SIZE 38 // 4
#define MAX_BLOCKS 25

struct Block {
    char header[HEADER_SIZE];
    char data[FILE_DATA_SIZE];
    int size;
    int free;
    int num;
};

struct Filesystem {
    struct Block blocks[MAX_BLOCKS];
};

void writeToBinaryFile(const char* fileName, struct Filesystem* fs){
    FILE* file = fopen(fileName, "wb");
    if(!file){
        fprintf(stderr, "Cannot open file\n");
        return;
    }

    fwrite(fs->blocks, sizeof(struct Block), MAX_BLOCKS, file);

    fclose(file);
}

void readFromBinaryFile(const char* fileName, struct Filesystem* fs){
    FILE* file = fopen(fileName, "rb");
    if(!file){
        fprintf(stderr, "Cannot open file\n");
        return;
    }

    fread(fs->blocks, sizeof(struct Block), MAX_BLOCKS, file);

    for (int i = 0; i < MAX_BLOCKS; ++i) {
        printf("Block %d - Header: %s,\n Data: \n%s\n", fs->blocks[i].num, fs->blocks[i].header, fs->blocks[i].data);
    }

    fclose(file);
}

struct Filesystem initFileSystem(){
    const char* fileName = "fs.bin";

    struct Filesystem fs;

    for (int i = 0; i < MAX_BLOCKS; ++i) {
        fs.blocks[i].free = 1;
        fs.blocks[i].num = i;
    }

    writeToBinaryFile(fileName, &fs);

    return fs;
}

void createFile(struct Filesystem* fs, const char* filePath, const char* fileData){
    int fileNum = -1;

    for (int i = 0; i < MAX_BLOCKS; ++i) {
        if (fs->blocks[i].free) {
            fileNum = i;
            break;
        }
    }

    if (fileNum == -1) {
        fprintf(stderr, "No free blocks available.\n");
        return;
    }

    int dataLength = strlen(fileData);
    int remainingData = dataLength;

    int dataIndex = 0;

    while (remainingData > 0){
        struct Block* newBlock = &(fs->blocks[fileNum]);

        if (!newBlock->free) {
            for (int i = fileNum + 1; i < MAX_BLOCKS; ++i) {
                if (fs->blocks[i].free) {
                    fileNum = i;
                    newBlock = &(fs->blocks[fileNum]);
                    break;
                }
            }
        }

        strncpy(newBlock->header, filePath, HEADER_SIZE - 1);
        newBlock->header[HEADER_SIZE - 1] = '\0';

        int blockSize = FILE_DATA_SIZE;
        if(remainingData < FILE_DATA_SIZE){
            blockSize = remainingData;
        }

        strncpy(newBlock->data, fileData + dataIndex, blockSize);
        newBlock->size = blockSize;
        newBlock->free = 0;

        remainingData -= blockSize;
        dataIndex += blockSize;

        fileNum++;
    }

    writeToBinaryFile("fs.bin", fs);
}

void deleteFile(struct Filesystem* fs, const char* filePath){
    for (int i = 0; i < MAX_BLOCKS; ++i) {
        if (!fs->blocks[i].free && strcmp(fs->blocks[i].header, filePath) == 0) {
            fs->blocks[i].free = 1;
            memset(fs->blocks[i].header, 0, HEADER_SIZE);
            memset(fs->blocks[i].data, 0, FILE_DATA_SIZE);
            fs->blocks[i].size = 0;
        }
    }

    writeToBinaryFile("fs.bin", fs);
}

void moveFile(struct Filesystem* fs, const char* oldPath, const char* newPath){
    for (int i = 0; i < MAX_BLOCKS; ++i) {
        if (!fs->blocks[i].free && strcmp(fs->blocks[i].header, oldPath) == 0) {
            strncpy(fs->blocks[i].header, newPath, HEADER_SIZE - 1);
            fs->blocks[i].header[HEADER_SIZE - 1] = '\0';
        }
    }

    writeToBinaryFile("fs.bin", fs);
}

void copyFile(struct Filesystem* fs, const char* sourcePath, const char* destinationPath){
    int fileNum = -1;

    for (int i = 0; i < MAX_BLOCKS; ++i) {
        if (fs->blocks[i].free) {
            fileNum = i;
            break;
        }
    }

    if (fileNum == -1) {
        fprintf(stderr, "No free blocks available.\n");
        return;
    }

    for (int i = 0; i < MAX_BLOCKS; ++i) {
        if (!fs->blocks[i].free && strcmp(fs->blocks[i].header, sourcePath) == 0) {
            struct Block* newBlock = &(fs->blocks[fileNum]);

            if (!newBlock->free) {
                for (int j = fileNum + 1; j < MAX_BLOCKS; ++j) {
                    if (fs->blocks[j].free) {
                        fileNum = j;
                        newBlock = &(fs->blocks[fileNum]);
                        break;
                    }
                }
            }

            strncpy(newBlock->header, destinationPath, HEADER_SIZE - 1);
            newBlock->header[HEADER_SIZE - 1] = '\0';

            strncpy(newBlock->data, fs->blocks[i].data, FILE_DATA_SIZE);
            newBlock->size = fs->blocks[i].size;
            newBlock->free = 0;

            fileNum++;
        }
    }

    writeToBinaryFile("fs.bin", fs);
}

void viewFile(struct Filesystem* fs, const char* filePath){
    printf("Data of %s:", filePath);
    for (size_t i = 0; i < MAX_BLOCKS; ++i) {
        if (!fs->blocks[i].free && strcmp(fs->blocks[i].header, filePath) == 0) {
            printf("%s", fs->blocks[i].data);
        }
    }
    printf("\n\n");
}

void listFiles(struct Filesystem* fs, const char* directory){
    char prev_header[HEADER_SIZE] = "";

    printf("Listing Files in %s:\n", directory);

    for (size_t i = 0; i < MAX_BLOCKS; ++i) {
        if (!fs->blocks[i].free && strncmp(fs->blocks[i].header, directory, strlen(directory)) == 0) {
            const char* fileName = fs->blocks[i].header + strlen(directory);
            if (*fileName == '/') {
                fileName++;
            }
            if(strncmp(prev_header, fs->blocks[i].header, HEADER_SIZE) != 0){
                printf("%s\n", fileName);
                strcpy(prev_header, fs->blocks[i].header);
            }
        }
    }
}

void listFileSystemDump(struct Filesystem* fs){
    char prev_header[HEADER_SIZE] = "";
    int indentLevel = 0;

    printf("\n\n\nListing File System Dump:\n");

    for (size_t i = 0; i < MAX_BLOCKS; ++i) {
        if (!fs->blocks[i].free && strcmp(fs->blocks[i].header, prev_header) != 0) {
            const char* fileName = fs->blocks[i].header;
            int slashes = 0;
            for (size_t j = 0; j < strlen(fileName); ++j) {
                if (fileName[j] == '/') {
                    slashes++;
                }
            }
            indentLevel = slashes - 1;
            for (size_t j = 0; j < indentLevel; ++j) {
                printf("  ");
            }
            printf("%s\n", fileName);
            strcpy(prev_header, fs->blocks[i].header);
        }
    }
}


int main(void){
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

    createFile(&fs, "/root/folder1/folder3/file1.txt",
            "Beer Recipe\n"
            "Ingredients:\n"
            "- 4 lbs malted barley\n"
            "- 1 oz hops\n"
            "- 1 packet yeast\n"
            "- 5 gallons water\n"
            "- Priming sugar for bottling\n");

    deleteFile(&fs, "/root/folder1/file.txt");

    createFile(&fs, "/root/folder1/fil___e.txt",
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
            "- greens, salt, pepper to tasteAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n");

    readFromBinaryFile("fs.bin", &fs);

    //deleteFile(&fs, "/root/folder1/folder3/file1.txt");
    readFromBinaryFile("fs.bin", &fs);

    viewFile(&fs, "/root/folder1/file.txt");

    viewFile(&fs, "/root/folder2/copied_file.txt");

    listFiles(&fs, "/root/folder2");

    listFileSystemDump(&fs);

    return 0;
}
