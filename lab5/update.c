#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define HEADER_SIZE 128
#define FILE_DATA_SIZE 38 // 4

struct Block {
    char header[HEADER_SIZE];
    char data[FILE_DATA_SIZE];
    int size;
    int free;
    int num;
    struct Block* next;
};

struct Filesystem {
    struct Block* head;
};

void writeToBinaryFile(const char* fileName, struct Filesystem* fs){
    FILE* file = fopen(fileName, "wb");
    if(!file){
        fprintf(stderr, "Cannot open file\n");
        return;
    }

    struct Block* current = fs->head;
    while(current != NULL){
        fwrite(current, sizeof(struct Block), 1, file);
        current = current->next;
    }

    fclose(file);
}


void readFromBinaryFile(const char* fileName, struct Filesystem* fs){
    FILE* file = fopen(fileName, "rb");
    if(!file){
        fprintf(stderr, "Cannot open file\n");
        return;
    }

    struct Block* current = fs->head;
    struct Block* prev = NULL;

    while(fread(current, sizeof(struct Block), 1, file) == 1){
        if(prev != NULL){
            prev->next = current;
        }

        printf("Block %d - Header: %s,\n Data: \n%s\n", current->num, current->header, current->data);

        prev = current;
        current = (struct Block*)malloc(sizeof(struct Block));
    }

    if(prev != NULL){
        prev->next = NULL;
    }

    fclose(file);
}

struct Filesystem initFileSystem(){
    const char* fileName = "fs.bin";

    struct Filesystem fs;
    fs.head = NULL;

    writeToBinaryFile(fileName, &fs);

    return fs;
}

void createFile(struct Filesystem* fs, const char* filePath, const char* fileData){
    struct Block* current = fs->head;
    struct Block* prev = NULL;

    while(current != NULL){
        prev = current;
        current = current->next;
    }

    int fileNum = 0;
    if(prev != NULL){
        fileNum = prev->num + 1;
    }

    int dataLength = strlen(fileData);
    int remainingData = dataLength;

    int dataIndex = 0;

    while(remainingData > 0){
        struct Block* newBlock = (struct Block*)malloc(sizeof(struct Block));
        newBlock->num = fileNum;

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

        newBlock->next = NULL;

        if(prev == NULL){
            fs->head = newBlock;
        } else {
            prev->next = newBlock;
        }

        prev = newBlock;
        fileNum++;
    }

    writeToBinaryFile("fs.bin", fs);
}


void freeFileSystem(struct Filesystem* fs){
    struct Block* current = fs->head;
    struct Block* next;

    while(current != NULL){
        next = current->next;
        free(current);
        current = next;
    }
}

void moveFile(struct Filesystem* fs, const char* oldPath, const char* newPath){
    struct Block* current = fs->head;

    while(current != NULL){
        if(strcmp(current->header, oldPath) == 0){
            strncpy(current->header, newPath, HEADER_SIZE - 1);
            current->header[HEADER_SIZE - 1] = '\0';
        }
        current = current->next;
    }

    writeToBinaryFile("fs.bin", fs);
}

void deleteFile(struct Filesystem* fs, const char* filePath){
    struct Block* current = fs->head;
    struct Block* prev = NULL;
    struct Block* temp;

    while(current != NULL){
        if(strcmp(current->header, filePath) == 0){
            if(prev == NULL){
                temp = current;
                fs->head = current->next;
                current = fs->head;
                free(temp);
            } else {
                temp = current;
                prev->next = current->next;
                current = current->next;
                free(temp);
            }
        } else {
            prev = current;
            current = current->next;
        }
    }

    writeToBinaryFile("fs.bin", fs);
}

void copyFile(struct Filesystem* fs, const char* sourcePath, const char* destinationPath){
    struct Block* current = fs->head;

    while(current != NULL){
        if(strcmp(current->header, sourcePath) == 0){
            createFile(fs, destinationPath, current->data);
        }
        current = current->next;
    }
}

void viewFile(struct Filesystem* fs, const char* filePath){
    struct Block* current = fs->head;
    int find = 0;

    while(current != NULL){
        if(strcmp(current->header, filePath) == 0){
            find = 1;
            printf("File Content for %s:\n", filePath);

            while(current != NULL && strcmp(current->header, filePath) == 0){
                printf("%s", current->data);
                current = current->next;
            }
        }
        current = current->next;
    }

    if(!find){
        printf("File %s not found\n", filePath);
    }
}

void listFiles(struct Filesystem* fs, const char* directory){
    struct Block* current = fs->head;
    char prev_header[HEADER_SIZE];

    while(current != NULL){
        if(strstr(current->header, directory) == current->header){
            if(strcmp(current->header, prev_header) != 0){
                const char* fileName = current->header + strlen(directory);
                if(*fileName == '/'){
                    fileName++;
                }
                printf("%s\n", fileName);
                strcpy(prev_header, current->header);
            }
        }
        current = current->next;
    }
}

void listFileSystemDump(struct Filesystem* fs){
    struct Block* current = fs->head;
    char prev_header[HEADER_SIZE] = "";
    int indentLevel = 0;

    while(current != NULL){
        if(strcmp(current->header, prev_header) != 0){
            const char* fileName = current->header;
            int slashes = 0;
            for(size_t i = 0; i < strlen(fileName); ++i){
                if(fileName[i] == '/'){
                    slashes++;
                }
            }
            indentLevel = slashes - 1;
            for(size_t i = 0; i < indentLevel; ++i){
                printf("  ");
            }
            printf("%s\n", fileName);
            strcpy(prev_header, current->header);
        }
        current = current->next;
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

    readFromBinaryFile("fs.bin", &fs);

    // moveFile(&fs, "/root/folder1/file.txt", "/root/folder2/file.txt");

    // printf("\n\n\nAfter Moving:\n");
    // readFromBinaryFile("fs.bin", &fs);

    // deleteFile(&fs, "/root/folder1/file.txt");

    // printf("\n\n\nAfter Deleting:\n");
    // readFromBinaryFile("fs.bin", &fs);

    copyFile(&fs, "/root/folder1/file.txt", "/root/folder2/file_copy.txt");

    printf("\n\n\nAfter Copying:\n");
    readFromBinaryFile("fs.bin", &fs);

    printf("\n\n\nViewing File:\n");
    viewFile(&fs, "/root/folder1/file.txt");

    printf("\n\n\nListing Files in /root/folder1:\n");
    listFiles(&fs, "/root/folder1");

    printf("\n\n\nFile System Dump:\n");
    listFileSystemDump(&fs);

    freeFileSystem(&fs);

    return 0;
}

