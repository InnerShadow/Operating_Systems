#include <iostream>
#include <fstream>
#include <cstring>
#include <execution>
#include <regex>
#include <algorithm>
#include <cmath>

#define HEADER_SIZE 128
#define FILE_DATA_SIZE 384
#define MAX_BLOCKS 25

struct Block {
    int size;
    bool free;
    char header[HEADER_SIZE];
    char data[FILE_DATA_SIZE];
};

struct Filesystem {
    Block blocks[MAX_BLOCKS];

    void writeToBinaryFile(const char* fileName) {
        std::ofstream file(fileName, std::ios::binary | std::ios::out);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file for writing");
        }

        file.write(reinterpret_cast<char*>(blocks), sizeof(Block) * MAX_BLOCKS);

        file.close();
    }

    void readFromBinaryFile(const char* fileName) {
        std::ifstream file(fileName, std::ios::binary | std::ios::in);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file for reading");
        }

        file.read(reinterpret_cast<char*>(blocks), sizeof(Block) * MAX_BLOCKS);

        for (int i = 0; i < MAX_BLOCKS; ++i) {
            std::cout << "Block " << i << " - Header: " << blocks[i].header << ",\n Data: \n" << blocks[i].data << "\n";
        }

        file.close();
    }
};

Filesystem initFileSystem() {
    const char* fileName = "fs.bin";
    Filesystem fs;

    memset(&fs, 0, sizeof(Filesystem));

    std::ofstream file(fileName, std::ios::binary | std::ios::out);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file");
    }

    for (int i = 0; i < MAX_BLOCKS; ++i) {
        fs.blocks[i].free = true;
        fs.blocks[i].size = sizeof(Block);
        file.write(reinterpret_cast<char*>(&fs.blocks[i]), sizeof(Block));
    }

    file.close();

    return fs;
}

int allocateBlocks(Filesystem& fs, int numBlocks) {
    int startIndex = -1;
    for (int i = 0; i < MAX_BLOCKS; ++i) {
        if (fs.blocks[i].free) {
            int availableBlocks = 1;
            int j = i + 1;
            while (j < MAX_BLOCKS && fs.blocks[j].free) {
                availableBlocks++;
                j++;
            }

            if (availableBlocks >= numBlocks) {
                startIndex = i;
                for (int k = 0; k < numBlocks; ++k) {
                    fs.blocks[i + k].free = false;
                }
                fs.writeToBinaryFile("fs.bin");
                i += numBlocks - 1;
                break;
            }
        }
    }
    return startIndex;
}


void createFile(Filesystem& fs, const std::string& filePath, const std::string& fileData) {
    int numBlocks = static_cast<int>(ceil(static_cast<double>(fileData.size()) / FILE_DATA_SIZE));

    int startIndex = allocateBlocks(fs, numBlocks);

    if (startIndex != -1) {
        int dataPos = 0;
        for (int i = startIndex; i < startIndex + numBlocks; ++i) {
            strncpy(fs.blocks[i].header, filePath.c_str(), HEADER_SIZE - 1);
            fs.blocks[i].header[HEADER_SIZE - 1] = '\0';

            int dataSize = std::min(FILE_DATA_SIZE, static_cast<int>(fileData.size()) - dataPos);
            strncpy(fs.blocks[i].data, fileData.c_str() + dataPos, dataSize);
            fs.blocks[i].data[dataSize] = '\0';

            dataPos += dataSize;

            if (dataPos >= fileData.size()) {
                break;
            }

            if (i == startIndex) {
                dataPos = 0;
            }
        }
    } else {
        throw std::runtime_error("No free blocks available for file creation");
    }
}


void deleteFile(Filesystem& fs, const std::string& filePath) {
    for (int i = 0; i < MAX_BLOCKS; ++i) {
        if (!fs.blocks[i].free && strcmp(fs.blocks[i].header, filePath.c_str()) == 0) {
            int numBlocks = static_cast<int>(ceil(static_cast<double>(strlen(fs.blocks[i].data)) / FILE_DATA_SIZE));

            for (int j = i; j < i + numBlocks; ++j) {
                fs.blocks[j].free = true;
                memset(fs.blocks[j].header, 0, sizeof(fs.blocks[j].header));
                memset(fs.blocks[j].data, 0, sizeof(fs.blocks[j].data));
            }

            fs.writeToBinaryFile("fs.bin");

            return;
        }
    }

    throw std::runtime_error("File not found for deletion");
}

void copyFile(Filesystem& fs, const std::string& sourcePath, const std::string& destinationPath) {
    int sourceIndex = -1;
    for (int i = 0; i < MAX_BLOCKS; ++i) {
        if (!fs.blocks[i].free && strcmp(fs.blocks[i].header, sourcePath.c_str()) == 0) {
            sourceIndex = i;
            break;
        }
    }

    if (sourceIndex == -1) {
        throw std::runtime_error("Source file not found");
    }

    std::cout << "Copying from block " << sourceIndex << ": " << fs.blocks[sourceIndex].header << "\n";

    std::string sourceData = fs.blocks[sourceIndex].data;

    createFile(fs, destinationPath, sourceData);
}


void moveFile(Filesystem& fs, const std::string& sourcePath, const std::string& destinationPath) {
    copyFile(fs, sourcePath, destinationPath);
    deleteFile(fs, sourcePath);
}

int main() {
    try {
        Filesystem fs = initFileSystem();

        createFile(fs, "/root/folder1/file.txt",
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

        createFile(fs, "/root/folder1/folder3/file3.txt",
                   "Krabby Patty Recipe\n"
                   "Ingredients:\n"
                   "- 1 lb Krabby Patty secret formula (shhh!)\n"
                   "- Krabby Patty bun\n"
                   "- Fresh lettuce\n"
                   "- Tomato slices\n"
                   "- Pickles\n"
                   "- Ketchup and Mustard\n"
                   "- Cheese (optional)\n");

        createFile(fs, "/root/folder1/file1.txt",
                   "Beer Recipe\n"
                   "Ingredients:\n"
                   "- 4 lbs malted barley\n"
                   "- 1 oz hops\n"
                   "- 1 packet yeast\n"
                   "- 5 gallons water\n"
                   "- Priming sugar for bottling\n");

        createFile(fs, "/root/folder2/file2.txt",
                   "Healthy Relationship Recipe\n"
                   "Ingredients:\n"
                   "- Communication\n"
                   "- Trust\n"
                   "- Mutual Respect\n"
                   "- Empathy\n"
                   "- Quality Time\n"
                   "- Shared Goals\n"
                   "- Individual Space\n");

        fs.readFromBinaryFile("fs.bin");

        copyFile(fs, "/root/folder1/file.txt", "/root/folder2/copied_file.txt");

        std::cout << "\n\nAfter copying:\n";
        fs.readFromBinaryFile("fs.bin");

        moveFile(fs, "/root/folder1/file.txt", "/root/folder2/moved_file.txt");

        std::cout << "\n\nAfter moving:\n";
        fs.readFromBinaryFile("fs.bin");

        deleteFile(fs, "/root/folder2/moved_file.txt");

        std::cout << "\n\nAfter deletion:\n";
        fs.readFromBinaryFile("fs.bin");

    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
