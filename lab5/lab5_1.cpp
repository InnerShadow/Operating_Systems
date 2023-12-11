#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <algorithm>
#include <iterator>
#include <set>

const int FILE_SIZE = 2048 * 10;
const int BLOCK_DATA_SIZE = 10;
const int HEADER_SIZE = 128;

const std::string FSNAME = "filesystem.bin";

struct index {
    std::size_t id;
    bool free;
};

struct fileSystem {
    std::list<index> list;
};

fileSystem fs;

void init(const std::string& filename, std::size_t fileSize) {
    std::ofstream file(filename, std::ios::binary | std::ios::out);
    if (file.is_open()) {
        file.seekp(fileSize - 1);
        file.put(0);

        file.close();
        std::cout << "File '" << filename << "' initialized successfully.\n\n";
    } else {
        std::cout << "Error: Unable to open file '" << filename << "' for writing.\n\n";
    }
}

void createFile(const std::string& filePath, const std::string& fileData) {
    std::size_t dataSize = fileData.size();
    std::size_t blocksNeeded = dataSize / BLOCK_DATA_SIZE + 1;

    std::size_t currentBlockId = 0;
    for (std::size_t i = 0; i < blocksNeeded; ++i) {
            auto it = std::find_if(std::begin(fs.list), std::end(fs.list), [](const index& idx) {
                return idx.free;
            });

        if (it != std::end(fs.list)) {
            currentBlockId = it->id;
            it->free = false;
        } else {
            currentBlockId = fs.list.back().id + BLOCK_DATA_SIZE + HEADER_SIZE;
            index newIndex = {currentBlockId, false};
            fs.list.push_back(newIndex);
        }

        std::ofstream binaryFile(FSNAME, std::ios::binary | std::ios::out | std::ios::in);
        if (binaryFile.is_open()) {
            binaryFile.seekp(currentBlockId);

            std::string header = filePath;
            header.resize(HEADER_SIZE, '\0');
            binaryFile.write(header.c_str(), HEADER_SIZE);

            std::size_t blockSize = std::min(static_cast<std::size_t>(BLOCK_DATA_SIZE), dataSize - i * BLOCK_DATA_SIZE);
            binaryFile.write(fileData.c_str() + i * BLOCK_DATA_SIZE, blockSize);

            binaryFile.close();
        } else {
            std::cout << "Error: Unable to open file 'filesystem.bin' for writing.\n\n";
            return;
        }
    }

    std::cout << "File '" << filePath << "' created successfully.\n\n";
}

void deleteFile(const std::string& filePath) {
    auto it = std::find_if(std::begin(fs.list), std::end(fs.list), [&](const index& idx) {
        std::ifstream binaryFile(FSNAME, std::ios::binary | std::ios::in);
        binaryFile.seekg(idx.id);
        std::string header(HEADER_SIZE, '\0');
        binaryFile.read(&header[0], HEADER_SIZE);

        return header.substr(0, filePath.size()) == filePath;
    });

    while (it != std::end(fs.list)) {
        it->free = true;

        it = std::find_if(++it, std::end(fs.list), [&](const index& idx) {
            std::ifstream binaryFile(FSNAME, std::ios::binary | std::ios::in);
            binaryFile.seekg(idx.id);
            std::string header(HEADER_SIZE, '\0');
            binaryFile.read(&header[0], HEADER_SIZE);

            return header.substr(0, filePath.size()) == filePath;
        });
    }

    std::cout << "File '" << filePath << "' deleted successfully.\n\n";
}

void readFile(const std::string& filePath) {
    std::cout << "READ " << filePath << " file.\n";

    auto it = std::find_if(std::begin(fs.list), std::end(fs.list), [&](const index& idx) {
        std::ifstream binaryFile(FSNAME, std::ios::binary | std::ios::in);
        binaryFile.seekg(idx.id);
        std::string header(HEADER_SIZE, '\0');
        binaryFile.read(&header[0], HEADER_SIZE);

        return header.substr(0, filePath.size()) == filePath;
    });

    if (it != std::end(fs.list)) {
        std::ifstream binaryFile(FSNAME, std::ios::binary | std::ios::in);

        while (it != std::end(fs.list)) {
            std::string data(BLOCK_DATA_SIZE, '\0');
            binaryFile.seekg(it->id + HEADER_SIZE);
            binaryFile.read(&data[0], BLOCK_DATA_SIZE);

            std::cout << data;

            it = std::find_if(++it, std::end(fs.list), [&](const index& idx) {
                std::ifstream binaryFile(FSNAME, std::ios::binary | std::ios::in);
                binaryFile.seekg(idx.id);
                std::string header(HEADER_SIZE, '\0');
                binaryFile.read(&header[0], HEADER_SIZE);

                return header.substr(0, filePath.size()) == filePath;
            });
        }

        std::cout << "\nFile '" << filePath << "' read successfully.\n\n";
    } else {
        std::cout << "Error: File '" << filePath << "' not found.\n\n";
    }
}

void listFiles(const std::string& directory) {
    std::cout << "LS of " << directory << ":\n";
    std::set<std::string> files;

    for (const auto& idx : fs.list) {
        std::ifstream binaryFile(FSNAME, std::ios::binary | std::ios::in);
        binaryFile.seekg(idx.id);
        std::string header(HEADER_SIZE, '\0');
        binaryFile.read(&header[0], HEADER_SIZE);

        size_t pos = header.find_last_of('/');
        std::string fileDir = header.substr(0, pos + 1);

        if (fileDir == directory) {
            files.insert(header.substr(pos + 1));
        }
    }

    for (const auto& file : files) {
        std::cout << file << "\n";
    }

    std::cout << "\n";
}

void moveFile(const std::string& oldPath, const std::string& newPath) {
    auto it = std::find_if(std::begin(fs.list), std::end(fs.list), [&](const index& idx) {
        std::ifstream binaryFile(FSNAME, std::ios::binary | std::ios::in);
        binaryFile.seekg(idx.id);
        std::string header(HEADER_SIZE, '\0');
        binaryFile.read(&header[0], HEADER_SIZE);

        return header.substr(0, oldPath.size()) == oldPath;
    });

    while (it != std::end(fs.list)) {
        std::ifstream binaryFile(FSNAME, std::ios::binary | std::ios::in);

        while (it != std::end(fs.list)) {
            std::string data(BLOCK_DATA_SIZE, '\0');
            binaryFile.seekg(it->id + HEADER_SIZE);
            binaryFile.read(&data[0], BLOCK_DATA_SIZE);

            std::ofstream outFile(FSNAME, std::ios::binary | std::ios::out | std::ios::in);
            outFile.seekp(it->id);

            std::string newHeader = newPath;
            newHeader.resize(HEADER_SIZE, '\0');
            outFile.write(newHeader.c_str(), HEADER_SIZE);
            outFile.write(data.c_str(), BLOCK_DATA_SIZE);

            it = std::find_if(++it, std::end(fs.list), [&](const index& idx) {
                std::ifstream binaryFile(FSNAME, std::ios::binary | std::ios::in);
                binaryFile.seekg(idx.id);
                std::string header(HEADER_SIZE, '\0');
                binaryFile.read(&header[0], HEADER_SIZE);

                return header.substr(0, oldPath.size()) == oldPath;
            });
        }
    }

    std::cout << "File '" << oldPath << "' moved to '" << newPath << "' successfully.\n\n";
}

void copyFile(const std::string& sourcePath, const std::string& destinationPath) {
    std::string fileData = "";

    auto it = std::find_if(std::begin(fs.list), std::end(fs.list), [&](const index& idx) {
        std::ifstream binaryFile(FSNAME, std::ios::binary | std::ios::in);
        binaryFile.seekg(idx.id);
        std::string header(HEADER_SIZE, '\0');
        binaryFile.read(&header[0], HEADER_SIZE);

        return header.substr(0, sourcePath.size()) == sourcePath;
    });

    if (it != std::end(fs.list)) {
        std::ifstream binaryFile(FSNAME, std::ios::binary | std::ios::in);

        while (it != std::end(fs.list)) {
            std::string data(BLOCK_DATA_SIZE, '\0');
            binaryFile.seekg(it->id + HEADER_SIZE);
            binaryFile.read(&data[0], BLOCK_DATA_SIZE);

            fileData += data;

            it = std::find_if(++it, std::end(fs.list), [&](const index& idx) {
                std::ifstream binaryFile(FSNAME, std::ios::binary | std::ios::in);
                binaryFile.seekg(idx.id);
                std::string header(HEADER_SIZE, '\0');
                binaryFile.read(&header[0], HEADER_SIZE);

                return header.substr(0, sourcePath.size()) == sourcePath;
            });
        }
    } else {
        std::cout << "Error: File '" << sourcePath << "' not found.\n\n";
    }

    createFile(destinationPath, fileData);

    std::cout << "File '" << sourcePath << "' copied to '" << destinationPath << "' successfully.\n\n";
}

int main(void) {
    init(FSNAME, FILE_SIZE);

    createFile("root/folder1/file1.txt", "file data. 1212))");
    createFile("root/folder1/file2.txt", "file data. 1212 13 12)))))");

    deleteFile("root/folder1/file1.txt");

    createFile("root/folder1/file4.txt", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

    readFile("root/folder1/file4.txt");

    listFiles("root/folder1/");

    moveFile("root/folder1/file4.txt", "root/folder2/moved_file.txt");

    copyFile("root/folder2/moved_file.txt", "root/copied_file.txt");

    return 0;
}

