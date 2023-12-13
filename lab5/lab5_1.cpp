#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <algorithm>
#include <iterator>
#include <set>
#include <sstream>

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
            currentBlockId = !fs.list.empty() ? fs.list.back().id + BLOCK_DATA_SIZE + HEADER_SIZE : BLOCK_DATA_SIZE + HEADER_SIZE;
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

std::vector<std::string> listFiles(const std::string& directory) {
    std::cout << "LS of " << directory << ":\n";
    std::set<std::string> files;

    for (const auto& idx : fs.list) {
        std::ifstream binaryFile(FSNAME, std::ios::binary | std::ios::in);
        binaryFile.seekg(idx.id);
        std::string header(HEADER_SIZE, '\0');
        binaryFile.read(&header[0], HEADER_SIZE);

        size_t pos = header.find('/');
        while (pos != std::string::npos) {
            std::string fileDir = header.substr(0, pos + 1);
            if (fileDir == directory) {
                files.insert(header.substr(pos + 1));
                break;
            }
            pos = header.find('/', pos + 1);
        }
    }

    int i = 0;

    for (const auto& file : files) {
        std::cout << ++i << ": " << file << "\n";
    }

    std::cout << "\n";
    std::vector<std::string> ret(files.begin(), files.end());
    return ret;
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

void writeIndexList() {
    std::ofstream binaryFile(FSNAME, std::ios::binary | std::ios::out | std::ios::in);
    if (binaryFile.is_open()) {
        std::size_t indexListOffset = FILE_SIZE - 400;
        binaryFile.seekp(indexListOffset);

        for (const auto& idx : fs.list) {
            binaryFile.write(reinterpret_cast<const char*>(&idx), sizeof(index));
        }

        binaryFile.close();
    } else {
        std::cout << "Error: Unable to open file '" << FSNAME << "' for writing.\n\n";
    }
}

void readIndexList() {
    std::ifstream binaryFile(FSNAME, std::ios::binary | std::ios::in);
    if (binaryFile.is_open()) {
        std::size_t indexListOffset = FILE_SIZE - 400;
        binaryFile.seekg(indexListOffset);

        fs.list.clear();

        while (binaryFile.tellg() < FILE_SIZE) {
            index idx;
            binaryFile.read(reinterpret_cast<char*>(&idx), sizeof(index));
            fs.list.push_back(idx);
        }

        binaryFile.close();
    } else {
        std::cout << "Error: Unable to open file '" << FSNAME <<  "' for reading.\n\n";
    }
}

int main(void) {
    bool b = false;

    std::size_t chose = 0;
    std::cout << "1 - init new file system;\n2 - load exsisting\n";
    std::cin >> chose;
    switch (chose) {
        case (1) : {
            init(FSNAME, FILE_SIZE);
            break;
        }

        case (2) : {
            readIndexList();
            break;
        }
    }

    std::cout << "\n";

    while (true) {
        std::cout << "1 - create file;\n2 - delete file;\n3 - move file;\n4 - copy file;\n5 - read file;\n6 - ls;\n7 - quit\n";
        std::cin >> chose;
        
        switch (chose) {
            case (1) : {
                std::string fName = "";
                std::cout << "Enter file name: ";
                std::cin >> fName;

                std::string fData = "";
                std::cout << "Enter file data: ";
                std::cin >> fData;

                createFile(fName, fData);
                break;
            }

            case (2) : {
                std::cout << "Enter file num: \n"; 
                std::vector<std::string> files = listFiles("root/");
                int file = 0;
                std::cin >> file;

                deleteFile("root/" + files[file - 1]);
                break;
            }

            case (3) : {
                std::cout << "Enter file num: \n"; 
                std::vector<std::string> files = listFiles("root/");
                int file = 0;
                std::cin >> file;

                std::string newPath = "";
                std::cout << "Enter new path: ";
                std::cin >> newPath;
                moveFile("root/" + files[file - 1], newPath);
                break;
            }

            case (4) : {
                std::cout << "Enter file num: \n"; 
                std::vector<std::string> files = listFiles("root/");
                int file = 0;
                std::cin >> file;

                std::string newPath = "";
                std::cout << "Enter new file path: ";
                std::cin >> newPath;
                copyFile("root/" + files[file - 1], newPath);
                break;
            }

            case (5) : {
                std::cout << "Enter file num: \n"; 
                std::vector<std::string> files = listFiles("root/");
                int file = 0;
                std::cin >> file;

                readFile("root/" + files[file - 1]);
                break;
            }

            case (6) : {
                listFiles("root/");
                break;
            }

            case (7) : {
                std::cout << "Have a good day!\n";
                b = true;
            }

            default : {
                std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAa\n";
            }
        }

        if (b) {
            break;
        }
    }

    writeIndexList();

    return 0;
}

