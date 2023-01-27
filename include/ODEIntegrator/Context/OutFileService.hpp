#pragma once

#include <string>
#include <sstream>
#include <filesystem>
#include <chrono>

class OutFileService {
    public:
        OutFileService(std::string outFolder = DFT_OUTFOLDER) {
            // Create directory if it does not exist
            if (std::filesystem::exists(outFolder)){
                // Remove if the path exists and it's not a directory
                if (!std::filesystem::is_directory(outFolder)) {
                    std::filesystem::remove_all(outFolder);
                    std::filesystem::create_directory(outFolder);
                }
            } else {
                std::filesystem::create_directories(outFolder);
            }

            // Inside the out directory, create another with the date of the execution
            auto const nowTimeT = std::chrono::system_clock::to_time_t(
                std::chrono::system_clock::now()
            );
            std::stringstream formatedDate;
            formatedDate << std::put_time(std::localtime(&nowTimeT), "%F_%T");
            std::filesystem::path executionPath = std::filesystem::path(outFolder) / formatedDate.str();
            std::filesystem::create_directory(executionPath);

            executionFolder = executionPath.string();
        }

        std::string getExecutionFolder() {
            return executionFolder;
        }

    private:
        inline static std::string const DFT_OUTFOLDER = "./out/";

        std::string executionFolder;
};