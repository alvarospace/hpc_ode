#include <string>

#include "ODEIntegrator/Logger/BaseLogger.hpp"

#define info(x) info(x, __LINE__, __FILE__)
#define debug(x) debug(x, __LINE__, __FILE__)
#define error(x) error(x, __LINE__, __FILE__)

// TODO: Implement the loggers
// TODO: research about copy constructors in derived classes

class FileLogger : public BaseLogger {
    protected:
        void log(std::string data) override;

    public:
        FileLogger();
        ~FileLogger();

        // Delete copy constructors
        FileLogger(const FileLogger&) = delete;
        FileLogger& operator=(const FileLogger&) = delete;

    private:
        std::string logFileName;
};

class ConsoleLogger : public BaseLogger {
    protected:
        void log(std::string data) override;

    public:
        ConsoleLogger();
        ~ConsoleLogger();

        // Delete copy constructors
        ConsoleLogger(const ConsoleLogger&) = delete;
        ConsoleLogger& operator=(const ConsoleLogger&) = delete;
};