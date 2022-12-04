#pragma once

// Abstract class that every writer should inherit
class Writer {
    public:
        // Pure virtual read function
        virtual void write() = 0;
};