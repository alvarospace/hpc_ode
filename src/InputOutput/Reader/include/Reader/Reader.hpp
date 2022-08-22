#pragma once

// Abstract class that every reader should inherit
class Reader {
    public:
        // Pure virtual read function. Does not returns anything
        // because data is saved in the Singleton "Mesh"
        virtual void read() = 0;
};