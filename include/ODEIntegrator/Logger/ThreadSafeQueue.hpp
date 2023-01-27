#pragma once

#include <queue>
#include <mutex>
#include <condition_variable>

template<typename T>
class ThreadSafeQueue {
    public:
        ThreadSafeQueue() {}
        void push(T item);
        void pop(T& item);

        // Delete copy constructors
        ThreadSafeQueue(const ThreadSafeQueue&) = delete;
        ThreadSafeQueue& operator=(const ThreadSafeQueue&) = delete;
    private:
        std::mutex mtx;
        std::condition_variable cond_var;
        std::queue<T> myqueue;
};
