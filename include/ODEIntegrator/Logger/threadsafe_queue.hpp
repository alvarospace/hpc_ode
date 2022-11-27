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

// Threadsafe push
template<typename T>
void ThreadSafeQueue<T>::push(T item) {
    std::lock_guard<std::mutex> lock(mtx);
    myqueue.push(item);
    cond_var.notify_one();
}

// Threadsafe pop waiting the queue until it has a item available
template<typename T>
void ThreadSafeQueue<T>::pop(T& item) {
    std::unique_lock<std::mutex> lock(mtx);
    cond_var.wait(lock, [this] {
        return !myqueue.empty();
    });
    item = myqueue.front();
    myqueue.pop();
}