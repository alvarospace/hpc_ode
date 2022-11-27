#include "ODEIntegrator/Logger/ThreadSafeQueue.hpp"

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