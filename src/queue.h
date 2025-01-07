#ifndef QUEUE_H_
#define QUEUE_H_

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <iostream>
#include <queue>

template <typename T>
class safe_queue {
 public:
  safe_queue(int capacity) : capacity_(capacity) {
    pthread_mutex_init(&mutex_, NULL);
    pthread_cond_init(&cond_not_full_, NULL);
    pthread_cond_init(&cond_not_empty_, NULL);
  }

  ~safe_queue() {
    pthread_mutex_destroy(&mutex_);
    pthread_cond_destroy(&cond_not_full_);
    pthread_cond_destroy(&cond_not_empty_);
  }

 public:
  void push(const T& value) {
    /*  Lock */
    pthread_mutex_lock(&mutex_);
    /*  Judge datapool is full? */
    while (full()) {
      /*  Notify other producer stop push */
      pthread_cond_wait(&cond_not_full_, &mutex_);
    }
    /*  Push data */
    queue_.push(value);
    /*  Unlock */
    pthread_mutex_unlock(&mutex_);
    /*  Notify consumers : datapool isn't empty */
    pthread_cond_signal(&cond_not_empty_);
  }

  void pop(T& value) {
    /*  Lock */
    pthread_mutex_lock(&mutex_);
    /*  Judge datapool is empty? */
    while (empty()) {
      /*  Notify other consumer stop pop */
      pthread_cond_wait(&cond_not_empty_, &mutex_);
    }
    /*  Pop data */
    value = queue_.front();
    queue_.pop();
    /*  Unlock */
    pthread_mutex_unlock(&mutex_);
    /*  Notify producer : datapool isn't full */
    pthread_cond_signal(&cond_not_full_);
  }

 private:
  bool full() const { return (queue_.size() == capacity_); }

  bool empty() const { return queue_.empty(); }

 private:
  std::queue<T> queue_;
  size_t capacity_;
  pthread_cond_t cond_not_full_;
  pthread_cond_t cond_not_empty_;
  pthread_mutex_t mutex_;
};

#endif  // QUEUE_H_
