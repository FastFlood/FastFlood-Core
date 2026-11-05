#ifndef MODELDATAFETCHER_H
 #define MODELDATAFETCHER_H
 
 //to help with the dependency of data processing,
 //we define a helper class
 //that takes a lambda to to the download/fetching
 //a bunch of input data (of this same class type)
 //then, we can wait for the other datasets to be done,
 //automatically start this one, and maximize concurrency
 #include <memory>
 #include <future>
 #include <vector>
 #include <functional>
 #include <chrono>
 
 //
 #include "../../common/threading/threadpool.h"
 
 // Base class to allow heterogeneous collections of fetchers for dependency waiting
 class ModelDataFetcherBase {
 public:
     virtual ~ModelDataFetcherBase() = default;
     
     // Wait for this specific data to finish fetching/processing
     virtual void wait() const = 0;
     
     // Check if the fetching is done without blocking
     virtual bool is_ready() const = 0;
 };
 
 // Templated fetcher class to handle arbitrary types (e.g. cTMapf, ShapeFile, RainfallData)
 template <typename T>
 class ModelDataFetcher : public ModelDataFetcherBase {
 public:
     // Wrapping the result in a shared_ptr as requested
     using DataType = std::shared_ptr<T>;
 
     // Constructor for already-loaded/immediate data (no fetching required)
     explicit ModelDataFetcher(DataType data) {
         std::promise<DataType> promise;
         promise.set_value(std::move(data));
         future_data_ = promise.get_future().share();
     }
 
     // Constructor taking a thread pool, explicit dependencies, and a lambda for fetching
     template <typename Func>
     inline ModelDataFetcher(BS::thread_pool& pool,
                      const std::vector<std::shared_ptr<ModelDataFetcherBase>>& dependencies, 
                      Func&& fetch_lambda) 
     {
         future_data_ = pool.submit_task([deps = dependencies, f = std::forward<Func>(fetch_lambda)]() -> DataType {
             // 1. Wait for all prerequisite datasets to finish
             for (const auto& dep : deps) {
                 if (dep) {
                     dep->wait();
                 }
             }
             
             // 2. Execute this fetcher's lambda and return the result
             return f();
         }).share();
     }
 
     // Convenience constructor for zero dependencies
     template <typename Func>
     inline ModelDataFetcher(BS::thread_pool& pool, Func&& fetch_lambda)
         : ModelDataFetcher(pool, std::vector<std::shared_ptr<ModelDataFetcherBase>>{}, std::forward<Func>(fetch_lambda)) 
     {}
 
     inline void wait() const override {
         if (future_data_.valid()) {
             future_data_.wait();
         }
     }
 
     inline bool is_ready() const override {
         if (!future_data_.valid()) return false;
         return future_data_.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
     }
 
     // Retrieves the fetched data, blocking until it is ready
     inline DataType get() const {
         if (future_data_.valid()) {
             return future_data_.get();
         }
         return nullptr;
     }
 
 private:
     std::shared_future<DataType> future_data_;
 };
 
 #endif //MODELDATAFETCHER_H