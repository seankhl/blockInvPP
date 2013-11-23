/* a very nice ThBarrier class found at http://stackoverflow.com/a/14360253 */

#include <atomic>
#include <condition_variable>
#include <mutex>

class ThBarrier
{
  public:
    ThBarrier(size_t P_in);
    bool sync();
  private:
    ThBarrier();
    const size_t P;
    
    std::atomic<size_t> step;
    std::atomic<size_t> waiting;
    std::condition_variable cond;
    std::mutex mutex;
};

