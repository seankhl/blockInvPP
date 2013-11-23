/* a very nice ThBarrier class found at http://stackoverflow.com/a/14360253 */

#include "Barrier.h"

ThBarrier::ThBarrier (size_t P_in) 
  : P(P_in), step(0), waiting(0)
{
    // nothing to do here!
}
bool ThBarrier::sync()
{
    size_t s = step.load();
    if (waiting.fetch_add(1) == P-1) {
        std::lock_guard<std::mutex> lock(mutex);
        waiting.store(0);
        ++step;
        cond.notify_all();
        return true;
    }
    else {
        std::unique_lock<std::mutex> lock(mutex);
        cond.wait(lock, [&]{return step != s;});
        return false;
    }
}
