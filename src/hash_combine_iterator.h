#ifndef HCI_HASH_COMBINER_ITERATOR_H_
#define HCI_HASH_COMBINER_ITERATOR_H_

#include <cstddef>
#include <boost/functional/hash.hpp>

namespace hci {

class HashCombineIterator {
  public:
    HashCombineIterator(std::size_t& seed) {
      this->seed = &seed;
    }
    void operator()(const std::size_t x) const {
      boost::hash_combine(*seed, x);
    }
  private:
    std::size_t* seed;
};

}

#endif