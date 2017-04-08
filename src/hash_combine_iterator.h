#ifndef HCI_HASH_COMBINER_ITERATOR_H_
#define HCI_HASH_COMBINER_ITERATOR_H_

#include <boost/functional/hash.hpp>
#include <cstddef>

namespace hci {

class HashCombineIterator {
 public:
  HashCombineIterator(std::size_t& seed) { this->seed = &seed; }
  void operator()(const std::size_t x) const { boost::hash_combine(*seed, x); }

 private:
  std::size_t* seed;
};
}

#endif