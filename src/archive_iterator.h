#ifndef HCI_ARCHIVE_ITERATOR_H_
#define HCI_ARCHIVE_ITERATOR_H_

#include <cstddef>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/functional/hash.hpp>

namespace hci {

template<class Archive, class T>
class ArchiveIterator {
  public:
    ArchiveIterator(Archive& ar) {
      this->ar = &ar;
      // cnt = 0;
    }
    void operator()(const T x) const {
      (*ar) << x;
      // cnt++;
    }
    int cnt;
  private:
    Archive* ar;
};

}

#endif