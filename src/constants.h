#ifndef HCI_CONSTANTS_H_
#define HCI_CONSTANTS_H_

#include <limits>

namespace hci {

class Constants {
  public:
    static constexpr double EPSILON = std::numeric_limits<double>::epsilon();
  private:
    Constants() {}
};

}

#endif