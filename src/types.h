#ifndef HCI_TYPES_H_
#define HCI_TYPES_H_

#include "det.h"

typedef std::pair<int, int> IntPair;
typedef std::array<int, 3> Int3;
typedef std::pair<Int3, Int3> Int3Pair;
typedef std::pair<Int3, double> Int3Double;
typedef std::pair<hci::Det, double> DetDouble;

#endif