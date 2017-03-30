#include "det.h"

#include <cstddef>
#include "spin_det.h"

namespace hci {
  
Det::Det(const Det& det) {
  this->up = det.up;
  this->dn = det.dn;
}

bool operator==(const Det& det1, const Det& det2) {
  return det1.up == det2.up && det1.dn == det2.dn;
}

Det& Det::operator=(const Det& det) {
  this->up = det.up;
  this->dn = det.dn;
  return *this;
}

void Det::resize(const int n_orb) {
  this->up.resize(n_orb);
  this->dn.resize(n_orb);
}

std::size_t hash_value(const Det& det) {
  return hash_value(det.up) ^ hash_value(det.dn);
}

}