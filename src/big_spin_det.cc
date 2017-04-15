#include "big_spin_det.h"

#include <cassert>
#include <boost/functional/hash.hpp>

namespace hci {

void BigSpinDet::from_eor(const BigSpinDet& lhs, const BigSpinDet& rhs) {
  // Find the orbitals where lhs and rhs differ from each other.
  // Store in ascending order.
  assert(lhs.n_orbs == rhs.n_orbs);
  n_orbs = lhs.n_orbs;
  const auto& lhs_elecs = lhs.elecs;
  const auto& rhs_elecs = rhs.elecs;
  const int lhs_size = lhs_elecs.size();
  const int rhs_size = rhs_elecs.size();
  int lhs_ptr = 0;
  int rhs_ptr = 0;
  elecs.clear();
  elecs.reserve(rhs_size + rhs_size);
  while (lhs_ptr < lhs_size && rhs_ptr < rhs_size) {
    if (lhs_elecs[lhs_ptr] < rhs_elecs[rhs_ptr]) {
      elecs.push_back(lhs_elecs[lhs_ptr]);
      lhs_ptr++;
    } else if (lhs_elecs[lhs_ptr] > rhs_elecs[rhs_ptr]) {
      elecs.push_back(rhs_elecs[rhs_ptr]);
      rhs_ptr++;
    } else {
      lhs_ptr++;
      rhs_ptr++;
    }
  }
  while (lhs_ptr < lhs_size) {
    elecs.push_back(lhs_elecs[lhs_ptr]);
    lhs_ptr++;
  }
  while (rhs_ptr < rhs_size) {
    elecs.push_back(rhs_elecs[rhs_ptr]);
    rhs_ptr++;
  }
}

BigSpinDet& BigSpinDet::set_orb(const int orb, const bool occ) {
  auto it = elecs.begin();
  const Orbital orb_orb = static_cast<Orbital>(orb);
  while (it != elecs.end() && *it < orb_orb) it++;
  if (occ == true) {
    if (it != elecs.end() && *it == orb_orb) return *this;
    elecs.insert(it, orb);
  } else {
    if (it == elecs.end() || *it > orb_orb) return *this;
    elecs.erase(it);
  }
  return *this;
}

std::size_t hash_value(const BigSpinDet& bigSpinDet) {
  std::size_t seed = 0;
  for (const auto elec: bigSpinDet.elecs) {
    boost::hash_combine(seed, elec);
  }
  return seed;
}

bool operator==(hci::BigSpinDet const& lhs, hci::BigSpinDet const& rhs) {
  return lhs.elecs == rhs.elecs && lhs.n_orbs == rhs.n_orbs;
}

}
