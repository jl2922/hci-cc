#ifndef HCI_BIG_SPIN_DET_H_
#define HCI_BIG_SPIN_DET_H_

#include <cstdint>
#include <vector>

namespace hci {
typedef uint32_t Orbital;
typedef uint32_t EncodeBlock;
typedef std::vector<EncodeBlock> EncodeType;

class BigSpinDet {
  public:
    BigSpinDet() {};
    BigSpinDet(const int n_orbs) { this->n_orbs = n_orbs; }
    BigSpinDet(const BigSpinDet &bigSpinDet) {
      elecs = bigSpinDet.elecs;
      n_orbs = bigSpinDet.n_orbs;
    }
    friend bool operator==(const BigSpinDet &, const BigSpinDet &);
    BigSpinDet& operator=(const BigSpinDet& rhs) {
      elecs = rhs.elecs;
      n_orbs = rhs.n_orbs;
      return *this;
    }
    int get_n_elecs() const { return elecs.size(); }
    std::vector<int> get_elec_orbs(const int n_elecs) const {
      std::vector<int> elec_orbs;
      elec_orbs.reserve(elecs.size());
      for (const auto elec: elecs) elec_orbs.push_back(static_cast<int>(elec));
      return elec_orbs;
    }
    int get_n_orbs() const { return n_orbs; }
    void resize(const int n_orbs) { this->n_orbs = n_orbs; }
    BigSpinDet& set_orb(const int, const bool occ = true);
    bool get_orb(const int orb) const {
      for (const Orbital elec: elecs) {
        if (elec == static_cast<Orbital>(orb)) return true;
      }
      return false;
    }
    void from_eor(const BigSpinDet &, const BigSpinDet &);
    bool is_zero() const { return get_n_elecs() == 0; }
    EncodeType encode() const { return elecs; }
    void decode(const EncodeType & elecs, const int n_orbs) {
      this->elecs = elecs;
      this->n_orbs = n_orbs;
    }
    friend std::size_t hash_value(const BigSpinDet &);
  private:
    std::vector<Orbital> elecs;
    int n_orbs;
};

bool operator==(const BigSpinDet& lhs, const BigSpinDet& rhs);
std::size_t hash_value(const BigSpinDet &);

}
#endif
