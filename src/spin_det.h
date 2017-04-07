#ifndef HCI_SPIN_DET_H_
#define HCI_SPIN_DET_H_

#include <cstddef>
#include <forward_list>
#include <iostream>
#include <vector>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/access.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/function_output_iterator.hpp>
#include "archive_iterator.h"

namespace hci {

typedef uint32_t BitsBlock;

class SpinDet {
  friend class boost::serialization::access;
  friend std::ostream& operator<<(std::ostream&, const SpinDet&);
  public:
    SpinDet() { };
    SpinDet(const int n_orb);
    SpinDet(const SpinDet& spinDet);
    friend bool operator==(const SpinDet&, const SpinDet&);
    SpinDet& operator=(const SpinDet&);
    int get_n_elecs() const;
    std::vector<int> get_elec_orbs(const int n_elecs = -1) const;
    int get_n_orbs() const;
    void resize(const int);
    SpinDet& set_orb(const int, const bool occ = true);
    bool get_orb(const int) const;
    void from_eor(const SpinDet&, const SpinDet&);
    bool is_zero() const;
    void print();
    friend std::size_t hash_value(const SpinDet&);
  private:
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const {
      const auto& output_iterator =
          boost::make_function_output_iterator(
              ArchiveIterator<Archive, BitsBlock>(ar));
      to_block_range(orbs, output_iterator);
      // std::cout << "<<" << orbs << std::endl;
    }
    template<class Archive>
    void load(Archive & ar, const unsigned int version) {
      int n = orbs.num_blocks();
      // int n_orbs = orbs.size();
      // std::cout << n; exit(1);
      orbs.clear();
      BitsBlock block = 0;
      for (int i = 0; i < n; i++) {
        ar >> block;
        orbs.append(block);
      }
      // std::cout << ">>" << orbs << std::endl;
      // orbs.resize(n_orbs);
      if (n != 3) exit(1);
      // if (orbs.size() != 81) exit(1);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER();
    boost::dynamic_bitset<BitsBlock> orbs;
};

bool operator==(const SpinDet&, const SpinDet&);
std::size_t hash_value(const SpinDet&);

std::ostream& operator<<(std::ostream&, const SpinDet&);

}

#endif // HCI_SPIN_DET_H_
