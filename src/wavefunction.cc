#include "wavefunction.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>
#include "det.h"

namespace hci {

Wavefunction::Wavefunction() {
  n = 0;
}

Det& Wavefunction::append_det(const Det& det, const double coef) {
  n++;
  dets.push_back(det);
  coefs.push_back(coef);
  return dets.back();
}

const std::list<hci::Det>& Wavefunction::get_dets() const {
  return dets;
}

const std::list<double>& Wavefunction::get_coefs() const {
  return coefs;
}

int Wavefunction::size() const {
  return n;
}

}