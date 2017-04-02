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
  dets.push_back(det);
  coefs.push_back(coef);
  n++;
  return dets.back();
}

void Wavefunction::load(const std::string& filename, const int n_orbs) {
  int n_in;
  int n_up, n_dn;
  double coef;
  const Det det_proto(n_orbs);
  int orb;
  
  // Read header line.
  std::ifstream wf_file(filename.c_str());
  wf_file >> n_in >> n_up >> n_dn;

  // Read each coef and det.
  for (int i = 0; i < n_in; i++) {
    wf_file >> coef;
    Det& det = append_det(det_proto, coef);
    for (int j = 0; j < n_up; j++) {
      wf_file >> orb;
      det.up.set_orb(orb - 1, true);
    }
    for (int j = 0; j < n_dn; j++) {
      wf_file >> orb;
      det.dn.set_orb(orb - 1, true);
    }
  }

  printf("Loaded wavefunction with %d dets.\n", n);
}

const std::list<hci::Det>& Wavefunction::get_dets() const {
  return dets;
}

const std::list<double>& Wavefunction::get_coefs() const {
  return coefs;
}

}