#include "wavefunction.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include "det.h"

namespace hci {

Wavefunction::Wavefunction() {
  n = 0;
}

void Wavefunction::append_det(const Det& det, const double coef) {
  dets.push_back(det);
  coefs.push_back(coef);
  n++;
}

void Wavefunction::load(const std::string& filename, const int n_orbs) {
  int n_up, n_dn;
  double coef;
  int orb;
  
  // Read header line.
  std::ifstream wf_file(filename.c_str());
  wf_file >> n >> n_up >> n_dn;
  coefs.resize(n);
  dets.resize(n);

  // Read each coef and det.
  for (int i = 0; i < n; i++) {
    wf_file >> coef;
    coefs[i] = coef;
    Det& det = dets[i];
    det.resize(n_orbs);
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

const hci::Det& Wavefunction::get_det(const int idx) const {
  return dets[idx];
}

void Wavefunction::set_det(const int idx, const Det& det) {
  dets[idx] = det;
}

double Wavefunction::get_coef(const int idx) const {
  return coefs[idx];
}

void Wavefunction::set_coef(const int idx, const double coef) {
  coefs[idx] = coef;
}

}