#include "wavefunction.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include "det.h"

namespace hci {

Wavefunction::Wavefunction() {
  this->n = 0;
}

void Wavefunction::load(std::string filename) {
  int n;
  int n_up, n_dn;
  double coef;
  int orb;
  int n_orb = 19; // TODO: From wf file.
  
  // Read header line.
  std::ifstream wf_file(filename.c_str());
  wf_file >> n >> n_up >> n_dn;
  this->n = n;
  this->coefs.resize(n);
  this->dets.resize(n);

  // Read each coef and det.
  for (int i = 0; i < n; i++) {
    wf_file >> coef;
    this->coefs[i] = coef;
    Det& det = this->dets[i];
    det.resize(n_orb);
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

const hci::Det& Wavefunction::get_det(int idx) const {
  return this->dets[idx];
}

double Wavefunction::get_coef(int idx) const {
  return this->coefs[idx];
}

}