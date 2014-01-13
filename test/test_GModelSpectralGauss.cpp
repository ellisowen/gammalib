#include <iostream>
#include "GammaLib.hpp"

void test_a() {
  GModelSpectralGauss model;
  GEnergy energy(1, "TeV");
  GTime time;
  std::cout << model.eval(energy, time) << std::endl;
}

void test_b() {
  GModelSpectralGauss model(42, GEnergy(3, "MeV"), GEnergy(1, "MeV"));
  GEnergy energy(3, "MeV");
  GTime time;
  std::cout <<  model.eval(energy, time) << std::endl;
  //std::cout <<  model.flux(GEnergy) << std::endl;
}


int main(void) {
  test_a();
  test_b();
}
