#include <iostream>
#include <string>
#include <cstdio>
#include "edda.h"
#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include "io/edda_vtk_reader.h"
#include "dataset/dataset.h"

using namespace std;
using namespace edda;

int main(int argc, char **argv)
{
  srand(time(NULL));  // random seeding

  cout << "Loading sample file" << endl;
  string filename = string(SAMPLE_DATA_PATH) + "/isabel_pressure_small.vti";

  // load data with random sampling
  shared_ptr<Dataset<Real> > dataset = loadEddaScalarDataset(filename, "");

  VECTOR3 pos;
  Real value;
  int i;

  pos = VECTOR3(10,10,10);
  dataset->at_phys(pos, value);
  cout << pos << ": " << value << endl;

  pos = VECTOR3(2.1,2.1,2.1);
  dataset->at_phys(pos, value);
  cout << pos << ": " << value << endl;

  value = dataset->at_comp(5, 5, 5);
  cout << "at_comp(5,5,5) : " << value << endl;

  return 0;
}