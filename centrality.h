#ifndef CENTRALITY_H
#define CENTRALITY_H

#include <iostream>
#include <fstream>
#include <math.h>

namespace centr{
  void centralityMultiplicity(int min, int max, int size);  //zapise do slozky centrality.txt jednu tridu centrality
  void centralitySpectator(int min, int max, int size);
  void centralityLinear(int min, int max, int size);
  void start(int language);
}
#endif
