#ifndef CENTRALITY_H
#define CENTRALITY_H

#include <iostream>
#include <fstream>
#include <math.h>

namespace centr{
  void centralityMultiplicity(int min, int max, int size);  //zapise do slozky centrality.txt jednu tridu centrality
  void centralitySpectator(int min, int max, int size);
  float centralityLinear(int min, int max, int size); //vraci posun po ose y v predpisu nizsi primky (nutne pro rekurzi)
  void start(int language);
  int getNumberEvents(float M1, float M2);
  int getSpectators(std::string line);
  int getMultiplicity(std::string line);
  void makeTempFile(int minLine);
  int getMinLine(float M1);
  void writeLinear(float M1, float M2);

}
#endif
