#ifndef SORT_H
#define SORT_H

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>

namespace srt{

  //trida pro razeni radku podle multiplicity (obsahuje obe hodnoty)
  class intAndString{

  public:
    int M;
    std::string line;

    intAndString(int mult, std::string l);
  };

  int partition(std::vector<intAndString*> &array, int leftBound, int rightBound); //vytvori dve oblasti podle vztahu k pivotu v rozsahu leftBound az rightBound
  void quickSort(std::vector<intAndString*> &array, int leftBound, int rightBound); //implementace quicksortu
  void copyFile(std::string inputFileName);  //zkopiruj data do slozky sorted.txt
  int getMultiplicity(std::string line); //ziskej multiplicitu z radku
  void sortMultiplicity();  //serad multiplicitu a zapis cele radky do sorted.txt
  void start(int language); //inicializace programu


}

#endif
