#ifndef GLAUBER_H
#define GLAUBER_H

#include <vector>
#include <iostream>
#include <gsl/gsl_rng.h>
#include  <gsl/gsl_randist.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include <fstream>

class Nucleus;

class Constants{

public:
  const float a = 0.535;  //skin thickness
  const float nucleusR = 6.38;   //polomer jadra ve fm
  const float nucleonR = 0.5;   //polomer jadra ve fm
  const float A = 1.5;  //parametr pro vypocet sily pro opravu prekryvu
  const float B = 10;  //parametr pro vypocet sily pro opravu prekryvu
  const float maxB = 2*nucleusR - 1;  //maximalni srazkovy parametr
  const int tableLength = 118;  //pocet prvku v table[]
  const std::string exceptions[4] = {"The element does not exists", "The nucleon number is smaller than proton number", "Input is not a number", "Invalid number of coordinates in Map node (INTERNAL ERROR)"}; //chybove hlasky
  const std::string table[118] = {"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"};  //periodicka soustava prvku


};

class Generator{

public:
  Generator();  //konstruktor
  ~Generator(); //destruktor
  double gen();  //rovnomerne rozdeleni <0 ; 1>
  int getSeed();  //generuj nahodny seed
  double genLinear();  //linearni rozdeleni <0 ; 13>
  double genInNuc(); //rovnomerne rozdeleni <0.5 ; 13.5>
  double genWoodSaxon(float a, float R0);
  double genNormWoodSaxon(float a, float R0);
  int poisson(float M_average);
  double genPosition(double min, double max); //rovnomerne rozdeleni od min do max

private:
  gsl_rng * generator;

};

class Map{

public:
  Map();  //konstruktor
  void writeCoords(double x, double y, double z);  //zapis souradnice do mapy
  std::vector<double> getCoords(int x, int y, int z);  //vrat souradnice z daneho uzlu
  bool isColliding(double x, double y, double z);
  void deleteCoords(int x, int y, int z, int start, int length); //smaz length souradnic z uzlu od pozice start

private:
  std::vector<std::vector<std::vector<std::vector<double>>>> map;  //3d pole obsahujici pole souradnic Nukleonu


};

class Nucleon{

public:
  //souradnice
  double x;
  double y;
  double z;
  int index;  //poradi v poli this->parent->nucleons tohoto Nukleonu
  float isospin;  // 1/2 => proton; -1/2 => neutron
  Nucleus * parent; //Nucleus obsahujici tento Nucleon
  bool isImpact = false;  //true => uz se srazil

  Nucleon(float I, Nucleus * par, int num); //konstruktor
  std::vector<double> makeNewCoords(double r);

};


class Nucleus{

public:
  Map * map;  //ukazatel na Mapu obsahujici 3d pole s poli souradnic Nukleonu
  int filled = 0;   //pocet Nukleonu v tomto objektu
  std::vector<Nucleon*> nucleons;   //pole ukazatelu na Nukleony v tomto objektu
  std::vector<Nucleon*> problematic;  //pole Nukleonu, ktere koliduji
  int problematicCounter = 0;   //pocet problematickych Nukleonu (delka this->problematic)
  int  binImp = 0;  //pocet binarnich srazek
  int imp = 0;    //pocet srazek
  int X;  //protonove cislo
  int Z;  //nukleonove cislo
  int N;  //neutronove cislo

  Nucleus(std::string symbol, int num);  //konstruktor
  ~Nucleus(); //destruktor

  void outputNucleons();  //zapise souradnice Nukleonu v tomto objektu do coordinates.txt (POUZE TESTOVACI METODA)
  void outputRad();  //zapise souradnice Nukleonu od stredu tohoto objektu do rads.txt (POUZE TESTOVACI METODA)
  void outputImp(float b, float impacts, double M_average, int M, int NA, int NnA, int NB, int NnB); //zapise impacts (soucet srazenych v obou jadrech), this->binImp, b, M_average, M, Na, Nna, Nb, Nnb do impacts.txt

private:

  void protonNumber(std::string symbol);  //ziskej ze symbolu prvku jeho protonove cislo

  void nucleonNumber(int num);   //ziskej ze zadaneho nukleonoveho cisla neutronove cislo, zkontroluj, jestli je this->Z >= this->X

  void createNucleons();  //vytvor dany pocet neutronu a protonu (dle this->X a this->N)

};

#endif
