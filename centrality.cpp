#include "centrality.h"

using namespace std;

void centr::centralityMultiplicity(int min, int max, int size){

    ofstream centralityFile("centrality.txt");
    ifstream allFile;

    allFile.open("sorted.txt");

    int maxIndex = floor(size * (100 - min) / 100);
    int minIndex = ceil(size * (100 - max) / 100);

    centralityFile << min << "-" << max << "%\n";

    string line;
    int i = 0;
    while(getline(allFile, line)){

      if( (i > minIndex) && (i < maxIndex) )  centralityFile << line << "\n";  //je v intervalu
      ++i;

    }

    centralityFile.close();

}


void centr::centralitySpectator(int min, int max, int size){

  ofstream centralityFile("centrality.txt");
  ifstream allFile;

  allFile.open("sorted.txt");

  int maxIndex = floor(size * max / 100);
  int minIndex = ceil(size * min / 100);

  centralityFile << min << "-" << max << "%\n";

  string line;
  int i = 0;
  while(getline(allFile, line)){

    if( (i > minIndex) && (i < maxIndex) )  centralityFile << line << "\n";  //je v intervalu
    ++i;

  }

  centralityFile.close();

}

void centr::centralityLinear(int min, int max, int size){

    ofstream centralityFile("centrality.txt");

}


void centr::start(int language){

  ifstream allFile;
  allFile.open("sorted.txt");

  int min, max;

  if(language == 1) cout << "Dolní hranice procent: ";
  else  cout << "Min percentage: ";
  cin >> min;

  if(language == 1) cout << "Horní hranice procent: ";
  else  cout << "Max percentage: ";
  cin >> max;

  int size = 0;
  string line;
  while(getline(allFile, line)) ++size; //pocet radku

  allFile.close();

  centralityMultiplicity(min, max, size);
  //centralitySpectator(min, max, size);

}
