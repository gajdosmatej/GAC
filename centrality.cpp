#include "centrality.h"

using namespace std;

void centr::centrality(int language){

  ifstream allFile;
  allFile.open("sorted.txt");
  ofstream centralityFile("centrality.txt");

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
  allFile.open("sorted.txt");

  int maxIndex = floor(size * (100 - min) / 100);
  int minIndex = ceil(size * (100 - max) / 100);

  centralityFile << min << "-" << max << "%\n";

  int i = 0;
  while(getline(allFile, line)){

    if( (i > minIndex) && (i < maxIndex) )  centralityFile << line << "\n";  //je v intervalu
    ++i;

  }

  centralityFile.close();

}
