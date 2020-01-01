#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

int main(){

  ifstream allFile;
  allFile.open("sorted.txt");
  ofstream centralityFile("centrality.txt");

  int min, max;
  cout << "Min percentage: ";
  cin >> min;
  cout << "Max percentage: ";
  cin >> max;

  int size = 0;
  string line;
  while(getline(allFile, line)) ++size; //pocet radku

  allFile.close();
  allFile.open("sorted.txt");

  int minIndex = ceil(size * min / 100);
  int maxIndex = floor(size * max / 100);

  centralityFile << min << "-" << max << "%\n";

  int i = 0;
  while(getline(allFile, line)){

    if( (i > minIndex) && (i < maxIndex) )  centralityFile << line << "\n";  //je v intervalu
    ++i;

  }

  centralityFile.close();
  return 0;
}
