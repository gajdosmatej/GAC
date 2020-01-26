#include "centrality.h"
#define SMERNICE 1  //smernice primek pro centralityLinear()
#define M_MAX 500 //maximalni hodnota multiplicity
#define POKLES 5  //o kolik ma y posun dolni primky klesat

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
  allFile.close();

}

int centr::getSpectators(string line){

  int c1 = 6;
  int c2 = 2;
  int spect1, spect2;

  for(int i = 0; i < c1; ++i)   line = line.substr(line.find(" ") + 1);
  spect1 = stoi(line.substr(0, line.find(" ")));

  for(int i = 0; i < c2; ++i) line = line.substr(line.find(" ") + 1);
  spect2 = stoi(line.substr(0, line.find(" ")));

  return spect1 + spect2;

}

int centr::getMultiplicity(string line){

  int c = 4;
  int M;

  for(int i = 0; i < c; ++i)   line = line.substr(line.find(" ") + 1);
  M = stoi(line.substr(0, line.find(" ")));

  return M;

}

int centr::getMinLine(float M1){

    ifstream sortedFile("sorted.txt");
    string line;

    int lineCounter = 0;

    while(getline(sortedFile, line)){

      int spect = getSpectators(line);
      int M = getMultiplicity(line);

      float M_line = SMERNICE * spect + M1;
      if(M_line > M)  return lineCounter;

      ++lineCounter;
    }

    sortedFile.close();
    return lineCounter;

}

int centr::getNumberEvents(float M1, float M2){

  ifstream filterFile("temp.txt");
  string line;

  int n = 0;

  while(getline(filterFile, line)){

    int M = getMultiplicity(line);
    int spect = getSpectators(line);

    float M_line = SMERNICE*spect + M2;

    if(M_line < M)  ++n;

  }

  filterFile.close();
  return n;
}

//vyfiltruj do temp.txt jen data pod vrchni primkou
void centr::makeTempFile(int minLine){

  ofstream tempFile("temp.txt");
  ifstream sortedFile("sorted.txt");
  string line;
  int i = 0;

  while(getline(sortedFile, line)){

    if(i >= minLine)  tempFile << line << "\n";
    ++i;
  }

  tempFile.close();
  sortedFile.close();
}

void centr::writeLinear(float M1, float M2){

  ifstream sortedFile("sorted.txt");
  ofstream centrFile("centrality.txt");
  string line;

  while(getline(sortedFile, line)){

      int spect = getSpectators(line);
      int M = getMultiplicity(line);

      int M_line1 = SMERNICE * spect + M1;
      int M_line2 = SMERNICE * spect + M2;

      if((M < M_line1) && (M > M_line2))  centrFile << line << "\n";

  }

  sortedFile.close();
  centrFile.close();
}

float centr::centralityLinear(int min, int max, int size){

    ofstream centralityFile("centrality.txt");
    ifstream allFile;

    allFile.open("sorted.txt");

    int maxIndex = floor(size * max / 100);
    int minIndex = ceil(size * min / 100);
    int n = maxIndex - minIndex;

    float M1; //posun po ose y pro horni primku

    //ziskej predpis horni primky
    if(min != 0) M1 = centralityLinear(0, min, size);
    else  M1 = M_MAX;


    float M2 = M1;
    makeTempFile(getMinLine(M1));

    int n0 = getNumberEvents(M1, M2);

    while(n0 < n){

     cout << M1 << " " << M2 << " - " << n0 << " " << n << endl;
     M2 -= POKLES;
     n0 = getNumberEvents(M1, M2); //zmensuj M2, dokud neni mezi primkami spravny pocet eventu
   }
    writeLinear(M1, M2);

    centralityFile.close();
    allFile.close();
    return M2;

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

  //centralityMultiplicity(min, max, size);
  //centralitySpectator(min, max, size);
  centralityLinear(min, max, size);

}
