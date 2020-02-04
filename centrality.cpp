#include "centrality.h"
#define SMERNICE 1  //smernice primek pro centralityLinear()
#define M_MAX 500 //maximalni hodnota multiplicity

using namespace std;

double POKLES = 100;  //o kolik ma y posun dolni primky klesat

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

int centr::getMinLine(double M1){

    ifstream sortedFile("sorted.txt");
    string line;

    int lineCounter = 0;

    while(getline(sortedFile, line)){

      int spect = getSpectators(line);
      int M = getMultiplicity(line);

      double M_line = SMERNICE * spect + M1;
      if(M_line > M)  return lineCounter;

      ++lineCounter;
    }

    sortedFile.close();
    return lineCounter;

}

int centr::getNumberEvents(double M1, double M2){

  ifstream filterFile("sorted.txt");
  string line;

  int n = 0;

  while(getline(filterFile, line)){

    int M = getMultiplicity(line);
    int spect = getSpectators(line);

    double M_lineDown = SMERNICE*spect + M2;
    double M_lineUp = SMERNICE*spect + M1;

    //cout << M << " - " << M_lineUp << " - " << M_lineDown << "\n";
    if((M_lineDown < M) && (M_lineUp > M))  ++n;

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

void centr::writeLinear(double M1, double M2, int n){

  ifstream sortedFile("sorted.txt");
  ofstream centrFile("centrality.txt");
  string line;
  int counter = 0;

  while(getline(sortedFile, line)){

      int spect = getSpectators(line);
      int M = getMultiplicity(line);

      double M_line1 = SMERNICE * spect + M1;
      double M_line2 = SMERNICE * spect + M2;

      if((M < M_line1) && (M > M_line2)){

          ++counter;

          if(counter <= n)  centrFile << line << "\n";
      }
  }

  sortedFile.close();
  centrFile.close();
}


float centr::centralityLinear(int min, int max, int size){

    int maxIndex = floor(size * max / 100);
    int minIndex = ceil(size * min / 100);
    int n = maxIndex - minIndex;

    double M1; //posun po ose y pro horni primku

    //ziskej predpis horni primky
    if(min != 0) M1 = centralityLinear(0, min, size);
    else  M1 = M_MAX;


    double M2 = M1;
    //makeTempFile(getMinLine(M1));

    int n0 = getNumberEvents(M1, M2);


   int previousUp;
   int previousDown;
   int sameCounterUp = 0;
   int sameCounterDown = 0;

   //dokud neni relativne presne
   while(!( (n0 > (n - 200)) && (n0 < (n + 200)) )){

     while(n0 < n){

       M2 -= POKLES;
       n0 = getNumberEvents(M1, M2); //zmensuj M2, dokud neni mezi primkami spravny pocet eventu

       if(previousDown == n0)  ++sameCounterDown;

       cout << M1 << " " << M2 << " - " << n0 << " " << n << " pokles: " << POKLES << endl;

       previousDown = n0;

     }

     POKLES /= 2;

     if((sameCounterUp > 2) || (sameCounterDown > 2)) break;

     while(n0 > n){

       M2 += POKLES;
       n0 = getNumberEvents(M1, M2); //zmensuj M2, dokud neni mezi primkami spravny pocet eventu

       if(previousUp == n0)  ++sameCounterUp;

       cout << M1 << " " << M2 << " - " << n0 << " " << n << " pokles: " << POKLES << endl;

       previousUp = n0;
     }

     POKLES /= 2;

   }

    if(min != 0)  writeLinear(M1, M2, n);

    POKLES = 100;

    return M2;

}


void centr::start(int language){

  ifstream allFile;
  allFile.open("sorted.txt");

  int min, max;

  if(language == 1) cout << "Určit dle [M]ultiplicity, [S]pektátorů, nebo [K]ombinace?\n";
  else  cout << "Define from [M]ultiplicity, [S]pectators, or [C]ombination?\n";

  string tech;
  cin >> tech;

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

  if((tech == "M") || (tech == "m")) centralityMultiplicity(min, max, size);
  else if((tech == "S") || (tech == "s")) centralitySpectator(min, max, size);
  //centralityMultiplicity(min, max, size);
  //centralitySpectator(min, max, size);
  else if((tech == "C") || (tech == "c") || (tech == "K") || (tech == "k")) centralityLinear(min, max, size);
  else  cout << "X\n";

}
