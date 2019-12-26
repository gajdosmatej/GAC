#include <iostream>
#include <fstream>
#include <vector>
#include <gsl/gsl_rng.h>
#include  <gsl/gsl_randist.h>
#include <chrono>
using namespace std;

float alpha;
float epsilon;
vector<vector<double>> table;  //pocet a prumer pro kazde mnozstvi srazenych

int seed(){

  //generate seed from current time
  time_t tt = chrono::system_clock::to_time_t(chrono::system_clock::now());
  tm t = *localtime(&tt);
  int seconds = t.tm_sec;
  int minutes = t.tm_min;
  int hours = t.tm_hour;

  return hours * minutes * seconds;

}

gsl_rng * init(){

  //create random number generator
  gsl_rng * generator = gsl_rng_alloc (gsl_rng_taus2);


  //initialize generator with seed
  gsl_rng_set(generator, seed());

  return generator;

}


void average(float n, double M){

  double ratio = M / n;

  table[n][0]++;
  //rekurze prumeru
  table[n][1] = ( (table[n][0] - 1) * table[n][1] + ratio ) / table[n][0];

}

void averageM(float n, double M){

  table[n][0]++;
  //rekurze prumeru
  table[n][1] = ( (table[n][0] - 1) * table[n][1] + M) / table[n][0];

}

void averageBin(float n, int bin){

  table[n][0]++;
  //rekurze prumeru
  table[n][1] = ( (table[n][0] - 1) * table[n][1] + bin) / table[n][0];

}

void writeToFile(){

  ofstream f;
  f.open("multiplicity.txt", ios::app);

  int j = 0;
  int size = table.size();
  while(j < size){

    double k = table[j][1];
    f << j << " " << k << endl;
    j++;

  }

  f.close();

}


int main(){

  gsl_rng * gen = init();

  string line;
  ifstream file ("impacts.txt");

  ofstream fileO;
  fileO.open ("multiplicityAll.txt");

  cout << "Enter alpha: ";
  cin >> alpha;

	cout << "Enter epsilon: ";
	cin >> epsilon;

  bool firstLine = true;
  float imp = 0;
  int bin;

  if (file.is_open())
  {
    while ( getline (file,line) )
    {

      //hlavička souboru -> důležité info z něj
      if(firstLine){

        line = line.substr(line.find(" ") + 1, string::npos);
        line = line.substr(line.find(" ") + 1, string::npos);

        int size = stoi(line.substr(0, line.find(" ")));

        table.resize(size + 1, vector<double>(2, 0));

        firstLine = false;

      //get impacts and bin impacts on line
      }else{

        imp = stoi(line.substr(0, line.find(" ")));

        line = line.substr(line.find(" ") + 1, string::npos);

        bin = stoi(line.substr(0, line.find(" ")));

      }

      if(imp != 0){
        double M_average = epsilon * (imp * (1 - alpha) + alpha * bin);
        float M = gsl_ran_poisson(gen, M_average);

        //VSECHNA DATA
        fileO << imp << " " << bin << " " << M_average << endl;

        average(imp, M);
        //averageM(imp, M);
        //average(imp, M_average);
        //averageM(imp, M_average);
        //averageBin(imp, bin);

      }
    }

    file.close();

    //delete content
    ofstream file;
    file.open ("multiplicity.txt");
    file << "";
    file.close();
    fileO.close();

    writeToFile();

  }
  else{

    cout << "Missing data! (Unable to find impacts.txt)";

  }

  return 0;

}
