#include "sort.h"
using namespace std;



  srt::intAndString::intAndString(int sortInteger, string l){

    this->sortInt = sortInteger;
    this->line = l;

  }


int srt::partition(vector<srt::intAndString*> &array, int leftBound, int rightBound){

  int pivot = array[rightBound]->sortInt;

  int index = leftBound - 1;
  for(int i = leftBound; i < rightBound; ++i){

    //i postupuje doprava a kazdou hodnotu mensi nez pivot prohazuje co nejvice doleva
    if(array[i]->sortInt < pivot){

      //posun index; prohod hodnoty, pokud nejsou shodne
      if(array[++index] != array[i]){

        swap(array[index], array[i]);
      }
    }
  }

  //kde index skonci, tam je hranice a polozi se tam pivot
  swap(array[++index], array[rightBound]);

  return index;

}

void srt::quickSort(vector<srt::intAndString*> &array, int leftBound, int rightBound){

  if(leftBound < rightBound){ //overeni

    int parIndex = partition(array, leftBound, rightBound);

    quickSort(array, leftBound, parIndex - 1);
    quickSort(array, parIndex + 1, rightBound);

  }

}

void srt::copyFile(string inputFileName){

  ifstream dataFile;
  dataFile.open(inputFileName);
  ofstream sortingFile("sorted.txt");
  string line;

  getline(dataFile, line);  //ignoruj prvni radek

  while(getline(dataFile, line)){

    sortingFile << line << endl;

  }

  dataFile.close();
  sortingFile.close();

}

//separuj multiplicitu z celého řádku dat
int srt::getMultiplicity(string line){

  const int inRowNumber = 4; //od 0
  for(int i = 0; i < inRowNumber; ++i)  line = line.substr(line.find(" ") + 1);  //odstran sloupce pred multiplicitou

  line = line.substr(0, line.find(" ")); //odstran sloupce za multiplicitou
  return stoi(line);

}

int srt::getSpectator(string line){

  const int inRowNumber1 = 6; //od 0
  const int inRowNumber2 = 8;

  for(int i = 0; i < inRowNumber1; ++i)  line = line.substr(line.find(" ") + 1);  //odstran sloupce pred spektatory1
  int number1 = stoi( line.substr(0, line.find(" ")) ); //odstran sloupce za multiplicitou

  for(int i = inRowNumber1; i < inRowNumber2; ++i)  line = line.substr(line.find(" ") + 1);  //odstran sloupce pred spektatory2
  int number2 = stoi( line.substr(0, line.find(" ")) );

  return number1 + number2;

}

void srt::sortMultiplicity(){

    vector<srt::intAndString*> data(0);
    ifstream sortingFile;
    sortingFile.open("sorted.txt");
    string line;

    //prepis hodnoty ze sorted.txt do pole
    while(getline(sortingFile, line))   data.push_back( new srt::intAndString(getMultiplicity(line), line) );

    sortingFile.close();

    //quicksort
    int length = data.size() - 1;
    quickSort(data, 0, length);

    //zapis serazena data do sorted.txt
    ofstream sortingFileW;
    sortingFileW.open("sorted.txt");

    for(int i = 0; i < length; ++i){

      sortingFileW << data[i]->line << endl;

    }

    sortingFileW.close();

    //smaz data
    for(int i = 0; i < data.size(); ++i)  delete data[i];
}


void srt::sortSpectator(){

    vector<srt::intAndString*> data(0);
    ifstream sortingFile;
    sortingFile.open("sorted.txt");
    string line;

    //prepis hodnoty ze sorted.txt do pole
    while(getline(sortingFile, line)){

      data.push_back( new srt::intAndString(getSpectator(line), line) );

    }

    sortingFile.close();

    //quicksort
    int length = data.size() - 1;
    quickSort(data, 0, length);

    //zapis serazena data do sorted.txt
    ofstream sortingFileW;
    sortingFileW.open("sorted.txt");

    for(int i = 0; i < length; ++i){

      sortingFileW << data[i]->line << endl;

    }

    sortingFileW.close();

    //smaz data
    for(int i = 0; i < data.size(); ++i)  delete data[i];
}


void srt::start(int language){

  string dataFileName;

  if(language == 1) cout << "Název souboru s daty: ";
  else  cout << "Enter data file name: ";
  cin >> dataFileName;

  if(language == 1) cout << "Kopíruji data do sorted.txt ...\n";
  else  cout << "Copying data to sorted.txt ...\n";
  copyFile(dataFileName);

  if(language == 1) cout << "Seřazuji sorted.txt ...\n";
  else  cout << "Sorting sorted.txt ...\n";

  sortSpectator();
  //sortMultiplicity();

}
