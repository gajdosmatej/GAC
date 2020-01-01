#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>
using namespace std;

class intAndString{

public:
  int M;
  string line;

  intAndString(int mult, string l){

    this->M = mult;
    this->line = l;

  }
};


/*vector<intAndString*> quickSort(vector<intAndString*> array, int leftBound, int rightBound){

  int pivotIndex = rightBound;
  int pivot = array[pivotIndex]->M;
  int leftIndex = leftBound;
  int rightIndex = rightBound;

  //dokud se indexy neprekroci, opakuj
  while(leftIndex < rightIndex){

    while((array[leftIndex]->M < pivot) && (leftIndex < rightBound)) ++leftIndex;  //posouvej se doprava, dokud jsou hodnoty mensi nez pivot
    while((array[rightIndex]->M > pivot) && (rightIndex > leftBound)) --rightIndex;

    //v tuto chvili jsou oba indexy na cislech, ktere jsou na opacne strane pole
    if(leftIndex < rightIndex){ //kontrola, jestli se indexy neprekrocily
      //iter_swap(array.begin() + leftIndex, array.begin() + rightIndex); //prohod hodnoty
      intAndString * temp = array[leftIndex];
      array[leftIndex] = array[rightIndex];
      array[rightIndex] = array[leftIndex];

    }

    ++leftIndex;
    --rightIndex;
  }

//cout << leftIndex << " " << leftBound << " " << rightIndex << " " << rightBound << "\n";
  //rekurze
  if(rightIndex > leftBound)  array = quickSort(array, leftBound, rightIndex);
  if(leftIndex < rightBound)  array = quickSort(array, leftIndex, rightBound);

  return array;

}*/

int partition(vector<intAndString*> &array, int leftBound, int rightBound){

  int pivot = array[rightBound]->M;

  int index = leftBound - 1;
  for(int i = leftBound; i < rightBound; ++i){

    //i postupuje doprava a kazdou hodnotu mensi nez pivot prohazuje co nejvice doleva
    if(array[i]->M < pivot){

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

void quickSort(vector<intAndString*> &array, int leftBound, int rightBound){

  if(leftBound < rightBound){ //overeni

    int parIndex = partition(array, leftBound, rightBound);

    quickSort(array, leftBound, parIndex - 1);
    quickSort(array, parIndex + 1, rightBound);

  }

}

void copyFile(string inputFileName){

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
int getMultiplicity(string line){

  const int inRowNumber = 4; //od 0
  for(int i = 0; i < inRowNumber; ++i)  line = line.substr(line.find(" ") + 1);  //odstran sloupce pred multiplicitou

  line = line.substr(0, line.find(" ")); //odstran sloupce za multiplicitou
  return stoi(line);

}

void sortMultiplicity(){

    vector<intAndString*> data(0);
    ifstream sortingFile;
    sortingFile.open("sorted.txt");
    string line;

    //prepis hodnoty ze sorted.txt do pole
    while(getline(sortingFile, line)){

      data.push_back( new intAndString(getMultiplicity(line), line) );

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


int main(){

  string dataFileName;

  cout << "Enter data file name: ";
  cin >> dataFileName;

  cout << "Writing to sorted.txt ..." << endl;
  copyFile(dataFileName);

  cout << "Sorting sorted.txt ..." << endl;
  sortMultiplicity();

  return 0;

}
