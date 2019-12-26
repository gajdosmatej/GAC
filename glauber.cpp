#include <cmath>
#include <chrono>
#include "glauber.h"
using namespace std;

//----------------
//GLOBAL VARIABLES
//----------------
Generator * generator;
Constants * konst;
ofstream impactsFile;

//---------------------
//FUNCTION DECLARATIONS
//---------------------
int symbolToNumber(string symbol);
bool check(Nucleon * nuc);
void smallestR(Nucleus * nuc);

//---------------
//CLASS FUNCTIONS
//---------------

//Class generator
//---------------------------------------------------------------------------------------------
Generator::Generator(){

  gsl_rng * g = gsl_rng_alloc (gsl_rng_taus);
  gsl_rng_set(g, this->getSeed());

  this->generator = g;

}

Generator::~Generator(){

  gsl_rng_free(this->generator);

}

//generuj z rovnomerneho rozdeleni v intervalu <0 ; 1>
float Generator::gen(){

   return gsl_rng_uniform(this->generator);

}

//vygeneruj seed z aktualniho casu
int Generator::getSeed(){

  time_t tt = chrono::system_clock::to_time_t(chrono::system_clock::now());
  tm t = *localtime(&tt);
  int seconds = t.tm_sec;
  int minutes = t.tm_min;
  int hours = t.tm_hour;

  return hours * minutes * seconds;

}

//generuj z linearniho rozdeleni v intervalu <0 ; maxB>
float Generator::genLinear(){

  return konst->maxB * sqrt(this->gen());

}

//generuj z rovnomerneho rozdeleni v intervalu <0.5 ; 13.5> (v intervalu velikosti jadra)
float Generator::genInNuc(){

    return 2*(konst->nucleusR - konst->nucleonR) * this->gen() + konst->nucleonR;

  }

//ziskej multiplicitu z Poissonova rozdeleni
int Generator::poisson(float M_average){

  return gsl_ran_poisson(this->generator, M_average);

}

//ziskej vzdalenost od centra jadra z Wood-Saxonova potencialu
float Generator::genWoodSaxon(float a, float R0){

  //normalizacni konstanta
  float N = 1 / (a * log( exp(R0 / a) + 1) );

  float y = -a * log( exp(-this->gen() / (N * a)) - exp( -R0 / a) );
  return y;

}

float Generator::genPosition(float min, float max){

  return this->gen() * (max - min) + min;

}


//Class Map
//---------------------------------------------------------------------------------------------

  //konstruktor Mapy
  Map::Map(){

    this->map = vector< vector< vector< vector<float> > > > (2*konst->nucleusR+1,vector<vector<vector<float>>>(2*konst->nucleusR+1,vector<vector<float>>(2*konst->nucleusR+1,vector<float>(0))));

  }

  //zapis do Mapy
  void Map::writeCoords(float x, float y, float z){

    //nejblizsi uzel
    int roundX = round(x);
    int roundY = round(y);
    int roundZ = round(z);

    this->map[roundX][roundY][roundZ].push_back(x);
    this->map[roundX][roundY][roundZ].push_back(y);
    this->map[roundX][roundY][roundZ].push_back(z);

  }

  //ziskej souradnice z Mapy
  vector<float> Map::getCoords(int x, int y, int z){

    return this->map[x][y][z];
  }

  //smaz souradnice z Mapy
  void Map::deleteCoords(int x, int y, int z, int start, int length){

    this->map[x][y][z].erase(this->map[x][y][z].begin() + start, this->map[x][y][z].begin() + start + length);

  }

//Class Nucleon
//---------------------------------------------------------------------------------------------

  //konstruktor Nukleonu
  Nucleon::Nucleon(float I, Nucleus * par, int num){

    this->index = num;
    this->parent = par;
    this->isospin = I;

   //vzdalenost od centra
   float r = generator->genWoodSaxon(konst->a, konst->nucleusR);

    if(r > konst->nucleusR){r = konst->nucleusR;}

    float one = generator->genPosition(konst->nucleusR - r, konst->nucleusR + r);

    //v jakem rozptylu muze two byt
    float border = sqrt(r * r - pow(one - konst->nucleusR, 2));

    float two = generator->genPosition(konst->nucleusR - border, konst->nucleusR + border);

    //threeR = three +- konst->nucleusR
    float threeR = sqrt(r * r - pow(one - konst->nucleusR, 2) - pow(two - konst->nucleusR, 2));

    //na 50% se this->z pricte k konst->nucleusR, na 50% se odecte
    float unit = generator->gen();
    if(unit > 0.5){ unit = 1; } else{ unit = -1; }
    float three = konst->nucleusR + unit * threeR;

    //nahodne pridel one, two, three k this->x, this->y, this->z
    float rand = generator->gen();
    if(rand > 0.6666){ this->x = one; this->y = two; this->z = three; }
    else if(rand < 0.3333){ this->y = one; this->z = two; this->x = three; }
    else{ this->z = one; this->x = two; this->y = three; }


    ofstream file;
    file.open ("coordinates.txt", ios::app);
    file << this->x << " " << this->y << " " << this->z << " " << r << endl;

//cout << r << endl;

//cout << this->x << " " << this->y << " " << this->z << endl;

    //STARA GENERACE
    /*this->x = generator->genInNuc();
    this->y = generator->genInNuc();
    this->z = generator->genInNuc();
    //vzdalenost od stredu jadra
    float r = sqrt(pow(this->x - konst->nucleusR, 2) + pow(this->y - konst->nucleusR, 2) + pow(this->z - konst->nucleusR, 2));

    //kontrola, jestli neni vzdalenost stredu vetsi nez polomer jadra
    while(r > (konst->nucleusR - konst->nucleonR)){

      //dokud je, opakuj generaci souradnic
      this->x = generator->genInNuc();
      this->y = generator->genInNuc();
      this->z = generator->genInNuc();

      r = sqrt(pow(x - konst->nucleusR, 2) + pow(y - konst->nucleusR, 2) + pow(z - konst->nucleusR, 2));

    }*/
    //zjisti, se kterymi Nukleony tento Nukleon koliduje
    this->collisions();

    //zapis souradnice tohoto Nukleonu do nejblizsiho uzlu Mapy this->parent->map
    this->parent->map->writeCoords(this->x, this->y, this->z);
  }


  //zjisti, jestli Nukleon koliduje s nejakym jinym Nukleonem
  void Nucleon::collisions(){

    //nevyjde zkouska -> zapise se do pole this->parent->problematic
    if(!check(this)){

      this->parent->problematic.push_back(this);
      this->parent->problematicCounter++;

    }


  }


//Class Nucleus
//---------------------------------------------------------------------------------------------

//destruktor
Nucleus::~Nucleus(){

  for(int i = 0; i < this->filled; i++){ delete this->nucleons[i]; }
  delete this->map;

}

//konstruktor
Nucleus::Nucleus(string symbol, int num){

    this->map = new Map();

    this->protonNumber(symbol);
    this->nucleonNumber(num);
    this->createNucleons();

    //dokud se nejake Nukleony prekryvaji, opakuj jejich posuny
    while(this->fix()){}

}

//posun Nukleonu v this->problematic
bool Nucleus::fix(){

/*  cout << this->problematicCounter << endl;
  cout << this->problematic[0]->x << endl;*/

  vector<float*> forces;  //pole sil na vsechny Nukleony
  bool moving = false;  //indikuje, jestli doslo k pohybu

  //projed cele pole this->problematic
  for(int c = 0; c < this->problematicCounter; c++){

    //problematicky Nukleon
    Nucleon nuc = *(this->problematic[c]);

    //uzel v Mape problematickeho Nukleonu
    int rX = round(nuc.x);
    int rY = round(nuc.y);
    int rZ = round(nuc.z);

    //pokud jsou zaokrouhlene souradnice v krajnich hodnotach, posunou se od nich (jinak dojde k chybe pri procitani okolnich hodnot pole)
    if(rX == 0){rX++;} if(rX == konst->nucleusR*2){rX--;}
    if(rY == 0){rY++;} if(rY == konst->nucleusR*2){rY--;}
    if(rZ == 0){rZ++;} if(rZ == konst->nucleusR*2){rZ--;}

//cout << rX << " " << rY << " " << rZ << endl;

    vector<float> coord;  //souradnice z mapy v uzlu kolem uzlu [rX, rY, rZ]
    int nucNum; //pocet Nukleonu v coord
    float * force = new float[3]; //sila, ktera bude pusobit na nuc

    force[0] = 0;
    force[1] = 0;
    force[2] = 0;

    //kotrola okolnich uzlu v Mape
    for(int i = -1; i < 2; i++){
      for(int j = -1; j < 2; j++){
        for(int k = -1; k < 2; k++){


          coord = this->map->getCoords(rX+i, rY+j, rZ+k);

          nucNum = coord.size() / 3;  //v coord jsou vzdy trojice souradnic -> pocet Nukleonu v uzlu je tretinovy delce pole

          //spocita se vzdalenost r od kazdeho Nukleonu v danem okolnim uzlu
          for(int m = 0; m < nucNum; m++){

            float r = sqrt(pow(coord[3*m] - nuc.x, 2) + pow(coord[3*m+1] - nuc.y, 2) + pow(coord[3*m+2] - nuc.z, 2));

            //pokud je jejich vzdalenost mensi nez dvojnasobek jejich polomeru, zapocni posun
            if((r < 2*konst->nucleonR) && (r != 0)){

              moving = true;  //dochazi k pohybovani

              //pricti silu k jiz existujici (pri prekryvu s vice Nukleony)
              force[0] += (nuc.x - coord[3*m]) / (konst->A * pow(r, 3) + konst->B) + (generator->gen() - 0.5) * 2;
              force[1] += (nuc.y - coord[3*m+1]) / (konst->A * pow(r, 3) + konst->B)  + (generator->gen() - 0.5) * 2;
              force[2] += (nuc.x - coord[3*m+2]) / (konst->A * pow(r, 3) + konst->B)  + (generator->gen() - 0.5) * 2;

            }
          }
        }
      }
    }

    //do pole vsech sil pridej silu na tento Nukleon
    forces.push_back(force);

  }

  //aplikuj sily
  for(int c = 0; c < this->problematicCounter; c++){

    float x = this->problematic[c]->x;
    float y = this->problematic[c]->y;
    float z = this->problematic[c]->z;

    //souradnice v Mape pred pohybem
    int xB = round(x);
    int yB = round(y);
    int zB = round(z);

    float tempX = x;

    //vypocet novych souradnic po aplikaci sily
    x += forces[c][0];
    y += forces[c][1];
    z += forces[c][2];

    //vzdalenost od stredu jadra
    float d = sqrt(pow(x - konst->nucleusR, 2) + pow(y - konst->nucleusR, 2) + pow(z - konst->nucleusR, 2));

    //kontrola, jestli pohyb nevyhodi z jadra
    if(d <= (konst->nucleusR - konst->nucleonR)){

      //aplikace sily
      this->problematic[c]->x += forces[c][0];
      this->problematic[c]->y += forces[c][1];
      this->problematic[c]->z += forces[c][2];

    }

    //smazani pouzite sily z pole vsech sil
    delete[] forces[c];

    //souradnice v Mape po pohybu


    //aktualizace mapy

    int index = 0;

    while(this->map->getCoords(xB, yB, zB)[index] != tempX){index++;} //vyhledani pozice stareho x v mape (->pozice zmeneneho Nukleonu)

    this->map->deleteCoords(xB, yB, zB, index, 3);  //smazani x, y, z posunuteho Nukleonu

    //zapsani novych souradnic do Mapy
    this->map->writeCoords(this->problematic[c]->x, this->problematic[c]->y, this->problematic[c]->z);


  }

  vector<Nucleon*> problematicTemp = this->problematic;
  int counterTemp = this->problematicCounter;

  this->problematicCounter = 0;

  //zkontroluj vsechny posunute, jestli s necim opet nekoliduji
  for(int u = 0; u < counterTemp; u++){

    problematicTemp[u]->collisions();

  }
//cout << moving;
  //navrat, jestli doslo k pohybu
  return moving;
}

//zapis vzdalenosti Nukleonu od jadra do rads.txt (POUZE TESTOVACI)
void Nucleus::outputRad(){

  ofstream file;
  file.open ("rads.txt", ios::app);

  float r;

  for(int i = 0; i < this->Z; i++){

      r = sqrt(pow(this->nucleons[i]->x - konst->nucleusR, 2) + pow(this->nucleons[i]->y - konst->nucleusR, 2) + pow(this->nucleons[i]->z - konst->nucleusR, 2));
      file << r << endl;

  }

  file.close();

}

//zapis souradnice Nukleonu do coordinates.txt (POUZE TESTOVACI)
void Nucleus::outputNucleons(){

  ofstream file;
  file.open ("coordinates.txt", ios::app);

  for(int i = 0; i < this->Z; i++){

    file << this->nucleons[i]->x << " ";
    file << this->nucleons[i]->y << " ";
    file << this->nucleons[i]->z << endl;

  }

  file.close();

}

//vypis pocet srazenych Nukleonu, pocet binarnich srazek a srazkovy parametr b do impacts.txt
void Nucleus::outputImp(float b, float impacts, double M_average, int M, int NA, int NnA, int NB, int NnB){


  impactsFile << impacts << " " << this->binImp << " " << b << " " << M_average << " " << M << " " << NA << " " << NnA<< " " << NB << " " << NnB << endl;

}

//ziskej protonove cislo
void Nucleus::protonNumber(string symbol){

  try{
    int output = symbolToNumber(symbol);  //ziskej protonove cislo ze znacky prvku

    //0 -> prvek nebyl nalezen
    if(output != 0){

      this->X = output;

    }
    else{

      throw 0;

    }
  }
  catch(int e){

    cout << "Error: " << konst->exceptions[e] << endl;
    cout << "Enter symbol again: ";

    string sym;
    cin >> sym;
    protonNumber(sym);   //repeat

  }

}

//proved kontrolu nukleonoveho cisla, zapis
void Nucleus::nucleonNumber(int input){

  while(true){

    //spatny uzivatelsky vstup
    if(cin.fail()){

      cin.clear();
      cin.ignore();
      cout << "Error: " << konst->exceptions[2] << endl;
      cout << "Enter nucleon number: ";
      cin >> input;

    }
    else{break;}

  }
  //zkontroluj, ze Z neni mensi nez X
  try{

    if( (input - this->X) < 0){
      throw 1;
    }
    else{

      this->Z = input;
      this->N = this->Z - this->X;

    }
  }
  catch(int e){

    cout << "Error: " << konst->exceptions[e] << endl;
    Nucleus::nucleonNumber(input);

  }

}

//vytvor Nukleony do jadra
void Nucleus::createNucleons(){

  for(int i = 0; i < this->Z; i++){

    if(i < this->X){
      this->nucleons.push_back(new Nucleon(0.5, this, filled));  //pridej Nukleon do pole this->nucleons
    }
    else{
      this->nucleons.push_back(new Nucleon(-0.5, this, filled));
    }

    filled++;
  }

}


//---------
//FUNCTIONS
//---------

//projdi vzdalenosti vsech Nukleonu a vypis tu nejmensi a nejvetsi (POUZE TESTOVACI)
void smallestR(Nucleus * n){

  int i = 0;
  float min = 100;
  float max = 0;

  //vsechny kombinace
  while(i < n->filled){

    int j = i+1;
    while(j < n->filled){

      float r = sqrt(pow(n->nucleons[i]->x - n->nucleons[j]->x, 2) + pow(n->nucleons[i]->y - n->nucleons[j]->y, 2) + pow(n->nucleons[i]->z - n->nucleons[j]->z, 2));
      if(r < min){min = r;}
      if(r > max){max = r;}
      j++;

    }
    i++;

  }

  cout << "Minimal distance: " << min << endl;
  cout << "Maximal distance: " << max << endl;

}

//zkontroluj, jestli Nukleon nekoliduje (true -> nekoliduje; false -> koliduje)
bool check(Nucleon * nuc){

  int x = round(nuc->x);
  int y = round(nuc->y);
  int z = round(nuc->z);

  if(x == 0){x++;} if(x == 2*konst->nucleusR){x--;}
  if(y == 0){y++;} if(y == 2*konst->nucleusR){y--;}
  if(z == 0){z++;} if(z == 2*konst->nucleusR){z--;}

  //blizke uzly k *nuc v Mape
  for(int i = -1; i < 2; i++){
    for(int j = -1; j < 2; j++){
      for(int k = -1; k < 2; k++){

        vector<float> coord = nuc->parent->map->getCoords(x+i, y+j, z+k);

        //pocet Nukleonu v uzlu (kazdy Nukleon zabira 3 pozice pole)
        int nucNum = coord.size() / 3;

        //zkontroluj vzdalenost s kazdym Nukleonem v danem uzlu
        for(int m = 0; m < nucNum; m++){

          float r = sqrt(pow(nuc->x - coord[m*3], 2) + pow(nuc->y - coord[m*3+1], 2) + pow(nuc->z - coord[m*3+2], 2));

          //dochazi k prekryvu -> vrat false
          if(r < 2*konst->nucleonR){
            return false;
          }
        }
      }
    }
  }
  return true;
}

//vrat protonove cislo podle znacky
int symbolToNumber(string symbol){

    int i = 0;
    bool found = false;

    //prohledej periodickou soustavu prvku (pole table)
    while(i < konst->tableLength){

      if(konst->table[i] == symbol){

        found = true;
        break;

      }
      i++;

    }

    if(found){
      return i + 1; //i starts at zero
    }
    else{
      return 0;
    }

}

//vytvor a sraz dve jadra (p1 a p2 znacky prvku, n1 a n2 nukleonova cisla, R polomer srazky ziskany z ucinneho prurezu, b srazkovy parametr)
void collide(string p1, int n1, string p2, int n2, float R, float b){

  Nucleus * nuc1 = new Nucleus(p1, n1);
  Nucleus * nuc2 = new Nucleus(p2, n2);

  for(int i = 0; i < nuc1->filled; i++){

    float x = nuc1->nucleons[i]->x;
    float y = nuc1->nucleons[i]->y + b; //posun po ose y o srazkovy parametr

    for(int j = 0; j < nuc2->filled; j++){

      float r = sqrt(pow(x - nuc2->nucleons[j]->x, 2) + pow(y - nuc2->nucleons[j]->y, 2));

      //sance srazky
      float p = exp(-(r/R)) / R;

      //aplikace sance srazky
      float rn = generator->gen();
      if(rn < p){

        //pokud se Nukleon jeste nesrazil, pricti pocet srazek do jeho jadra
        if(!nuc1->nucleons[i]->isImpact){nuc1->imp++;}
        if(!nuc2->nucleons[j]->isImpact){nuc2->imp++;}

        nuc1->nucleons[i]->isImpact = true;
        nuc2->nucleons[j]->isImpact = true;
        nuc1->binImp++; //pricti jadru binarni srazku
      }
    }
  }

  float epsilon = 0.25;
  float alpha = 0.145;

  //trefene nukleony v obou jadrech
  int impacts = nuc1->imp + nuc2->imp;

  //multiplicita
  double M_average = epsilon * (impacts * (1 - alpha) / 2 + alpha * nuc1->binImp);
  float M = generator->poisson(M_average);

  int NA = 0, NnA = 0, NB = 0, NnB = 0;

  //spektator nukleony a neutrony v nuc1
  for(int j = 0; j < nuc1->Z; j++){

    if(!nuc1->nucleons[j]->isImpact){

      ++NA;
      if(nuc1->nucleons[j]->isospin == -0.5){ ++NnA; }
    }
  }

  for(int j = 0; j < nuc2->Z; j++){

    if(!nuc2->nucleons[j]->isImpact){

      ++NB;
      if(nuc2->nucleons[j]->isospin == -0.5){ ++NnB; }
    }
  }

  nuc1->outputImp(b, impacts, M_average, M, NA, NnA, NB, NnB);  //zapis srazky do impacts.txt

	//testovací procedury

  //smallestR(nuc1);
	//nuc1->outputNucleons();
  nuc1->outputRad();

//smaz Nukleony
  delete nuc1;
  delete nuc2;

}


//priblizna doba vypoctu jedne srazky [mus]
float executionTime(string input1, int n1, string input2, int n2, float R){

  auto t1 = chrono::high_resolution_clock::now();
  collide(input1, n1, input2, n2, R, generator->genLinear());
  auto t2 = chrono::high_resolution_clock::now();

  auto duration1 = chrono::duration_cast<chrono::microseconds>( t2 - t1 ).count();

  t1 = chrono::high_resolution_clock::now();
  collide(input1, n1, input2, n2, R, generator->genLinear());
  t2 = chrono::high_resolution_clock::now();

  auto duration2 = chrono::duration_cast<chrono::microseconds>( t2 - t1 ).count();

  t1 = chrono::high_resolution_clock::now();
  collide(input1, n1, input2, n2, R, generator->genLinear());
  t2 = chrono::high_resolution_clock::now();

  auto duration3 = chrono::duration_cast<chrono::microseconds>( t2 - t1 ).count();

  return (duration1 + duration2 + duration3) / 3;

}

int main(){

  cout << "brum";
  //initializuj generator a konstanty
  konst = new Constants;
  generator = new Generator;

  ofstream cord;
  cord.open ("coordinates.txt");
  cord << " ";
  cord.close();

  cord.open ("rads.txt");
  cord << " ";
  cord.close();


  string input1;
  cout << "Enter symbol of the first element: ";
  cin >> input1;

  int n1;
  cout << "Enter nucleon number: ";
  cin >> n1;

  string input2;
  cout << "Enter symbol of the second element: ";
  cin >> input2;

  int n2;
  cout << "Enter nucleon number: ";
  cin >> n2;

  //ziskej polomer z ucinneho prurezu
  float s;

  cout << "Enter cross section [mb]" << endl;
  cout << "Sigma = ";
  cin >> s;

  s = s / 10; //mb -> fm^2

  float R = sqrt(s / M_PI);

  //vytvor hlavicku obsahujici maximalni mozny pocet srazenych Nukleonu a ucinny prurez v mb
  impactsFile.open("impacts.txt");
  impactsFile << "nuc = " << n1 + n2 << " sigma = " << 10*s << "mb " << n1 << input1 << " + " << n2 << input2 << endl;
  impactsFile.close();
  impactsFile.open("impacts.txt" , ios::app);

  int iter;
  cout << "Number of iterations: " << endl;
  cin >> iter;

  cout << "Processing..." << endl;


  //vypocti pribliznou dobu trvani srazeni
  cout << "Estimated time: " << iter * executionTime(input1, n1, input2, n2, R) / 1000000 << "s" << endl;
  cout << "---------------------------" << endl;

  //setina poctu iteraci
  int rat = round(iter / 100);

  for(int i = 1; i < iter; i++){

    //pokud je nynejsi iterace nasobek rat, vypis tento nasobek jakozto procento a zbyvajici cas
    if((i % rat) == 0){

      cout << (i / rat) << "%  ";
      cout << "Estimated time: " << (iter - i) * executionTime(input1, n1, input2, n2, R) / 1000000 << "s" << endl;
      continue;

    }

    //vytvor a sraz dve jadra
    collide(input1, n1, input2, n2, R, generator->genLinear());
  }

  cout << "Done. Results are in impacts.txt" << endl;
  cout << "First column - nucleons that hit at least one nucleon" << endl;
  cout << "Second column - number of all impacts" << endl;
  cout << "Third column - impact parameter" << endl;
  cout << "Fourth column - average multiplicity" << endl;
  cout << "Fifth column - multiplicity" << endl;
  cout << "Sixth column - spectator nucleons from nucleus A" << endl;
  cout << "Seventh column - spectator neutrons from nucleus A" << endl;
  cout << "Eighth column - spectator nucleons from nucleus B" << endl;
  cout << "Nineth column - spectator neutrons from nucleus B" << endl;

  impactsFile.close();
  delete konst;
  delete generator;

  return 0;

}
