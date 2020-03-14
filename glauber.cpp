#include "glauber.h"
using namespace std;

//----------------
//GLOBAL VARIABLES
//----------------
Generator * generator;
Constants * konst;
ofstream impactsFile;
bool glob_returnRads;
bool glob_returnCoords;

//---------------------
//FUNCTION DECLARATIONS
//---------------------
void glaub::smallestR(Nucleus * nuc);

//---------------
//CLASS FUNCTIONS
//---------------

//Class UI
//---------------------------------------------------------------------------------------------
UI::UI(int language){

  if(language == 1) this->english = false;
  else  this->english = true;

}

void UI::englishInput(){

  cout << "Enter cross section [mb]" << endl;
  cout << "Sigma = ";
  cin >> this->sigma;

  cout << "Parameter for the calculation of multiplicity alpha: ";
  cin >> this->alpha;

  cout << "Number of iterations: " << endl;
  cin >> this->iter;

  cout << "Processing..." << endl;

}
void UI::czechInput(){

  cout << "Účinný průřez [mb]" << endl;
  cout << "Sigma = ";
  cin >> this->sigma;

  cout << "Parametr pro výpočet multiplicity alpha: ";
  cin >> this->alpha;

  cout << "Počet opakování: " << endl;
  cin >> this->iter;

  cout << "Počítám..." << endl;

}

void UI::englishTime(float time){

  cout << "Estimated execution time: " << time << "s" << endl;
  cout << "____________________________________" << endl;

}

void UI::czechTime(float time){

  cout << "Předpokládaná doba výpočtu: " << time << "s" << endl;
  cout << "____________________________________" << endl;

}

void UI::englishPercent(int percent, float time){

  cout << percent << "%  ";
  cout << "Estimated time: " << time << "s" << endl;

}

void UI::czechPercent(int percent, float time){

  cout << percent << "%  ";
  cout << "Předpokládaný zbývající čas: " << time << "s" << endl;

}

void UI::englishOutput(){

  cout << "Done. Results are in impacts.txt" << endl;
  cout << "First column - nucleons that hit at least one nucleon" << endl;
  cout << "Second column - number of all impacts" << endl;
  cout << "Third column - impact parameter" << endl;
  cout << "Fourth column - average multiplicity" << endl;
  cout << "Fifth column - multiplicity" << endl;
  cout << "Sixth column - spectator nucleons from the first nucleus" << endl;
  cout << "Seventh column - spectator neutrons from the first nucleus" << endl;
  cout << "Eighth column - spectator nucleons from the second nucleus" << endl;
  cout << "Nineth column - spectator neutrons from the second nucleus" << endl;
  cout << "Tenth column - excentricity" << endl;

}
void UI::czechOutput(){

  cout << "Hotovo. Výsledky jsou v impacts.txt" << endl;
  cout << "První sloupec - Počet sražených nukleonů" << endl;
  cout << "Druhý sloupec - Počet všech srážek" << endl;
  cout << "Třetí sloupec - Parametr srážky [fm]" << endl;
  cout << "Čtvrtý sloupec - Průměrná multiplicita" << endl;
  cout << "Pátý sloupec - Multiplicita" << endl;
  cout << "Šestý sloupec - Spektátoři (nukleony) z prvního jádra" << endl;
  cout << "Sedmý sloupec - Spektátoři (neutrony) z prvního jádra" << endl;
  cout << "Osmý sloupec - Spektátoři (nukleony) z druhého jádra" << endl;
  cout << "Devátý sloupec - Spektátoři (neutrony) z druhého jádra" << endl;
  cout << "Desátý sloupec - Excentricita" << endl;

}

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
double Generator::gen(){

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
double Generator::genLinear(){

  return konst->maxB * sqrt(this->gen());

}

//generuj z rovnomerneho rozdeleni v intervalu <0.5 ; 13.5> (v intervalu velikosti jadra)
double Generator::genInNuc(){

    return 2*(konst->nucleusR - konst->nucleonR) * this->gen() + konst->nucleonR;

  }

//ziskej multiplicitu z Poissonova rozdeleni
int Generator::poisson(float M_average){

  return gsl_ran_poisson(this->generator, M_average);

}

double Generator::genNormWoodSaxon(float a, float R0){

  double trilogarithm = -gsl_sf_fermi_dirac_2(R0 / a);
  double N = -1 / (2 * pow(a, 3) * trilogarithm);

  //rejection method
  while(true){

    double x = this->genPosition(0, konst->maxR);
    double validator = this->gen();

    double y = N * x * x / (1 + exp( (x - R0) / a ));
    double max = 0.5; //snad se nepresahne... najit maximum funkce vyzaduje Lambertovu W funkci
    if((y / max) > validator){ return x; } //vhodna hodnota
  }
}

//ziskej vzdalenost od centra jadra z Wood-Saxonova potencialu
double Generator::genWoodSaxon(float a, float R0){

  //normalizacni konstanta
  double N = 1 / (a * log( exp(R0 / a) + 1) );

  double temp = exp(-this->gen() / (N * a));
  double kon = exp( -R0 / a);

  while(temp <= kon){ temp = exp(-this->gen() / (N * a)); }

  double y = -a * log( temp - kon );
  //while(isnan(y)){ y = this->genWoodSaxon(a, R0); }
  return y;

}

//z rozdeleni sin(x)/2
double Generator::genSinus(){

  return acos(1 - 2*this->gen());

}


double Generator::genPosition(double min, double max){

  return ( this->gen() * (max - min) + min );

}


//Class Map
//---------------------------------------------------------------------------------------------

  //konstruktor Mapy
  Map::Map(){

    int size = 2*ceil(konst->nucleusR) + 1;
    this->map = vector< vector< vector< vector<double> > > > (size,vector<vector<vector<double>>>(size,vector<vector<double>>(size,vector<double>(0))));

  }

  //zapis do Mapy
  void Map::writeCoords(double x, double y, double z){

    //nejblizsi uzel
    int roundX = round(x);
    int roundY = round(y);
    int roundZ = round(z);

    if(roundX <= 0){ roundX = 1; } if(roundX >= 2*konst->nucleusR){ roundX = 2*ceil(konst->nucleusR) - 1; }
    if(roundY <= 0){ roundY = 1; } if(roundY >= 2*konst->nucleusR){ roundY = 2*ceil(konst->nucleusR) - 1; }
    if(roundZ <= 0){ roundZ = 1; } if(roundZ >= 2*konst->nucleusR){ roundZ = 2*ceil(konst->nucleusR) - 1; }

    this->map[roundX][roundY][roundZ].push_back(x);
    this->map[roundX][roundY][roundZ].push_back(y);
    this->map[roundX][roundY][roundZ].push_back(z);

  }

  //ziskej souradnice z Mapy
  vector<double> Map::getCoords(int x, int y, int z){

    return this->map[x][y][z];
  }

  //smaz souradnice z Mapy
  void Map::deleteCoords(int x, int y, int z, int start, int length){

    this->map[x][y][z].erase(this->map[x][y][z].begin() + start, this->map[x][y][z].begin() + start + length);

  }

  //zjisti, jestli Nukleon na pozici x, y, z koliduje s jinym nukleonem
  bool Map::isColliding(double x, double y, double z){

    int roundX = round(x);
    int roundY = round(y);
    int roundZ = round(z);

    if(roundX <= 0){ roundX = 1; } if(roundX >= 2*konst->nucleusR){ roundX = 2*ceil(konst->nucleusR) - 1; }
    if(roundY <= 0){ roundY = 1; } if(roundY >= 2*konst->nucleusR){ roundY = 2*ceil(konst->nucleusR) - 1; }
    if(roundZ <= 0){ roundZ = 1; } if(roundZ >= 2*konst->nucleusR){ roundZ = 2*ceil(konst->nucleusR) - 1; }

    //projdi vsechny okolni uzly
    for(int i = -1; i < 2; ++i){  //x
      for(int j = -1; j < 2; ++j){  //y
        for(int k = -1; k < 2; ++k){  //z

          //vsechny souradnice nukleonu v uzlu
          vector<double> coords = this->getCoords(roundX + i, roundY + j, roundZ + k);

          int size = coords.size();
          if(size % 3 != 0){ cout << konst->exceptions[3] << endl; }
          int n = size / 3;  //souradnice jsou vzdy 3 -> pocet nukleonu je delen tremi

          for(int u = 0; u < n; ++u){

            double nearX = coords[3*u];
            double nearY = coords[3*u+1];
            double nearZ = coords[3*u+2];

            double r = sqrt( pow(x - nearX, 2) + pow(y - nearY, 2) + pow(z - nearZ, 2) );

            if(r < 1){ return true; }

          }
        }
      }
    }

    return false;
  }

//Class Nucleon
//---------------------------------------------------------------------------------------------

  //konstruktor Nukleonu
  Nucleon::Nucleon(float I, Nucleus * par, int num){

    this->index = num;
    this->parent = par;
    this->isospin = I;

    do{
    vector<double> coords = this->makeSphericalCoords();

    double r = coords[0];
    double phi = coords[1];
    double theta = coords[2];

    //nahodne pricti nebo odecti od nucleusR (stred neni v pocatku souradnic)
    float ranX = generator->gen(); float ranY = generator->gen(); float ranZ = generator->gen();
    int unitX, unitY, unitZ;
    if(ranX > 0.5){ unitX = 1;  }else{  unitX = -1; }
    if(ranY > 0.5){ unitY = 1;  }else{  unitY = -1; }
    if(ranZ > 0.5){ unitZ = 1;  }else{  unitZ = -1; }

    //transformace
    this->x = konst->nucleusR + unitX * r * sin(theta) * cos(phi);
    this->y = konst->nucleusR + unitY * r * sin(theta) * sin(phi);
    this->z = konst->nucleusR + unitZ * r * cos(theta);

    }//opakuj, pokud potreba
    while(this->parent->map->isColliding(this->x, this->y, this->z));

  }


  vector<double> Nucleon::makeSphericalCoords(){

    vector<double> coords(3);

    double r = generator->genNormWoodSaxon(konst->a, konst->nucleusR);
    double phi = generator->genPosition(0, 2*M_PI);
    double theta = generator->genSinus();

    coords[0] = r;  coords[1] = phi;  coords[2] = theta;

    return coords;

  }


//Class Nucleus
//---------------------------------------------------------------------------------------------

//destruktor
Nucleus::~Nucleus(){

  for(int i = 0; i < this->filled; i++){ delete this->nucleons[i]; }
  delete this->map;

}

//konstruktor
Nucleus::Nucleus(){

    this->map = new Map();

    this->createNucleons();

}

int Nucleus::getBin(float x, float y){

  int bin = 0;

  for(int i = 0; i < this->filled; ++i){

    Nucleon * nucleon = this->nucleons[i];
    float r = sqrt( pow(nucleon->x - x, 2) + pow(nucleon->y - y, 2) );

    if(r < konst->nucleonR){  bin += nucleon->binImp; }

  }

  return bin;
}

float Nucleus::excentricity(){

  float lowBound = konst->nucleusR - konst->maxR;
  float highBound = konst->nucleusR + konst->maxR;
  const int iterations = 50;

  auto funcX = [this](float x, float y)->float{ return (float) this->getBin(x, y) * x * x;  };
  auto funcY = [this](float x, float y)->float{ return (float) this->getBin(x, y) * y * y;  };

  float integral_x = glaub::integral(funcX, lowBound, highBound, iterations);
  float integral_y = glaub::integral(funcY, lowBound, highBound, iterations);

  return (integral_y - integral_x) / (integral_y + integral_x);

}

//zapis vzdalenosti Nukleonu od jadra do rads.txt (POUZE TESTOVACI)
void Nucleus::outputRad(){

  ofstream file;
  file.open ("rads.txt", ios::app);

  float r;

  for(int i = 0; i < konst->A; i++){

      r = sqrt(pow(this->nucleons[i]->x - konst->nucleusR, 2) + pow(this->nucleons[i]->y - konst->nucleusR, 2) + pow(this->nucleons[i]->z - konst->nucleusR, 2));
      file << r << "\n";

  }

  file.close();

}

//zapis souradnice Nukleonu do coordinates.txt (POUZE TESTOVACI)
void Nucleus::outputNucleons(){

  ofstream file;
  file.open ("coordinates.txt", ios::app);

  for(int i = 0; i < konst->A; i++){

    file << this->nucleons[i]->x << " ";
    file << this->nucleons[i]->y << " ";
    file << this->nucleons[i]->z << "\n";

  }

  file.close();

}

//vypis pocet srazenych Nukleonu, pocet binarnich srazek a srazkovy parametr b do impacts.txt
void Nucleus::outputImp(float b, float impacts, double M_average, int M, int NA, int NnA, int NB, int NnB, float excentricity){


  impactsFile << impacts << " " << this->binImp << " " << b << " " << M_average << " " << M << " " << NA << " " << NnA<< " " << NB << " " << NnB << " " << excentricity << "\n";

}

//vytvor Nukleony do jadra
void Nucleus::createNucleons(){

  for(int i = 0; i < konst->A; i++){

    if(i < konst->Z){
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
template<typename T>
float glaub::integral(T&& lambda, float low, float high, int iterations){

  float x = low;
  float y = low;
  float jump = (high - low) / (float) iterations;
  float output = 0;

  for(int i = 0; i < iterations; ++i){

    x = low;
    for(int j = 0; j < iterations; ++j){

        output += jump * lambda(x, y);
        x += jump;

    }

    y += jump;
  }

  return output;
}

//projdi vzdalenosti vsech Nukleonu a vypis tu nejmensi a nejvetsi (POUZE TESTOVACI)
void glaub::smallestR(Nucleus * n){

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

  cout << "Minimal distance: " << min << "\n";
  cout << "Maximal distance: " << max << "\n";

}

//vytvor a sraz dve jadra (p1 a p2 znacky prvku, n1 a n2 nukleonova cisla, R polomer srazky ziskany z ucinneho prurezu, b srazkovy parametr)
bool glaub::collide(float R, float alpha){

  float b = generator->genLinear();

  Nucleus * nuc1 = new Nucleus;
  Nucleus * nuc2 = new Nucleus;

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
        if(!nuc1->nucleons[i]->isImpact){++nuc1->imp;}
        if(!nuc2->nucleons[j]->isImpact){++nuc2->imp;}

        nuc1->nucleons[i]->isImpact = true;
        nuc2->nucleons[j]->isImpact = true;
        ++nuc1->nucleons[i]->binImp;
        ++nuc1->binImp; //pricti jadru binarni srazku
      }
    }
  }


  //trefene nukleony v obou jadrech
  int impacts = nuc1->imp + nuc2->imp;

  //multiplicita
  double M_average = impacts * (1 - alpha) / 2 + alpha * nuc1->binImp;
  int M = generator->poisson(M_average);

  if(M == 0){

    delete nuc1;
    delete nuc2;
    return false; //minimum bias

  }

  int NA = 0, NnA = 0, NB = 0, NnB = 0;

  //spektator nukleony a neutrony v nuc1
  for(int j = 0; j < konst->A; j++){

    if(!nuc1->nucleons[j]->isImpact){

      ++NA;
      if(nuc1->nucleons[j]->isospin == -0.5){ ++NnA; }
    }
  }

  for(int j = 0; j < konst->A; j++){

    if(!nuc2->nucleons[j]->isImpact){

      ++NB;
      if(nuc2->nucleons[j]->isospin == -0.5){ ++NnB; }
    }
  }

  nuc1->outputImp(b, impacts, M_average, M, NA, NnA, NB, NnB, nuc1->excentricity());  //zapis srazky do impacts.txt

	//testovací procedury
  //smallestR(nuc1);
	if(glob_returnCoords)  nuc1->outputNucleons();
  if(glob_returnRads) nuc1->outputRad();

//smaz Nukleony
  delete nuc1;
  delete nuc2;

  return true;
}


//priblizna doba vypoctu jedne srazky [mus]
float glaub::executionTime(float R, float alpha){

  auto t1 = chrono::high_resolution_clock::now();
  collide(R, alpha);
  collide(R, alpha);
  collide(R, alpha);
  auto t2 = chrono::high_resolution_clock::now();

  auto duration1 = chrono::duration_cast<chrono::microseconds>( t2 - t1 ).count();



  return duration1 / 3;

}

void glaub::start(int language, bool returnCoords, bool returnRads){

  glob_returnRads = returnRads;
  glob_returnCoords = returnCoords;

  //initializuj generator a konstanty
  konst = new Constants;
  generator = new Generator;
  UI * ui = new UI(language);

  //smaz minule zaznamy
  ofstream cord;
  cord.open ("coordinates.txt");
  cord << " ";
  cord.close();

  cord.open ("rads.txt");
  cord << " ";
  cord.close();

  //ziskej vstupni hodnoty
  if(ui->english){ ui->englishInput(); }else{ ui->czechInput(); }

  ui->sigma = ui->sigma / 10; //mb -> fm^2

  //polomer nutny pro srazku z ucinneho prurezu
  float R = sqrt(ui->sigma / M_PI);

  //vytvor hlavicku obsahujici maximalni mozny pocet srazenych Nukleonu a ucinny prurez v mb
  impactsFile.open("impacts.txt");
  impactsFile << "nuc = " << 2*konst->A << " sigma = " << 10*ui->sigma << "mb " << " alpha = " << ui->alpha << " " << konst->A << "Au + " << konst->A << "Au" << endl;
  impactsFile.close();
  impactsFile.open("impacts.txt" , ios::app);


  //vypocti pribliznou dobu trvani srazeni
  if(ui->english){ ui->englishTime(ui->iter * executionTime(R, ui->alpha) / 1000000); }
  else{ ui->czechTime(ui->iter * executionTime(R, ui->alpha) / 1000000); }


  //setina poctu iteraci
  int rat = round(ui->iter / 100);
  if(rat == 0){ ++rat; }

  for(int i = 1; i < ui->iter; i++){

    //pokud je nynejsi iterace nasobek rat, vypis tento nasobek jakozto procento a zbyvajici cas
    if((i % rat) == 0){

      if(ui->english){ ui->englishPercent(i / rat, (ui->iter - i) * executionTime(R, ui->alpha) / 1000000); }
      else{ ui->czechPercent(i / rat, (ui->iter - i) * executionTime(R, ui->alpha) / 1000000); }
      continue;

    }

    //vytvor a sraz dve jadra
    if( !collide(R, ui->alpha) )  --i;  //pokud multiplicita nulova, nezapocitej

  }

  //popisky sloupcu
  if(ui->english){ ui->englishOutput(); }else{ ui->czechOutput(); }

  impactsFile.close();
  delete konst;
  delete generator;


}
