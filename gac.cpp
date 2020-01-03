#include <iostream>
#include <fstream>
#include "centrality.h"
#include "sort.h"
#include "glauber.h"
using namespace std;

//tento program funguje jako jednotne uzivatelske rozhrani pro programy
//glauber.cpp, centrality.cpp a sort.cpp;
//sam o sobe nepridava zadnou funkcionalitu


//PRIKAZY
//-----------------------------------------------
void help(int lan){

  if(lan == 1){

    cout << " !help ... Vypiš seznam příkazů\n";
    cout << " !collide ... Spusť srážení jader\n";
    cout << " !prepare ... Zpracuj data z příkazu !collide pro následné použití v příkazu !centrality\n";
    cout << " !centrality ... Vyfiltruj z dat připravených příkazem !prepare pouze data obsahující multiplicitu v uživatelem daných hranicích\n";
    cout << " !exit ... Opusť program\n";
    cout << " !language ... Změn jazyk\n";

  }
  else{

    cout << " !help ... print list of commands\n";
    cout << " !collide ... Initialize the nuclear collision procedure\n";
    cout << " !prepare ... Process data from !collide command for later use in command !centrality\n";
    cout << " !centrality ... Filter from !prepare data only those that contain multiplicity in user defined interval\n";
    cout << " !exit ... Exit program\n";
    cout << " !language ... Change the language\n";

  }

}

void collide(int lan, string input){

  input = input.substr(input.find(" ") + 1);
  if(input == "-c")  glaub::start(lan, true, false);
  else if(input == "-r")  glaub::start(lan, false, true);
  else if((input == "-c -r") || (input == "-r -c")) glaub::start(lan, true, true);
  else glaub::start(lan, false, false);


}

void prepare(int lan){

  srt::start(lan);

}

void centrality(int lan){

  centr::centrality(lan);

}

//OSTATNI FUNKCE
//-----------------------------------------------
void setLanguage(int czech){

    ofstream lanFile("lan.txt");
    lanFile << czech;
    lanFile.close();

}

int getLanguage(){

  ifstream lanFile("lan.txt");
  string line;
  getline(lanFile, line);
  lanFile.close();

  try{ int input = stoi(line); return input; }
  catch(const invalid_argument& ia){  return -1;  }

}


void commands(int lan){

  cout << "> ";
  string input;
  getline(cin, input);

  bool e = false;

  if(input == "!help") help(lan);
  else if(input.substr(0, 8) == "!collide") collide(lan, input);
  else if(input == "!prepare") prepare(lan);
  else if(input == "!centrality") centrality(lan);
  else if(input == "!exit")  e = true;
  else if(input == "!language"){ lan = abs(lan - 1); setLanguage(lan); }
  else  cout << "Not valid command\n";

  if(!e){
    commands(lan);
  }
}


void first(int lan){

  if(lan == 1){

    cout << "Vítejte v programu GAC, počítačové simulaci srážek těžkých jader\n";
    cout << "Pro seznam příkazů zadejte !help\n";
    commands(1);

  }
  else{

    cout << "Welcome in GAC, program simulating heavy nuclear collisions\n";
    cout << "For the command list enter !help\n";
    commands(0);

  }
}


int main(){

    int language = getLanguage();
    if(language == -1){

      cout << "Which language would you like to use? (english [E] or czech [C]) \n/ Jaký jazyk chcete používat? (angličtina [A], nebo čeština [C])\n";
      cout << "> ";

      char language;
      cin >> language;

      if((language == 'C') || (language == 'c')){ setLanguage(1); first(1); }
      else{ setLanguage(0); first(0); }

    }
  else if(language == 1)  first(1);
  else  first(0);
  return 0;
}
