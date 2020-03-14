# GAC
Program pro simulování atomových srážek a definování tříd centrality podle Glauberova modelu. 
Program je součástí práce SOČ "Simulace ultra-relativistických srážek atomových jader podle Glauberova modelu". 

gac(Linux).exe a gac(Windows).exe jsou zkompilované spustitelné verze programu.<br/>
gac.cpp tvoří rámec programu, lze z něj volat ostatní části<br/>
glauber.cpp a glauber.h - simulace srážení jader<br/>
sort.cpp a sort.h - příprava dat pro definování centrality<br/>
centrality.cpp a centrality.h - filtrování dat odpovídající třídě centrality<br/><br/>
Pro vlastní kompilaci je potřeba mít staženou knihovnu GNU Scientific Library.<br/>

### Kompilace
g++ -I<CESTA KE SLOŽCE>/gsl/include -c gac.cpp<br/>
g++ -c centrality.cpp<br/>
g++ -c sort.cpp<br/>
g++ -Wall -I<CESTA KE SLOŽCE>/gsl/include -c glauber.cpp<br/>
g++ gac.o centrality.cpp sort.cpp glauber.o -L<CESTA KE SLOŽCE>/gsl/lib -lgsl -lgslcblas -lm -o gac.exe
