# GAC
Program pro simulování atomových srážek a definování tříd centrality podle Glauberova modelu. 
Program je součástí práce SOČ "Simulace ultra-relativistických srážek atomových jader podle Glauberova modelu". 

gac(Linux).exe a gac(Windows).exe jsou zkompilované spustitelné verze programu.
gac.cpp tvoří rámec programu, lze z něj volat ostatní části
glauber.cpp a glauber.h - simulace srážení jader
sort.cpp a sort.h - příprava dat pro definování centrality
centrality.cpp a centrality.h - filtrování dat odpovídající třídě centrality

Pro vlastní kompilaci je potřeba mít staženou knihovnu GNU Scientific Library.
