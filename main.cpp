#include <iostream>
#include <iomanip>
#include <algorithm>
#include <math.h>
const double EPS = 10e-11;

int Menu();
/*
s - rozmiar tablicy, na ktorej dzialamy
i - krok operacji, dla ktorej liczymy
tab/results - tablice z danymi/wnikami
row/col - ref. od kolumn/wierszy do metody z max ele. w wierszu/kolumnie
*/

void GetMatrix(int s, double ** tab);
void PrintfMatrix(int s, double ** tab);
void PrintfResults(int s, double * t);
bool ReductionX(int i, int s, double ** tab);
bool CalX(int s, double ** tab, double * results);
bool MaxElementMatrix(int s, int kr, int & w, int & k, double ** tab);
bool MaxElementRow(int s, int w, int & kol, double ** tab);
bool MaxElementColumn(int s, int k, int & wer, double ** tab);
void ReplaceRow(int s, int & row, int i, double ** tab);
void ReplaceCol(int s, int & col, int i, double ** tab, double * results);
void Gauss(int s, double ** tab, double * results);
bool GaussFull(int s, double ** tab, double * results);
void GaussRow(int s, double ** tab, double * results);
void GaussCol(int s, double ** tab, double * results);


int main()
{
    //Mozna by tu dodac do-while do powatrzania programu by przetestowac metody wszystkie;
    //ew dac opcje wczytania sztywnych liczb dla przetestowania lub wczytywac je z pliku
    //ale to chyba nie o to chodzi;
    int n, w;
    bool t;
    double test1 = 0.0000000001;
    double test2 = 0.0000000001;
    std::cout << "1)" << test1 << " a " << test2 << "Podaj liczbe rownan: ";
    test1 = 10e-11;
    test2 = 1e-010;
    std::cout << "1)" << test1 << " a " << test2 << "Podaj liczbe rownan: ";

    std::cin >> n;
    double ** tab = new double * [n];
    double * results = new double[n];
    for(int i = 0; i < n; i++) tab[i] = new double[n + 1];

    GetMatrix(n, tab);
    std::cout << std::endl;
    PrintfMatrix(n, tab);

    do
    {
        w = Menu();
        switch(w)
        {
            case 0: break;
            case 1:
            {
                Gauss(n, tab, results);
                PrintfResults(n, results);
                break;
            }
            case 2:
            {
                t = GaussFull(n, tab, results);
                if(t == false)
                {
                    std::cout << "Obliczenie ukladu jest niemozliwe";
                    break;
                }
                else PrintfResults(n, results);
                break;
            }
            case 3:
            {
                std::cout << "Metoda z wyborem elementu maksymalnego w wierszu:";
                GaussRow(n, tab, results);
                PrintfResults(n, results);
                break;
            }
            case 4:
            {
                std::cout << "Metoda z wyborem elementu maksymalnego w kolumnie";
                GaussCol(n, tab, results);
                PrintfResults(n, results);
                break;
            }
        }
    } while(w < 0 && w > 3);

    for(int i = 0; i < n; i++) delete[] tab[i];
    delete[] tab;
    tab = NULL;

    delete[] results;
    results = NULL;
    return 0;
}

int Menu()
{
    int x;
    do
    {
        std::cout << "\tMENU:\n";
        std::cout << "1 - Metoda podstawowa \n2 - Metoda z pelnym wyborem elementu max\n"
              << "3 - Metoda z wyborem max. elementu w wierszu\n4 - Metoda z wyborem max. elementu w kolumnie\n";
        std::cin >> x;
    } while(!(x >= 0 && x <= 4));
    return x;
}

void GetMatrix(int s, double ** tab)
{
    for(int i = 0; i < s; i++)
    {
        for(int j = 0; j < s + 1; j++)
        {
            if(j != s)
            {
                std::cout << "A[" << i << ", " << j << "'] = ";
                std::cin >> tab[i][j];
            }
            else
            {
                 std::cout << "wyraz wolny: A[" << i << ", " << j << "] = ";
                 std::cin >> tab[i][j];
            }
        }
    }
}

void PrintfMatrix(int s, double ** tab)
{
    std::cout << "Uklad rownan:\n";
    for(int i = 0; i < s; i++)
    {
        for(int j = 0; j < s + 1; j++)
        {
            if(j != s)
            {
                std::cout << std::setw(5);
                std::cout << tab[i][j];

            }
            else
            {
                std::cout << " |";
                std::cout << std::setw(4);
                std::cout << tab[i][j];
            }
        }
        std::cout << "\n";
    }
}

void PrintfResults(int s, double * t)
{
    std::cout << "\nX = [ ";

    for(int i = 0; i < s - 1; i++) std::cout << std::setw(3) << t[i] << "; ";
    std::cout << std::setw(3);
    std::cout << t[s - 1] << " ]";
}

bool ReductionX(int i, int s, double ** tab)
{
    double z = 0;
    for(int j = i + 1; j < s; j++)
    {
        if(fabs(tab[i][i]) < EPS)
        {
            return false;
        }

        z = tab[j][i] / tab[i][i];
        for(int a = 0; a < s + 1; a++) tab[j][a] = tab[j][a] - (z * tab[i][a]);
    }
    return true;
}

bool CalX(int s, double ** tab, double * results)
{
    int j;
    bool t = true;
    double w, vec[s], v[s];
    std::copy(results, results + s, vec);
    for(int i = s - 1 ; i > -1; --i)
    {
        if(fabs(tab[i][i]) < EPS)
        {
            return false;
        }
        w = tab[i][s];
        for(int j = i + 1; j < s; j++) w -= (tab[i][j] * vec[j]);
        vec[i] = w / tab[i][i];
    }

    for(int i = 0; i < s; i++)
    {
        j = 0;
        while(j != results[i]) j += 1;
        v[j] = vec[i];
    }
    std::copy(v, v + s, results);
    return true;
}

void Gauss(int s, double ** tab, double * results)
{
    std::cout << "Podstawowa metoda eliminacji Gaussa\n";
    for(int i = 0; i < s; i++) results[i] = i;
    for(int j = 0; j < s - 1; j++) ReductionX(j, s, tab);
    CalX(s, tab, results);
}

bool MaxElementMatrix(int s, int kr, int & w, int & k, double ** tab)
{
    double max_;
    max_ = fabs(tab[kr][kr]);
    w = kr;
    k = kr;
    for(int i = kr; i < s; i++)
    {
        for(int j = kr; j < s; j++)
        {
            if(max_ < fabs(tab[i][j]))
            {
                max_ = tab[i][j];
                w = i;
                k = j;
            }
        }
    }
    if(fabs(max_) < EPS) return false;
    else return true;
}

void ReplaceRow(int s, int & row, int i, double ** tab)
{
    for(int j = i; j < s + 1; j++) std::swap(tab[i][j], tab[row][j]);
}


void ReplaceCol(int s, int & col, int i, double ** tab, double * results)
{
    for(int j = 0; j < s; j++) std::swap(tab[j][i], tab[j][col]);

    std::swap(results[i], results[col]);
}

bool GaussFull(int s, double ** tab, double * results)
{
    std::cout << "Metoda eliminacji Gaussa z wyborem max. elementu\n";
    int row = 0, column = 0;
    bool t = true;
    for(int i = 0; i < s; i++) results[i] = i;

    for(int j = 0; j < s - 1; j++)
    {
        t = MaxElementMatrix(s, j, row, column, tab);
        if(t == false)
        {
            return false;
            break;
        }

        if(row != j) ReplaceRow(s, row, j, tab);

        MaxElementMatrix(s, j, row, column, tab);

        if(column != j) ReplaceCol(s, column, j, tab, results);

        t = t & ReductionX(j, s, tab);
        if(t == false)
        {
            return false;
            break;
        }

    }
    t = t & CalX(s, tab, results);
    return t;
}

bool MaxElementRow(int s, int w, int & kol, double ** tab)
{
    double max_;
    max_ = tab[w][w];
    kol = w;
    for(int i = w; i < s; i++)
    {
        if( max_ < fabs(tab[w][i]))
        {
            max_ = tab[w][i];
            kol = i;
        }
    }
    if(fabs(max_) < EPS) return false;
    return true;
}

void GaussRow(int s, double ** tab, double * results)
{
    bool t;
    int column = 0;
    for(int i = 0; i < s; i++) results[i] = i;

    for(int j = 0; j < s - 1; j++)
    {
        t = MaxElementRow(s, j, column, tab);
        if(t == false) break;
        if(column != j) ReplaceCol(s, column, j, tab, results);
        ReductionX(j, s, tab);
    }
    CalX(s, tab, results);
}

bool MaxElementColumn(int s, int k, int & wer, double ** tab)
{
    double max_;
    max_ = tab[k][k];
    wer = k;
    for(int i = k; i < s; i++)
    {
        if( max_ < fabs(tab[i][k]))
        {
            max_ = tab[i][k];
            wer = i;
        }
    }
    if(fabs(max_) < EPS) return false;
    return true;
}

void GaussCol(int s, double ** tab, double * results)
{
    bool t;
    int row  = 0;
    for(int i = 0; i < s; i++) results[i] = i;

    for(int j = 0; j < s - 1; j++)
    {
        t = MaxElementColumn(s, j, row, tab);
        if(t == false) break;
        if(row != j) ReplaceRow(s, row, j, tab);
        ReductionX(j, s, tab);
    }
    CalX(s, tab, results);
}

