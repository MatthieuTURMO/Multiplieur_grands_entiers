#include <iostream>
#include <string.h>
#include <math.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <complex>
#include <complex.h>

using namespace std;

int calculM(string a, string b);
int inv(int j, int m);
vector<complex<double> > decouper(string a, string b);
vector<complex<double> > cooleyTukey(vector<complex<double> > vect, int m, int inverse);
vector<complex<double> > multipliePolynome(vector<complex<double> > vectA, vector<complex<double> > vectB);
void ssa(string A, string B);
complex<double> nettoyer(complex<double> A);
int arrondi(double a);


int main()
{
    
    cout << "Algorithme de Schonhage-Strassen"<<endl;
    cout << "Realise par Matthieu TURMO et Nicolas Wacquant"<<endl<<endl;
    string sa, sb;
    cout<<"Entrez le premier nombre entier a multiplier"<<endl;
    cin >> sa;
    cout << "Entrez maintenant le second nombre"<<endl;
    cin >> sb;
    
    ssa(sa,sb);
    
    
}


//methode Stronhage-Strassen
//affichera le reultat final
void ssa(string A, string B){
    
    vector<complex<double> > vectA;
    vector<complex<double> > vectB;
    vector<complex<double> > resultat;
    vector<complex<double> > vraiResultat;
    int tmp;
    int cpt;
    stringstream ss;
    string resultString = "";
    string resultatFinal = "";
    string nombreFinal ="";
    
    
    //on decoupe les deux nombre passes en parametres
    //et on les stocke dans deux vector de complexe (de type double)
    vectA = decouper(A, B);
    vectB = decouper(B,A);
    
    //m sera la puissance de deux directement superieur
    //a la somme des ordres des deux nombres
    int m = calculM(A,B);
    
    //pour chaque nombre A et B, on calcule leur transformee de fourier
    vectA = cooleyTukey(vectA, m, 1);
    vectB = cooleyTukey(vectB, m, 1);
    
    //puis on  multiplie le polynome resultant
    resultat = multipliePolynome(vectA, vectB);
    
    
    
    int j;
    
    //vraiResultat est le resultat de la multiplication des deux polynomes
    //auquel on applique la methode inv, qui inverse les indices
    for(int i=0; i< resultat.size(); i++){
        j = inv(i,m);//on applique la methode inv a l'indice i
        //et on met la valeur correspondante dans le vector
        vraiResultat.push_back(resultat.at(j));
    }
    
    
    //puis on lui applique la transformee de fourier inverse
    vraiResultat = cooleyTukey(vraiResultat, m, -1);
    
    //on divise ensuite tous les nombres contenus dans ce vector par N = 2^m
    for(int i=0; i<vraiResultat.size(); i++){
        
        vraiResultat.at(i)/=pow(2,m);
    }
    
    //on initialise maintenant un tableau d'entier,
    //contenant toutes les valeurs de vraiResultat, mais en valeur entiere
    int tab[vraiResultat.size()];
    
    for(int i=0; i<vraiResultat.size(); i++){
        //on appelle la methode arrondi(), pour que les valeurs soient bien conservees
        //sinon, etant donne que vraiResultat contient des doubles, il y aura un chiffre de difference
        tab[i] = arrondi(real(vraiResultat.at(i)));
        
    }
    
    //la boucle ci-dessous permet de ne laisser dans le vector que des chiffres
    //le tableau tab contiendra donc maintenant les chiffres du resultat total
    for(int i=0; i<vraiResultat.size()-1; i++){
        
        cpt = 0;
        tmp = tab[i];
        if(tmp >= 10){
            //si la valeur de tab[i] est superieure a 10
            while(tmp >= 10){ //on boucle
                //on lui retranche 10 a chaque tour
                //cpt contient le nombre de fois que l'on a rentranche 10 au nombre
                cpt++;
                tmp = tmp - 10;
            }
        }
        
        //on met tmp dans tab[i], si tmp n'est pas superieur a 10, la valeur ne change pas
        tab[i] = tmp;
        //on ajoute a la case suivante le nombre de dizaines qui etait dans tab[i]
        tab[i+1] += cpt;
    }
    
    //on met dans le flux ss, tous les chiffres du resultat, dans le bon ordre
    for(int i=vraiResultat.size()-1; i>=0; i--){
        ss << tab[i];
    }
    //on met le contenu du flux dans la chaine de caracetere resultString
    ss >> resultString;
    
    
    //resultString contient maintenant le resultat, mais on va enlever les 0
    //pour plus de clarte
    
    //on recupere le cpt qui correcspond au premier entier non nul
    for(int i=0; i<resultString.size(); i++){
        if(resultString.at(i)!='0'){
            //cpt contient maintenant la valeur de la premiere occurence
            //d'un nombre different de 0
            cpt = i;
            break;
        }
    }
    
    
    
    //resultatFinal va maintenant contenir le nombre, sans les 0 inutiles
    for(int i = cpt; i<resultString.size(); i++){
        resultatFinal += resultString.at(i);
    }
    
    //on remet le compteur a 0
    cpt = 0;
    
    //la methode ci-dessous va maintenant decouper le nombre en bloc de 3 chiffres
    //en partant de la fin du nombre
    //toujours pour plus de lisibilite
    for(int i=resultatFinal.size()-1; i>=0; i--){
        
        //si cpt est un multiple de 4
        if(cpt%4==0){
            //on affiche un espace, et on incremente l'indice i,
            //afin de ne pas oublier de chiffre lors de l'affichage
            nombreFinal+=" ";
            i++;
        }
        else //sinon on stocke tout simplement le chiffre
            nombreFinal+=resultatFinal.at(i);
        
        //et on incremente le compteur
        cpt++;
        
    }
    
    cout<<endl<<"Voici le resultat de la multiplication :"<<endl;
    //finalement, on affiche le resultat, ce resultat est le nombre final de la multiplication.
    for(int i=nombreFinal.size()-1; i>=0; i--){
        cout<<nombreFinal.at(i);
    }
    cout<<endl;
    
}//FIN METHODE SSA

//cette methode retourne un entier a partie d'un double
//l'entier retourne sera l'arrondi du double passe en parametre
int arrondi(double a){
    //on stocke la partie entiere du nombre
    int partieEntiere = floor(a);
    //si la difference est superieure a 0.2, on incremente la valeur retournee
    if(a-partieEntiere >= 0.2){
        return partieEntiere+1;
    }
    else    //sinon on retourne simplement la valeur de la partie entiere
        return partieEntiere;
}


//cette methode servira dans l'algorithme de CooleyTukey
//etant donne l'encodage des nombres, si une petite valeur subsiste
//telle que 0.000000000001, pour ne pas fausser les calculs, on la mettra a 0
//on le fera pour la partie relle et pour la partie imaginaire du complexe passe en parametre
complex<double> nettoyer(complex<double> A){
    complex<double> res;
    if(abs(real(A)) < pow(10,-3)){
        res.real(0);
       // real(res) = 0;
    }
    else{
        res.real(real(A));
        //real(res) = real(A);
    }
    if(abs(imag(A)) < pow(10,-3)){
        res.imag(0);
//        imag(res) = 0;
    }
    else{
        res.imag(imag(A));
//        imag(res) = imag(A);
    }
    return res;
}


//methode qui multiplie 2 polynome est qui le retourne dans un vector
vector<complex<double> > multipliePolynome(vector<complex<double> > vectA, vector<complex<double> > vectB){
    
    vector<complex<double> > resultat;
    for(int i=0; i<vectA.size(); i++){
        //on multiplie simplement les cases correspondantes
        resultat.push_back(vectA.at(i) * vectB.at(i));
    }
    return resultat;
    
}

//cette methode correspond a l'algorithme de Cooley-Tukey
//qui calcule la transformee rapide de fourier d'un nombre,
//ou la transformee rapide de fourier inverse
//suivant si l'on rentre 1 ou -1 en parametre (voir le paramere inverse)
vector<complex<double> > cooleyTukey(vector<complex<double> > vect, int m, int inverse){
    int N = pow(2,m);
    int Ncourant = 2;
    int NcourantMoitie = 1;
    double nbBlocs =(double) N/2;
    complex<double> omega = 0;
    omega.real(0);
//    real(omega)=0;
    omega.imag(0);
//    imag(omega)=0;
    //définition de la variable I imaginaire pur
    complex<double> I(0.0,1.0);
    complex<double> iTheta = I*M_PI;
    complex<double> jCplx;
    complex<double> deux = 2;
    complex<double> x1;
    complex<double> x2;
    complex<double> inverseCplx = 0;
    inverseCplx.real(inverse);
//    real(inverseCplx) = inverse;
    
    
    //ici debute l'algorithme, avec les 3 boucles for
    //explique dans notre compte-rendu
    for(int cpt=1; cpt <=m; cpt++){
        for(int k=0; k<nbBlocs; k++){
            for(int j=0; j<NcourantMoitie; j++){
                //jCplx contient la valeur de j
                jCplx.real(j);
//                real(jCplx) = j;
//                imag(jCplx)=0;
                jCplx.imag(0);
                //on multiplie par inverse, qui est 1 ou -1.
                omega = exp(inverseCplx*jCplx*iTheta);
                
                
                x1 = vect.at(j+k*Ncourant);
                x2 = vect.at(j+k*Ncourant + NcourantMoitie);
                
                
                //vect contiendra les coefficients de la fft du nombre passe en parametre
                //ou ceux de la fft inverse
                vect.at(j + k*Ncourant) = x1 + omega*x2;
                vect.at(j + k*Ncourant + NcourantMoitie) = x1 - omega*x2;
                
                //on nettoie les valeurd si un residu subsiste
                vect.at(j + k*Ncourant + NcourantMoitie) = nettoyer(vect.at(j + k*Ncourant + NcourantMoitie));
                vect.at(j + k*Ncourant) = nettoyer(vect.at(j + k*Ncourant));
            }//fin sous sous for
        }//fin sous for
        
        //on modifie les variables, et on reboucle
        Ncourant = 2*Ncourant;
        NcourantMoitie = 2*NcourantMoitie;
        nbBlocs = nbBlocs/2;
        iTheta = iTheta /deux;
    }//fin for
    
    //on retourne alors le vector contenant les coefficients de la fft (ou fft inverse) du nombre
    return vect;
}





//un string representant un nombre va etre insere
//on le decoupe a l'envers, et on utilise la methode inverser
//si l'on rentre 123, on va l'inverser en 321 puis lui appliquer la methode inv.
vector<complex<double> > decouper(string s, string s2){
    
    int m = calculM(s,s2);
    int puiss2M = (int)pow(2,m);
    vector<complex<double> > vect;
    complex<double> tab[puiss2M];
    complex<double> tabInv[puiss2M];
    stringstream ss;
    int courant;
    int cpt = 0;
    int j;
    
    //boucle qui remplit les chiffres dans l'ordre inverse
    for(int i=s.length()-1; i>=0; i--){
        ss.clear();
        ss<<s.at(i);
        ss >> courant;
        tab[cpt] = (double)courant;
        cpt++;
    }
    
    //on ajoute des 0 au nombre jusqu'a 2^m
    //afin de bien realiser l'algorithme cooley-tukey
    for(int i = s.length(); i<puiss2M; i++){
        tab[i] = 0;
    }
    
    //on rempli maintenant le vector, avec la methode inverser
    //si on rentre 123 a la base, le vector contient maintenant 30201000....
    for(int i=0; i<puiss2M; i++){
        j = inv(i,m);
        tabInv[i] = tab[j];
        vect.push_back(tabInv[i]);
    }
    
    //et on retourne ce vector
    return vect;
    
}//FIN DECOUPER


//renvoie l'entier m directement superieur a la somme des degres de chaque nombre
int calculM(string a, string b){
    int m = 0;
    //m2 est 2^m
    int m2 = pow(2,m);
    int sommeDegre = a.length() + b.length();
    
    //tant que la somme des degres est superieur ou egale a m2
    while(sommeDegre >= m2){
        m++;    //on incremente m
        m2 = pow(2,m);//et on actualise la puissance
    }
    //on retourne m2
    return m;
}//FIN CALCULM


//retourne un indice apres la valeur de l'indice passe en parametre
int inv(int j, int m){
    
    int n = 0;
    for(int i=1; i<=m; i++){
        n = 2*n + j%2;
        j = j*0.5;
    }
    
    return n;
    
}//FIN INV


