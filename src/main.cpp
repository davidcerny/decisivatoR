/*******************************R-callable function***************************************/
/*decisivatoR
Software realisation of Mareike Fischer's algorithms.
*/

#include <stdio.h>
#include <iostream>
#include <string.h>
#include "comb.h"
#include "setlib.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

RcppExport SEXP IsDecisive(SEXP taxa, SEXP s, SEXP n, SEXP _k) {
    
    short N = as<int>(n); //Number of taxons
    //std::vector<string> X(taxa.size());   
    Rcpp::CharacterVector cx = Rcpp::CharacterVector(taxa);  
    string *X = new string[N];
    for (int i=0; i<cx.size(); i++) {  
      X[i] = cx[i];  
    } 
    
    //Triplet
    string triplet[3];

    //Actual data set
    Rcpp::List SS = Rcpp::List(s);   
    int k= SS.size(); //Number of taxon subsets
    string **S = NULL;
    try {
        S = new string*[k];
    } catch (bad_alloc& ba) {
        cerr << "bad_alloc caught: " << ba.what() << endl;
    }
    for (int i=0; i<k; i++) {  
        Rcpp::CharacterVector SSS = SS[i];
        S[i] = new string[SSS.size()];
        for (int ii=0; ii<SSS.size(); ii++) {  
          S[i][ii] = SSS[ii];  
        }
    }  
    
    
    //Unrooted tree case
    short m = 3; //We consider triplets.
    
    //Number of possible combinations:
    unsigned long combnum = Res(N,m);
    cout << "Number of combinations is: " << combnum << endl;
    int *C = NULL;
    C = new int[m];
    for(int i=0; i<m; i++)
    {
      C[i] = i;
      triplet[i] = X[C[i]];
      
    }
    
    //Decisiveness. Main loop.
    //First, check the very first triplet:
    bool decisive = false;
    for(int j = 0; j < k; j++)
    {
       if(issubset(triplet, S[j], 3, 4) == true) {decisive = true;}
    }
    if(decisive == false) 
    {
      cout << "The given set is not decisive." << endl;
      return(R_NilValue);
    }
    
    for(int ii=0; ii < combnum-1; ++ii)
    {
      decisive = false;
      GetNext( C, N, m );
      
      //Printing result
      for(int jj=0; jj < m; jj++)
      {
	      triplet[jj] = X[C[jj]];
      }
      
      for(int j = 0; j < k; j++)
      {
        if(issubset(triplet, S[j], 3, 4) == true) {decisive = true;}
      }
      if(decisive == false) 
      {
        cout << "The given set is not decisive." << endl;
        return(R_NilValue);
      }
    }

    free(C);
    cout << "Computation is done. The given data set is decisive.\n";
    
    return(R_NilValue);
}
