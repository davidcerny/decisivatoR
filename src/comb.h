/*
  Combination: generates all combinations from given N and M.
  Source: "Programmirovanie v algoritmakh" by S. Okulov
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <streambuf>
#include <exception>

using namespace std;

unsigned int Res(unsigned int N, unsigned int k);
void GetNext(int *C, unsigned int N, unsigned int k );
