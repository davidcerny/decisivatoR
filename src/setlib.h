/*
  This file contains some useful functions that are used in set operations:
  1. issubset(set1, set2) - checks if a set set1 is a subset of set2
  2.
*/
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <streambuf>
#include <exception>
#include <map>
#include <vector>


using namespace std;

bool issubset(string set1[], string set2[], short m, int n);
vector<int> diff(int set1[], int set2[], short m, int n);
vector<string> diff2(string set1[], string set2[], short m, int n) ;
