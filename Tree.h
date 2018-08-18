#ifndef Tree_h
#define Tree_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "define.h"

class Tree
{
public:
  Tree();
  ~Tree();

public:
  double runnum, evnum;             // run info
  int NDataA00, NDataA01, NDataA02, NDataA03, NDataA04, NDataA05, NDataA06, NDataA07, NDataA08, NDataA09, NDataA10, NDataA11, NDataA12, NDataA13, NDataA14, NDataA15; 
  int NDataT00, NDataT01, NDataT02, NDataT03, NDataT04, NDataT05, NDataT06, NDataT07, NDataT08, NDataT09, NDataT10, NDataT11, NDataT12, NDataT13, NDataT14, NDataT15; 
  double A00[1], A01[1], A02[1], A03[1], A04[1], A05[1], A06[1], A07[1], A08[1], A09[1], A10[1], A11[1], A12[1], A13[1], A14[1], A15[1]; 
  double T00[10], T01[10], T02[10], T03[10], T04[10], T05[10], T06[10], T07[10], T08[10], T09[10], T10[10], T11[10], T12[10], T13[10], T14[10], T15[10]; 

  ///useful variables//
  double toftdc_t[nToF], toftdc_b[nToF], tofadc_t[nToF], tofadc_b[nToF];
  double actdc_t[nAC], actdc_b[nAC], acadc_t[nAC], acadc_b[nAC];
  double wctdc_t[nWC], wctdc_b[nWC], wcadc_t[nWC], wcadc_b[nWC];
  double noise_adc;

public:
  TChain *tree;

  void readtree();
  void read_one_event(int n);
  void add_tree(string ifname);
  int GetEntries()    const { return tree->GetEntries(); }

private:

};

#endif

