#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;

#include "Tree.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
Tree::Tree()
{
  tree = new TChain("T");
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
Tree::~Tree()
{
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void Tree::add_tree(string ifname)
{
  tree->Add(ifname.c_str());
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void Tree::readtree()
{
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("g.runnum"         ,1);  tree->SetBranchAddress("g.runnum"         ,&runnum     );
  tree->SetBranchStatus("g.evnum"          ,1);  tree->SetBranchAddress("g.evnum"          ,&evnum      );
  tree->SetBranchStatus("Ndata.D.A00"      ,1);  tree->SetBranchAddress("Ndata.D.A00"      ,&NDataA00   );
  tree->SetBranchStatus("Ndata.D.A01"      ,1);  tree->SetBranchAddress("Ndata.D.A01"      ,&NDataA01   );
  tree->SetBranchStatus("Ndata.D.A02"      ,1);  tree->SetBranchAddress("Ndata.D.A02"      ,&NDataA02   );
  tree->SetBranchStatus("Ndata.D.A03"      ,1);  tree->SetBranchAddress("Ndata.D.A03"      ,&NDataA03   );
  tree->SetBranchStatus("Ndata.D.A04"      ,1);  tree->SetBranchAddress("Ndata.D.A04"      ,&NDataA04   );
  tree->SetBranchStatus("Ndata.D.A05"      ,1);  tree->SetBranchAddress("Ndata.D.A05"      ,&NDataA05   );
  tree->SetBranchStatus("Ndata.D.A06"      ,1);  tree->SetBranchAddress("Ndata.D.A06"      ,&NDataA06   );
  tree->SetBranchStatus("Ndata.D.A07"      ,1);  tree->SetBranchAddress("Ndata.D.A07"      ,&NDataA07   );
  tree->SetBranchStatus("D.A00"            ,1);  tree->SetBranchAddress("D.A00"            , A00        );
  tree->SetBranchStatus("D.A01"            ,1);  tree->SetBranchAddress("D.A01"            , A01        );
  tree->SetBranchStatus("D.A02"            ,1);  tree->SetBranchAddress("D.A02"            , A02        );
  tree->SetBranchStatus("D.A03"            ,1);  tree->SetBranchAddress("D.A03"            , A03        );
  tree->SetBranchStatus("D.A04"            ,1);  tree->SetBranchAddress("D.A04"            , A04        );
  tree->SetBranchStatus("D.A05"            ,1);  tree->SetBranchAddress("D.A05"            , A05        );
  tree->SetBranchStatus("D.A07"            ,1);  tree->SetBranchAddress("D.A07"            , A07        );
  tree->SetBranchStatus("Ndata.D.TH00"     ,1);  tree->SetBranchAddress("Ndata.D.TH00"     ,&NDataT00   );
  tree->SetBranchStatus("Ndata.D.TH01"     ,1);  tree->SetBranchAddress("Ndata.D.TH01"     ,&NDataT01   );
  tree->SetBranchStatus("Ndata.D.TH02"     ,1);  tree->SetBranchAddress("Ndata.D.TH02"     ,&NDataT02   );
  tree->SetBranchStatus("Ndata.D.TH03"     ,1);  tree->SetBranchAddress("Ndata.D.TH03"     ,&NDataT03   );
  tree->SetBranchStatus("Ndata.D.TH04"     ,1);  tree->SetBranchAddress("Ndata.D.TH04"     ,&NDataT04   );
  tree->SetBranchStatus("Ndata.D.TH05"     ,1);  tree->SetBranchAddress("Ndata.D.TH05"     ,&NDataT05   );
  tree->SetBranchStatus("Ndata.D.TH06"     ,1);  tree->SetBranchAddress("Ndata.D.TH06"     ,&NDataT06   );
  tree->SetBranchStatus("Ndata.D.TH07"     ,1);  tree->SetBranchAddress("Ndata.D.TH07"     ,&NDataT07   );
  tree->SetBranchStatus("D.TH00"           ,1);  tree->SetBranchAddress("D.TH00"           , T00        );
  tree->SetBranchStatus("D.TH01"           ,1);  tree->SetBranchAddress("D.TH01"           , T01        );
  tree->SetBranchStatus("D.TH02"           ,1);  tree->SetBranchAddress("D.TH02"           , T02        );
  tree->SetBranchStatus("D.TH03"           ,1);  tree->SetBranchAddress("D.TH03"           , T03        );
  tree->SetBranchStatus("D.TH04"           ,1);  tree->SetBranchAddress("D.TH04"           , T04        );
  tree->SetBranchStatus("D.TH05"           ,1);  tree->SetBranchAddress("D.TH05"           , T05        );
  tree->SetBranchStatus("D.TH07"           ,1);  tree->SetBranchAddress("D.TH07"           , T07        );
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void Tree::read_one_event(int n)
{
  tree->GetEntry(n);
  toftdc_t[0]=T00[0];tofadc_t[0]=A00[0];
  toftdc_b[0]=T03[0];tofadc_b[0]=A01[0];
  toftdc_t[1]=T02[0];tofadc_t[1]=A02[0];
  toftdc_b[1]=T01[0];tofadc_b[1]=A03[0];

  actdc_t[0]=T04[0];acadc_t[0]=A04[0];
  actdc_b[0]=T05[0];acadc_b[0]=A05[0];
  noise_adc = A07[0];

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
