#include "define.h"
#include "ParamMan.cc"

//____________________________________________________________________________________________
void FitGaus(TH1F *h, double &gamin, double &gamax, double range=2,int itr=5){
  gamin = h  ->GetBinCenter(h ->GetMaximumBin())-100;
  gamax = h  ->GetBinCenter(h ->GetMaximumBin())+100;
  for(Int_t l=0; l<itr; l++){
    TF1 *ga = new TF1("ga","gaus");
    ga->SetParameter(2,(gamin+gamax)/2.);
    h  ->Fit(ga,"0QR","",gamin,gamax);
    gamin = ga->GetParameter(1) - ga->GetParameter(2)*range;
    gamax = ga->GetParameter(1) + ga->GetParameter(2)*range;
    ga->Clear();
  }
  return;
}
/////////////////////////////////////////
double poisgaus(double *x, double *par){
  //par[0]=Amplitude
  //par[1]=mu of poisson
  //par[2]=Width (sigma) of convoluted Gaussian function
  //par[3]=Shift of Function Peak (offset val)
  //par[4]=Scale of Function Peak (QDCch/n.p.e. val)
  double np = 1000.0;      // number of convolution steps
  double sc =   8.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, fland, sum = 0.0, xlow, xupp, step, i, x_pois;
  double val;

  // Range of convolution integral
  //xlow = 0.;
  xlow = x[0] - sc * par[2];
  xupp = x[0] + sc * par[2];
  step = (xupp-xlow) / np;
  // Convolution integral
  for(i=1.0; i<=np/2; i++){
     xx = xlow + (i-0.5) * step;
     fland = TMath::Poisson((xx-par[3])/par[4],par[1]);
     sum += fland * TMath::Gaus(x[0],xx,par[2]);

     xx = xupp - (i-0.5) * step;
     fland = TMath::Poisson((xx-par[3])/par[4],par[1]);
     sum += fland * TMath::Gaus(x[0],xx,par[2]);
  }
  val = par[0] * step * sum;

  return val;
}
//____________________________________________________________________________________________
double fit_1pe(TH1F *h_1pe){
  double min,max,mean;
  FitGaus(h_1pe,min,max,1.,5);
  mean = 0.5*(min+max);
  cout<<mean<<endl;

  //TF1 *f_poisgaus = new TF1("f_poisgaus",poisgaus,0,1000,5);
  //f_poisgaus ->SetParNames("amp","mu","#sigma","offs","scale");
  //f_poisgaus ->SetParameter(0,100);
  //f_poisgaus ->SetParameter(1,1);
  //f_poisgaus ->SetParameter(2,24.1);
  //f_poisgaus ->SetParameter(3,mean);
  //f_poisgaus ->SetParameter(4,100.);
  //f_poisgaus ->SetParLimits(0,1,99999);
  //f_poisgaus ->FixParameter(1,1);
  //f_poisgaus ->SetParLimits(2,21,1000);
  //f_poisgaus ->SetParLimits(3,mean - 200,mean + 500);
  //f_poisgaus ->SetParLimits(4,50,1000);
  //h_1pe -> Fit(f_poisgaus,"QE","",mean - 70,mean + 500);
  //mean =f_poisgaus -> GetParameter(3);
  //cout<<mean<<f_poisgaus->GetMaximumX(mean - 70,mean + 500)<<endl;

  TF1 *f_ga = new TF1("f_ga","gaus",0,1000);
  h_1pe ->Fit(f_ga,"QR","",min,max);
  mean = f_ga->GetParameter(1);

  return mean;
}
//____________________________________________________________________________________________
double fit_ped(TH1F *h_ped){
  double min,max;
  TF1 *f_ga = new TF1("f_ga","gaus",0,1000);
  
  FitGaus(h_ped,min,max);
  h_ped ->Fit(f_ga,"QR","",min,max);
  double mean = f_ga->GetParameter(1);
  return mean;
}
//____________________________________________________________________________________________
void mk_param(){
  int runnum=1454;
  //ParamMan *paramMan = new ParamMan("aa");
  ParamMan *paramMan = new ParamMan(Form("param/%d.param",runnum));
  paramMan->Initialize();
  TString ifname1 = Form("./rootfiles/calib/%d.root",runnum);
  TFile *ifp1  = new TFile(ifname1);
  TH1F *h_ped_AC1t = (TH1F*)ifp1->Get("h_ped_AC1t");
  TH1F *h_ped_AC1b = (TH1F*)ifp1->Get("h_ped_AC1b");
  TH1F *h_1pe_AC1t = (TH1F*)ifp1->Get("h_1pe_AC1t");
  TH1F *h_1pe_AC1b = (TH1F*)ifp1->Get("h_1pe_AC1b");
  
  //paramMan->SetAdcOffset(CID_AC,1,0,fit_ped(h_ped_AC1t));
  //paramMan->SetAdcOffset(CID_AC,1,1,fit_ped(h_ped_AC1b));
  //paramMan->SetAdcPeak(CID_AC,1,0,fit_1pe(h_1pe_AC1t));
  paramMan->SetAdcPeak(CID_AC,1,1,fit_1pe(h_1pe_AC1b));
  paramMan->WriteToFile(Form("param/%d.param",runnum));
 
}
