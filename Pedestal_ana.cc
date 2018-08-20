#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
using namespace std;

#include "TApplication.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TCut.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLatex.h"
#include "TText.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TGaxis.h"

#include "TRandom.h"

#include "define.h"
#include "Tree.h"

#define Calibration

static const double PI = 4.0*atan(1.);
static const double mrad_to_deg = 1./1000*180./PI;
const double Mp = 938.272046;          // proton       mass (MeV/c2)
const double c = 0.299792458;          // speed of light in vacuum (m/ns)

///////////////////////////////////////
double exp_ln(double *x, double *par){
  //par[0]=peak height
  //par[1]=assym
  //par[2]=width2
  //par[3]=offset

  double xx=x[0];
  double val = par[0]*TMath::Exp(-0.5*TMath::Power( TMath::Log(xx/par[1])/par[2],2. ) )+par[3];
  return val;
}

///////////////////////////////////////
class Ped_ana : public Tree
{
 public:
         Ped_ana();
        ~Ped_ana();
    void makehist();
    void loop();
    void fit();
    void draw(); 
    void savecanvas(string ofname); 
    void SetMaxEvent( int N )  { ENumMax = N; }
    void SetNoiseFlag(bool flag){ noise_flag = flag;}
    void SetRoot(string ifname);
    void SetOutputRoot(string ofname);
    bool CheckADCRange(double adc);

  private:
    TFile *ofp;
    int GetMaxEvent() { return ENumMax; }
    int ENumMax;
    int ENum;
    bool noise_flag;

    TH1F *h_ctime;
    TH1F *h_TDC_Trig1t, *h_TDC_Trig1b, *h_TDC_Trig2t, *h_TDC_Trig2b;
    TH1F *h_Q_Trig1t, *h_Q_Trig1b, *h_Q_Trig2t, *h_Q_Trig2b;
    TH2F *h2_Q_Trig1t_ev, *h2_Q_Trig1b_ev, *h2_Q_Trig2t_ev, *h2_Q_Trig2b_ev;
    TH2F *h2_Q_Trig1t2t, *h2_Q_Trig1b2b;

    TH1F *h_TDC_AC1t, *h_TDC_AC1b, *h_Q_AC1t, *h_Q_AC1b;
    TH1F *h_1pe_AC1t, *h_1pe_AC1b, *h_ped_AC1t, *h_ped_AC1b;
    TH1F *h_Qdiff_Trig1t2t, *h_Qdiff_AC1tTr1t;
    TH2F *h2_Q_AC1tb, *h2_Q_AC1tTr1t;

    TH2F *h2_QI_Trig1t, *h2_QI_Trig1b, *h2_QI_Trig2t, *h2_QI_Trig2b;
    TH2F *h2_QT_Trig1t, *h2_QT_Trig1b, *h2_QT_Trig2t, *h2_QT_Trig2b;


    int run_num;
    TCanvas *c1,*c2,*c3,*c4,*c5;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Ped_ana::Ped_ana()
{

  gErrorIgnoreLevel = kError;
  gROOT->SetStyle("Plain");
  gROOT->SetBatch(1);

  gStyle->SetOptDate(0);
  gStyle->SetOptFit(1);
  gStyle->SetHistFillStyle(3002);
  gStyle->SetHistFillColor(0);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetLineWidth(1);
  gStyle->SetOptDate(0);
//  gStyle->SetStatW(0.15);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatTextColor(1);
  gStyle->SetTitleX(0.15);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleTextColor(1);
  gStyle->SetGridWidth(0);
  gStyle->SetNdivisions(510); // tertiary*10000 + secondary*100 + first
  gStyle->SetOptStat("iMen");
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.13);

  //const Int_t NRGBs = 5;
  //const Int_t NCont = 255;
  //Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  //Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  //Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  //Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  //TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  //gStyle->SetNumberContours(NCont);
      

  c1= new TCanvas("c1","c1",1400,800 );
  c2= new TCanvas("c2","c2",1400,800 );
  c3= new TCanvas("c3","c3",1400,800 );
  c4= new TCanvas("c4","c4",1400,800 );
  c5= new TCanvas("c5","c5",1400,800 );

}
////////////////////////////////////////////////////////////////////////////
Ped_ana::~Ped_ana(){
}
////////////////////////////////////////////////////////////////////////////
void Ped_ana::SetRoot(string ifname){
  std::cout<<"SetRoot"<<std::endl;
  add_tree(ifname);
  readtree();
  ENum = GetEntries();
}
////////////////////////////////////////////////////////////////////////////
void Ped_ana::SetOutputRoot(string ofname){
  std::cout<<"SetOutputRoot"<<std::endl;
  ofp = new TFile(ofname.c_str(),"recreate");
}
////////////////////////////////////////////////////////////////////////////
void Ped_ana::makehist(){
  std::cout<<"makehist"<<std::endl;
  h_ctime    = new TH1F("h_ctime"   ,"h_ctime"   ,2400,-30,30);

  h_TDC_Trig1t    = new TH1F("h_TDC_Trig1t"   ,"h_TDC_Trig1t" ,4000,   0, 4000);
  h_TDC_Trig1b    = new TH1F("h_TDC_Trig1b"   ,"h_TDC_Trig1b" ,4000,   0, 4000);
  h_TDC_Trig2t    = new TH1F("h_TDC_Trig2t"   ,"h_TDC_Trig2t" ,4000,   0, 4000);
  h_TDC_Trig2b    = new TH1F("h_TDC_Trig2b"   ,"h_TDC_Trig2b" ,4000,   0, 4000);
  h_Q_Trig1t      = new TH1F("h_Q_Trig1t"     ,"h_Q_Trig1t"   , 500, 400,  900);
  h_Q_Trig1b      = new TH1F("h_Q_Trig1b"     ,"h_Q_Trig1b"   , 500, 400,  900);
  h_Q_Trig2t      = new TH1F("h_Q_Trig2t"     ,"h_Q_Trig2t"   , 500, 400,  900);
  h_Q_Trig2b      = new TH1F("h_Q_Trig2b"     ,"h_Q_Trig2b"   , 500, 400,  900);
  h_TDC_AC1t      = new TH1F("h_TDC_AC1t"     ,"h_TDC_AC1t"   ,4000,   0, 4000);
  h_TDC_AC1b      = new TH1F("h_TDC_AC1b"     ,"h_TDC_AC1b"   ,4000,   0, 4000);
  h_Q_AC1t        = new TH1F("h_Q_AC1t"       ,"h_Q_AC1t"     ,2000, -1000,  1000);
  h_Q_AC1b        = new TH1F("h_Q_AC1b"       ,"h_Q_AC1b"     ,2000, -1000,  1000);
  h_ped_AC1t      = new TH1F("h_ped_AC1t"       ,"h_ped_AC1t" ,2000, -1000,  1000);
  h_ped_AC1b      = new TH1F("h_ped_AC1b"       ,"h_ped_AC1b" ,2000, -1000,  1000);
  h_1pe_AC1t      = new TH1F("h_1pe_AC1t"       ,"h_1pe_AC1t" ,2000, -1000,  1000);
  h_1pe_AC1b      = new TH1F("h_1pe_AC1b"       ,"h_1pe_AC1b" ,2000, -1000,  1000);
  h_Qdiff_Trig1t2t= new TH1F("h_Qdiff_Trig1t2t" , "h_Qdiff_Trig1t2t" ,   800, -400,   400);
  h_Qdiff_AC1tTr1t= new TH1F("h_Qdiff_AC1tTr1t" , "h_Qdiff_AC1tTr1t" ,   800, -400,   400);

  h2_Q_Trig1t_ev  = new TH2F("h2_Q_Trig1t_ev", "h2_Q_Trig1t_ev", 10000,   0, 10000, 500, 400, 900);
  h2_Q_Trig1b_ev  = new TH2F("h2_Q_Trig1b_ev", "h2_Q_Trig1b_ev", 10000,   0, 10000, 500, 400, 900);
  h2_Q_Trig2t_ev  = new TH2F("h2_Q_Trig2t_ev", "h2_Q_Trig2t_ev", 10000,   0, 10000, 500, 400, 900);
  h2_Q_Trig2b_ev  = new TH2F("h2_Q_Trig2b_ev", "h2_Q_Trig2b_ev", 10000,   0, 10000, 500, 400, 900);
  h2_Q_Trig1t2t   = new TH2F("h2_Q_Trig1t2t" , "h2_Q_Trig1t2t" ,   500, 400,   900, 500, 400, 900);
  h2_Q_Trig1b2b   = new TH2F("h2_Q_Trig1b2b" , "h2_Q_Trig1b2b" ,   500, 400,   900, 500, 400, 900);
  h2_Q_AC1tTr1t   = new TH2F("h2_Q_AC1tTr1t" , "h2_Q_AC1tTr1t" ,   500, 400,   900, 500, 400, 900);
  h2_Q_AC1tb      = new TH2F("h2_Q_AC1tb"    , "h2_Q_AC1tb"    ,   500, 400,   900, 500, 400, 900);
  h2_QI_Trig1t    = new TH2F("h2_QI_Trig1t","h2_QI_Trig1t",200,   0, 50,1501,5000,15000);
  h2_QI_Trig1b    = new TH2F("h2_QI_Trig1b","h2_QI_Trig1b",200,   0, 50,1501,5000,15000);
  h2_QI_Trig2t    = new TH2F("h2_QI_Trig2t","h2_QI_Trig2t",200,   0, 50,1501,5000,15000);
  h2_QI_Trig2b    = new TH2F("h2_QI_Trig2b","h2_QI_Trig2b",200,   0, 50,1501,5000,15000);
  h2_QT_Trig1t    = new TH2F("h2_QT_Trig1t","h2_QT_Trig1t",200,-100,100,1501,5000,15000);
  h2_QT_Trig1b    = new TH2F("h2_QT_Trig1b","h2_QT_Trig1b",200,-100,100,1501,5000,15000);
  h2_QT_Trig2t    = new TH2F("h2_QT_Trig2t","h2_QT_Trig2t",200,-100,100,1501,5000,15000);
  h2_QT_Trig2b    = new TH2F("h2_QT_Trig2b","h2_QT_Trig2b",200,-100,100,1501,5000,15000);

  
}
////////////////////////////////////////////////////////////////////////////
void Ped_ana::loop(){

  if( GetMaxEvent()>0 && GetMaxEvent()<ENum) ENum = GetMaxEvent();
  for(int n=0;n<ENum;n++){
    read_one_event(n);
    if(n%100==0) cout<<n<<" / "<<ENum<<endl;
    if(noise_flag && (NDataA07!=1 || NDataA04!=1 || NDataA05!=1) && !CheckADCRange(noise_adc))continue;
    if(!noise_flag && (NDataA04!=1 || NDataA05!=1))continue;

    ////////
    //Trig//
    ////////
    h_TDC_Trig1t    ->Fill(toftdc_t[0]);
    h_TDC_Trig1b    ->Fill(toftdc_b[0]);
    h_TDC_Trig2t    ->Fill(toftdc_t[1]);
    h_TDC_Trig2b    ->Fill(toftdc_b[1]);
    if(NDataT00==0&&NDataA00==1){h_Q_Trig1t      ->Fill(tofadc_t[0]);h2_Q_Trig1t_ev ->Fill(evnum, tofadc_t[0]);}
    if(NDataT01==0&&NDataA01==1){h_Q_Trig1b      ->Fill(tofadc_b[0]);h2_Q_Trig1b_ev ->Fill(evnum, tofadc_b[0]);}
    if(NDataT02==0&&NDataA02==1){h_Q_Trig2t      ->Fill(tofadc_t[1]);h2_Q_Trig2t_ev ->Fill(evnum, tofadc_t[1]);}
    if(NDataT03==0&&NDataA03==1){h_Q_Trig2b      ->Fill(tofadc_b[1]);h2_Q_Trig2b_ev ->Fill(evnum, tofadc_b[1]);}
    if(NDataT00==0&&NDataT02==0&&NDataA00==1&&NDataA02==1){h2_Q_Trig1t2t      ->Fill(tofadc_t[0],tofadc_t[1]);h_Qdiff_Trig1t2t      ->Fill(tofadc_t[0]-tofadc_t[1]);}
    if(NDataT01==0&&NDataT03==0)h2_Q_Trig1b2b      ->Fill(tofadc_b[0],tofadc_b[1]);
    //h2_QI_Trig1t   ->Fill(peak[0]-offset[0],Trig_sum[0]);
    //h2_QI_Trig1b   ->Fill(peak[1]-offset[1],Trig_sum[1]);
    //h2_QI_Trig2t   ->Fill(peak[2]-offset[2],Trig_sum[2]);
    //h2_QI_Trig2b   ->Fill(peak[3]-offset[3],Trig_sum[3]);

    //////
    //AC//
    //////
    //h_TDC_AC1t    ->Fill(actdc_t[0]);
    //h_TDC_AC1b    ->Fill(actdc_b[0]);
    if(noise_flag){
      if(CheckADCRange(acadc_t[0])){
        h_Q_AC1t      ->Fill(acadc_t[0] - noise_adc);
        if(acadc_t[0] - noise_adc>102. &&acadc_t[0] - noise_adc<104.)cout<<acadc_t[0] <<" "<<noise_adc<<endl;
        if(NDataT04==0)h_ped_AC1t      ->Fill(acadc_t[0] - noise_adc);
        if(NDataT04==1)h_1pe_AC1t      ->Fill(acadc_t[0] - noise_adc);
      }
      if(CheckADCRange(acadc_b[0])){
        h_Q_AC1b      ->Fill(acadc_b[0] - noise_adc);
        if(NDataT05==0)h_ped_AC1b      ->Fill(acadc_b[0] - noise_adc);
        if(NDataT05==1)h_1pe_AC1b      ->Fill(acadc_b[0] - noise_adc);
      }
    }
    else{
      if(NDataT04==0 && CheckADCRange(acadc_t[0]))h_Q_AC1t      ->Fill(acadc_t[0]);
      if(NDataT05==0 && CheckADCRange(acadc_b[0]))h_Q_AC1b      ->Fill(acadc_b[0]);
    }

    if(NDataT00==0&&NDataT04==0){h2_Q_AC1tTr1t      ->Fill(acadc_t[0],tofadc_t[0]);h_Qdiff_AC1tTr1t      ->Fill(acadc_t[0]-tofadc_t[0]);}
    if(NDataT04==0&&NDataT05==0)h2_Q_AC1tb      ->Fill(acadc_t[0],acadc_b[0]);

  }//event loop

  ofp->Write();
}
////////////////////////////////////////////////////////////////////////////
void Ped_ana::fit(){
  //double top_min = h_Q_AC1t->GetBinCenter(h_Q_AC1t->GetMaximumBin())-100;
  //double top_max = h_Q_AC1t->GetBinCenter(h_Q_AC1t->GetMaximumBin())+300;
  //double bot_min = h_Q_AC1b->GetBinCenter(h_Q_AC1b->GetMaximumBin())-100;
  //double bot_max = h_Q_AC1b->GetBinCenter(h_Q_AC1b->GetMaximumBin())+300;
  double top_min = 0;
  double top_max = 500;
  double bot_min = 0;
  double bot_max = 500;
  h_Q_AC1t   ->GetXaxis() -> SetRangeUser(top_min, top_max);
  h_ped_AC1t ->GetXaxis() -> SetRangeUser(top_min, top_max);
  h_1pe_AC1t ->GetXaxis() -> SetRangeUser(top_min, top_max);
  h_Q_AC1b   ->GetXaxis() -> SetRangeUser(bot_min, bot_max);
  h_ped_AC1b ->GetXaxis() -> SetRangeUser(bot_min, bot_max);
  h_1pe_AC1b ->GetXaxis() -> SetRangeUser(bot_min, bot_max);
}
////////////////////////////////////////////////////////////////////////////
void Ped_ana::draw(){

  c1->Clear();
  c1->Divide(4,4);
  c1->cd(1); gPad->SetLogy(1);h_TDC_Trig1t    ->Draw();
  c1->cd(2); gPad->SetLogy(1);h_TDC_Trig1b    ->Draw();
  c1->cd(3); gPad->SetLogy(1);h_TDC_Trig2t    ->Draw();
  c1->cd(4); gPad->SetLogy(1);h_TDC_Trig2b    ->Draw();
  c1->cd(5); gPad->SetLogy(1);h_Q_Trig1t      ->Draw();
  c1->cd(6); gPad->SetLogy(1);h_Q_Trig1b      ->Draw();
  c1->cd(7); gPad->SetLogy(1);h_Q_Trig2t      ->Draw();
  c1->cd(8); gPad->SetLogy(1);h_Q_Trig2b      ->Draw();
  c1->cd(9); gPad->SetLogy(1);h_Q_Trig2b      ->Draw();
  c1->cd(10); gPad->SetLogy(0);h2_Q_Trig1t_ev ->Draw();
  c1->cd(11); gPad->SetLogy(0);h2_Q_Trig1b_ev ->Draw();
  c1->cd(12); gPad->SetLogy(0);h2_Q_Trig2t_ev ->Draw();
  c1->cd(13); gPad->SetLogy(0);h2_Q_Trig2b_ev ->Draw();
  c1->cd(14); gPad->SetLogy(0);h2_Q_Trig1t2t  ->Draw();
  c1->cd(15); gPad->SetLogy(0);h2_Q_Trig1b2b  ->Draw();
  c1->cd(15); gPad->SetLogy(1);h_Qdiff_Trig1t2t  ->Draw();
  c2->Clear();
  c2->Divide(4,4);
  c2->cd(1); gPad->SetLogy(1);h_TDC_AC1t   ->Draw();
  c2->cd(2); gPad->SetLogy(1);h_TDC_AC1b   ->Draw();
  c2->cd(3); gPad->SetLogy(1);h_Q_AC1t     ->Draw();
  c2->cd(4); gPad->SetLogy(1);h_Q_AC1b     ->Draw();
  c2->cd(5); h2_Q_AC1tTr1t      ->Draw("");
  c2->cd(5); gPad->SetLogy(1);h_Qdiff_AC1tTr1t      ->Draw("");
  c2->cd(7); gPad->SetLogy(1);h_ped_AC1t     ->Draw();
  c2->cd(8); gPad->SetLogy(1);h_ped_AC1b     ->Draw();
  c2->cd(11);gPad->SetLogy(1);h_1pe_AC1t     ->Draw();
  c2->cd(12);gPad->SetLogy(1);h_1pe_AC1b     ->Draw();
  c3->Clear();
  c4->Clear();
  c5->Clear();

}
////////////////////////////////////////////////////////////////////////////
void Ped_ana::savecanvas(string ofname){
  //string ofname = "test.root";//tmp
  string ofname_pdf = ofname;
  ofname_pdf.erase(ofname_pdf.size()-5);
  ofname_pdf.append(".pdf");
  c1->Print(Form("%s[",ofname_pdf.c_str()));
  c1->Print(Form("%s" ,ofname_pdf.c_str()));
  c2->Print(Form("%s" ,ofname_pdf.c_str()));
  c2->Print(Form("%s]",ofname_pdf.c_str()));
  cout<<ofname_pdf.c_str()<<" saved"<<endl;
}
////////////////////////////////////////////////////////////////////////////
bool Ped_ana::CheckADCRange(double adc){
  if(adc>400. && adc<4000.)return true;
  else return false;
}
////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  string ifname = "rootfiles/fadcAERO_117.root";
  string ofname = "rootfiles/conv/hoge.root";
  int ch;
  int MaxNum = 0;
  bool output_flag = false;
  bool draw_flag = true;
  bool noise_ext_flag = false;
  extern char *optarg;
  while((ch=getopt(argc,argv,"hf:w:n:be"))!=-1){
    switch(ch){
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      break;
    case 'w':
      output_flag = true;
      draw_flag = true;
      ofname = optarg;
      cout<<"output filename : "<<ofname<<endl;
      break;
    case 'n':
      MaxNum = atoi(optarg);
      break;
    case 'b':
      draw_flag = false;
      cout<<"BACH MODE!"<<endl;
      break;
    case 'e':
      noise_ext_flag = true;
      cout<<"Noise cancel mode!"<<endl;
      break;
    case 'h':
      cout<<"-f : input root file name"<<endl;
      cout<<"-w : output root file name"<<endl;
      cout<<"-n : maximum number of analysed events"<<endl;
      return 0;
      break;
    case '?':
      cout<<"unknown option...."<<endl;
      return 0;
      break;
    default:
      cout<<"type -h to see help!!"<<endl;
      return 0;
    }
  }

  cout<<"before TApplication"<<endl;
  TApplication *theApp = new TApplication("App", &argc, argv);
  cout<<"after TApplication"<<endl;
  Ped_ana *ana = new Ped_ana();
  cout<<"after class Ped_ana constructor"<<endl;

  ana->SetMaxEvent(MaxNum);
  ana->SetRoot(ifname);
  ana->SetOutputRoot(ofname);
  ana->SetNoiseFlag(noise_ext_flag); 
  ana->makehist();
  ana->loop();
  ana->fit();
  if(draw_flag){
    ana->draw();
    ana->savecanvas(ofname);
  }
  delete ana;

  gSystem->Exit(1);
  theApp->Run();
  return 0;
}

