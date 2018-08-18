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
#include "ParamMan.h"

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
struct TreeBranch{
  int runnum, evnum;
  double trig_time_t[2], trig_time_b[2];
  double trig_adc_t[2], trig_adc_b[2];
  double trig_tof;
  double AC_time_t, AC_time_b;
  double AC_adc_t, AC_adc_b;
  double AC_npe_t, AC_npe_b, AC_npe_sum; 
  int NDataA[6],NDataT[6];
};

static TreeBranch tr;         

///////////////////////////////////////
class AC_ana : public Tree
{
 public:
         AC_ana();
        ~AC_ana();
    void makehist();
    void loop();
    void fit();
    void draw(); 
    void savecanvas(string ofname); 
    void SetMaxEvent( int N )  { ENumMax = N; }
    void SetNoiseFlag(bool flag){ noise_flag = flag;}
    void SetRoot(string ifname);
    void SetParamFile(string paramname);
    void SetOutputRoot(string ofname);

  private:
    ParamMan *paramMan;
    TFile *ofp;
    TTree *tree_out;
    string input_param;
    int GetMaxEvent() { return ENumMax; }
    int ENumMax;
    int ENum;
    bool noise_flag;

    TH1F *h_ctime;
    TH1F *h_TDC_Trig1t, *h_TDC_Trig1b, *h_TDC_Trig2t, *h_TDC_Trig2b;
    TH1F *h_Q_Trig1t, *h_Q_Trig1b, *h_Q_Trig2t, *h_Q_Trig2b;

    TH1F *h_TDC_AC1t, *h_TDC_AC1b, *h_Q_AC1t, *h_Q_AC1b;
    TH1F *h_npe_AC1t, *h_npe_AC1b, *h_npe_AC1sum;

    TH2F *h2_QI_Trig1t, *h2_QI_Trig1b, *h2_QI_Trig2t, *h2_QI_Trig2b;
    TH2F *h2_QT_Trig1t, *h2_QT_Trig1b, *h2_QT_Trig2t, *h2_QT_Trig2b;

    int run_num;
    TCanvas *c1,*c2,*c3,*c4,*c5;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AC_ana::AC_ana()
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
AC_ana::~AC_ana(){
}
////////////////////////////////////////////////////////////////////////////
void AC_ana::SetRoot(string ifname){
  std::cout<<"SetRoot"<<std::endl;
  add_tree(ifname);
  readtree();
  ENum = GetEntries();
}
////////////////////////////////////////////////////////////////////////////
void AC_ana::SetParamFile(string paramname){
  input_param = paramname;
  cout<<"input param name: "<<input_param<<endl;
  paramMan= new ParamMan(input_param.c_str());
  paramMan->Initialize();
}
////////////////////////////////////////////////////////////////////////////
void AC_ana::SetOutputRoot(string ofname){
  std::cout<<"SetOutputRoot"<<std::endl;
  ofp = new TFile(ofname.c_str(),"recreate");
  tree_out = new TTree("tree","tree");
  tree_out->Branch("runnum"      ,&tr.runnum      ,"runnum/I"         );
  tree_out->Branch("evnum"       ,&tr.evnum       ,"evnum/I"          );
  tree_out->Branch("trig_time_t" , tr.trig_time_t ,"trig_time_t[2]/D" );
  tree_out->Branch("trig_time_b" , tr.trig_time_b ,"trig_time_b[2]/D" );
  tree_out->Branch("trig_adc_t"  , tr.trig_adc_t  ,"trig_adc_t[2]/D"  );
  tree_out->Branch("trig_adc_b"  , tr.trig_adc_b  ,"trig_adc_b[2]/D"  );
  tree_out->Branch("trig_tof"    ,&tr.trig_tof    ,"trig_tof/D"     );
  tree_out->Branch("AC_time_t"   ,&tr.AC_time_t   ,"AC_time_t/D"    );
  tree_out->Branch("AC_time_b"   ,&tr.AC_time_b   ,"AC_time_b/D"    );
  tree_out->Branch("AC_adc_t"    ,&tr.AC_adc_t    ,"AC_adc_t/D"     );
  tree_out->Branch("AC_adc_b"    ,&tr.AC_adc_b    ,"AC_adc_b/D"     );
  tree_out->Branch("AC_npe_t"    ,&tr.AC_npe_t    ,"AC_npe_t/D"     );
  tree_out->Branch("AC_npe_b"    ,&tr.AC_npe_b    ,"AC_npe_b/D"     );
  tree_out->Branch("AC_npe_sum"  ,&tr.AC_npe_sum  ,"AC_npe_sum/D"   );
  tree_out->Branch("NDataA"      , tr.NDataA      ,"NDataA[6]/I"   );
  tree_out->Branch("NDataT"      , tr.NDataT      ,"NDataT[6]/I"   );
}
////////////////////////////////////////////////////////////////////////////
void AC_ana::makehist(){
  std::cout<<"makehist"<<std::endl;
  h_ctime    = new TH1F("h_ctime"   ,"h_ctime"   ,2400,-30,30);

  h_TDC_Trig1t    = new TH1F("h_TDC_Trig1t"   ,"h_TDC_Trig1t" ,4000,   0, 4000);
  h_TDC_Trig1b    = new TH1F("h_TDC_Trig1b"   ,"h_TDC_Trig1b" ,4000,   0, 4000);
  h_TDC_Trig2t    = new TH1F("h_TDC_Trig2t"   ,"h_TDC_Trig2t" ,4000,   0, 4000);
  h_TDC_Trig2b    = new TH1F("h_TDC_Trig2b"   ,"h_TDC_Trig2b" ,4000,   0, 4000);
  h_Q_Trig1t      = new TH1F("h_Q_Trig1t"     ,"h_Q_Trig1t"   ,1500,   0, 1500);
  h_Q_Trig1b      = new TH1F("h_Q_Trig1b"     ,"h_Q_Trig1b"   ,1500,   0, 1500);
  h_Q_Trig2t      = new TH1F("h_Q_Trig2t"     ,"h_Q_Trig2t"   ,1500,   0, 1500);
  h_Q_Trig2b      = new TH1F("h_Q_Trig2b"     ,"h_Q_Trig2b"   ,1500,   0, 1500);
  h_TDC_AC1t      = new TH1F("h_TDC_AC1t"     ,"h_TDC_AC1t"   ,4000,   0, 4000);
  h_TDC_AC1b      = new TH1F("h_TDC_AC1b"     ,"h_TDC_AC1b"   ,4000,   0, 4000);
  h_Q_AC1t        = new TH1F("h_Q_AC1t"       ,"h_Q_AC1t"     ,1500,   0, 1500);
  h_Q_AC1b        = new TH1F("h_Q_AC1b"       ,"h_Q_AC1b"     ,1500,   0, 1500);
  h_npe_AC1t      = new TH1F("h_npe_AC1t"     ,"h_npe_AC1t"   , 100,  -2,   20);
  h_npe_AC1b      = new TH1F("h_npe_AC1b"     ,"h_npe_AC1b"   , 100,  -2,   20);
  h_npe_AC1sum    = new TH1F("h_npe_AC1sum"   ,"h_npe_AC1sum" , 250,  -2,   50);

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
void AC_ana::loop(){

  if( GetMaxEvent()>0 && GetMaxEvent()<ENum) ENum = GetMaxEvent();
  for(int n=0;n<ENum;n++){
    //tree->GetEntry(n);
    read_one_event(n);
    if(n%100==0) cout<<n<<" / "<<ENum<<endl;
    double ac_npe_t , ac_npe_b , ac_npe_sum;
    double ac_time_t , ac_time_b;
    double tr_time_t[2], tr_time_b[2];

    ac_npe_t = ac_npe_b = ac_npe_sum = -99.;
    ac_time_t = ac_time_b = -99.;
    for(int i=0;i<2;i++){//seg
      tr_time_t[i] = tr_time_b[i] = -99.;
    }
   
    ////////
    //Trig//
    ////////
    h_TDC_Trig1t    ->Fill(toftdc_t[0]);
    h_TDC_Trig1b    ->Fill(toftdc_b[0]);
    h_TDC_Trig2t    ->Fill(toftdc_t[1]);
    h_TDC_Trig2b    ->Fill(toftdc_b[1]);
    h_Q_Trig1t      ->Fill(tofadc_t[0]);
    h_Q_Trig1b      ->Fill(tofadc_b[0]);
    h_Q_Trig2t      ->Fill(tofadc_t[1]);
    h_Q_Trig2b      ->Fill(tofadc_b[1]);
    //h2_QI_Trig1t   ->Fill(peak[0]-offset[0],Trig_sum[0]);
    //h2_QI_Trig1b   ->Fill(peak[1]-offset[1],Trig_sum[1]);
    //h2_QI_Trig2t   ->Fill(peak[2]-offset[2],Trig_sum[2]);
    //h2_QI_Trig2b   ->Fill(peak[3]-offset[3],Trig_sum[3]);
    for(int i=0;i<2;i++){//seg
      tr_time_t[i] = paramMan->time(CID_ToF, i+1,0,toftdc_t[i]);
      tr_time_b[i] = paramMan->time(CID_ToF, i+1,1,toftdc_b[i]);
    }

    //////
    //AC//
    //////
    if(noise_flag){//noise ch(A07 ) is available after run1425
      ac_npe_t = paramMan->npe(CID_AC, 1,0, acadc_t[0] - noise_adc);
      ac_npe_b = paramMan->npe(CID_AC, 1,1, acadc_b[0] - noise_adc);
    }
    else{
      ac_npe_t = paramMan->npe(CID_AC, 1,0, acadc_t[0]);
      ac_npe_b = paramMan->npe(CID_AC, 1,1, acadc_b[0]);
    }
    ac_npe_sum = ac_npe_t + ac_npe_b;
    ac_time_t = paramMan->time(CID_AC, 1,0, actdc_t[0]);
    ac_time_b = paramMan->time(CID_AC, 1,1, actdc_b[0]);
    h_TDC_AC1t    ->Fill(actdc_t[0]);
    h_TDC_AC1b    ->Fill(actdc_b[0]);
    h_Q_AC1t      ->Fill(acadc_t[0]);
    h_Q_AC1b      ->Fill(acadc_b[0]);
    //if(NDataT04>0 && NDataT05>0){
      h_npe_AC1t    ->Fill(ac_npe_t);
      h_npe_AC1b    ->Fill(ac_npe_b);
      h_npe_AC1sum  ->Fill(ac_npe_sum);
    //}

    ///////////////
    //Remake tree//
    ///////////////
    if(tree_out!=0){ //this is tentative. you can change conditions to remake tree
      tr.runnum = runnum;
      tr.evnum = evnum;
      for(int i=0;i<2;i++){
        tr.trig_time_t[i] = tr_time_t[i];
        tr.trig_time_b[i] = tr_time_b[i];
        tr.trig_adc_t[i]  = tofadc_t[i];
        tr.trig_adc_b[i]  = tofadc_b[i];
      }
      tr.trig_tof   = (tr_time_t[0]+tr_time_b[0])/2. -(tr_time_t[1]+tr_time_b[1])/2.;  
      tr.AC_time_t  = ac_time_t;
      tr.AC_time_b  = ac_time_b;
      tr.AC_adc_t   = acadc_t[0];
      tr.AC_adc_b   = acadc_b[0];
      tr.AC_npe_t   = ac_npe_t;
      tr.AC_npe_b   = ac_npe_b;
      tr.AC_npe_sum = ac_npe_sum;
      tr.NDataA[0]  = NDataA00;  tr.NDataT[0]  = NDataT00;
      tr.NDataA[1]  = NDataA01;  tr.NDataT[1]  = NDataT01;
      tr.NDataA[2]  = NDataA02;  tr.NDataT[2]  = NDataT02;
      tr.NDataA[3]  = NDataA03;  tr.NDataT[3]  = NDataT03;
      tr.NDataA[4]  = NDataA04;  tr.NDataT[4]  = NDataT04;
      tr.NDataA[5]  = NDataA05;  tr.NDataT[5]  = NDataT05;
      tree_out -> Fill();
    }
  }//event loop

  ofp->Write();
}
////////////////////////////////////////////////////////////////////////////
void AC_ana::fit(){
}
////////////////////////////////////////////////////////////////////////////
void AC_ana::draw(){

  c1->Clear();
  c1->Divide(4,4);
  c1->cd(1); gPad->SetLogy(1);h_TDC_Trig1t ->Draw();
  c1->cd(2); gPad->SetLogy(1);h_TDC_Trig1b ->Draw();
  c1->cd(3); gPad->SetLogy(1);h_TDC_Trig2t ->Draw();
  c1->cd(4); gPad->SetLogy(1);h_TDC_Trig2b ->Draw();
  c1->cd(5); gPad->SetLogy(1);h_Q_Trig1t   ->Draw();
  c1->cd(6); gPad->SetLogy(1);h_Q_Trig1b   ->Draw();
  c1->cd(7); gPad->SetLogy(1);h_Q_Trig2t   ->Draw();
  c1->cd(8); gPad->SetLogy(1);h_Q_Trig2b   ->Draw();
  c2->Clear();
  c2->Divide(4,4);
  c2->cd(1); gPad->SetLogy(1);h_TDC_AC1t   ->Draw();
  c2->cd(2); gPad->SetLogy(1);h_TDC_AC1b   ->Draw();
  c2->cd(3); gPad->SetLogy(1);h_Q_AC1t     ->Draw();
  c2->cd(4); gPad->SetLogy(1);h_Q_AC1b     ->Draw();
  c2->cd(5); gPad->SetLogy(1);h_npe_AC1t   ->Draw();
  c2->cd(6); gPad->SetLogy(1);h_npe_AC1b   ->Draw();
  c2->cd(7); gPad->SetLogy(1);h_npe_AC1sum ->Draw();
  c3->Clear();
  c4->Clear();
  c5->Clear();

}
////////////////////////////////////////////////////////////////////////////
void AC_ana::savecanvas(string ofname){
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
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  string ifname = "rootfiles/fadcAERO_117.root";
  string ofname = "rootfiles/conv/hoge.root";
  string paramname="param/default.param";
  int ch;
  int MaxNum = 0;
  bool output_flag = false;
  bool draw_flag = true;
  bool noise_ext_flag = false;
  extern char *optarg;
  while((ch=getopt(argc,argv,"hf:w:n:p:be"))!=-1){
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
    case 'p':
      paramname = optarg;
      break;
    case 'e':
      noise_ext_flag = true;
      cout<<"Noise cancel mode!"<<endl;
      break;
    case 'h':
      cout<<"-f : input root file name"<<endl;
      cout<<"-w : output root file name"<<endl;
      cout<<"-n : maximum number of analysed events"<<endl;
      cout<<"-p : input parameter file name"<<endl;
      cout<<"-e : noise cancel mode"<<endl;
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

  //cout<<"before TApplication"<<endl;
  TApplication *theApp = new TApplication("App", &argc, argv);
  //cout<<"after TApplication"<<endl;
  AC_ana *ana = new AC_ana();
  //cout<<"after class AC_ana constructor"<<endl;

  ana->SetMaxEvent(MaxNum);
  ana->SetParamFile(paramname);
  ana->SetRoot(ifname);
  ana->SetOutputRoot(ofname); 
  ana->SetNoiseFlag(noise_ext_flag); 
  ana->makehist();
  
  ana->loop();
  //ana->fit();
  if(draw_flag){
    ana->draw();
    ana->savecanvas(ofname);
  }
  delete ana;

  gSystem->Exit(1);
  theApp->Run();
  return 0;
}

