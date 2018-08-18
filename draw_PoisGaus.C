//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double expgaus(double *x, double *par) {
  //par[0]=Total area
  //par[1]=tau of exp function
  //par[2]=Width (sigma) of convoluted Gaussian function
  //par[3]=Shift of Function Peak
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double np = 500.0;      // number of convolution steps
  double sc =   8.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, fland, sum = 0.0, xlow, xupp, step, i;
  double val;

// Range of convolution integral
  xlow = 0.;
  xupp = x[0] + sc * par[2];
  step = (xupp-xlow) / np;
// Convolution integral
  for(i=1.0; i<=np/2; i++){
     xx = xlow + (i-0.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);

     xx = xupp - (i-.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);
  }
  val = par[0] * step * sum * invsq2pi / (par[2]*par[1]*exp(par[3]/par[1]));

  return val;
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
//////////////////////////
void draw_PoisGaus(){

    TF1* f = new TF1("f",poisgaus,500,1200,5);
    f->SetParNames("amp","mu","#sigma","offs","scale");
    f->SetParameter(0,100);
    f->SetParameter(1,1);
    f->SetParameter(2,24.1);
    f->SetParameter(3,650.);
    f->SetParameter(4,100.);
    f->SetParLimits(0,1,99999);
    f->FixParameter(1,1);
    f->SetParLimits(2,21,1000);
    f->SetParLimits(3,500,1000);
    f->SetParLimits(4,50,1000);

    //TFile *file = new TFile("replay_rootfiles/hkstest_01439_01440.root");
    TFile *file = new TFile("replay_rootfiles/hkstest_01439.root");
    TTree *tree = (TTree*)file->Get("T");
    TH1F* h = new TH1F("h","h",50,500,1200);
    tree->Project("h","D.A04","","");

    TCanvas *c=new TCanvas("c1","title",800,800);
    c->Divide(1,2);
    c->cd(1);
    //TF1* f_pois = new TF1("f_pois","TMass::Poisson(x,par[0])");
    f->SetNpx(1000);
    h->Fit(f);
    h->Draw();
    c->cd(2);
    f->Draw("");
}
