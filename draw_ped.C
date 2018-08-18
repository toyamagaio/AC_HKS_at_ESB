
void SetTH2(TH2 *h, TString name, TString xname, TString yname, double min=1, double MStyle=1, double MSize=0.5){
  h->SetTitle(name);
  h->SetMinimum(min);
  h->SetLineWidth(0);
  h->SetTitleSize(0.05,"");
  h->SetMarkerStyle(1);
  h->SetMarkerSize(0.1);
  h->SetMarkerColor(1);

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(1.0);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(0.8);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(3);

  h->SetStats(0);
}

//____________________________________________________________________________________________
void SetTH1(TH1 *h, TString hname, TString xname, TString yname, int LColor, int FStyle, int FColor){
  h->SetTitle(hname);
  h->SetLineColor(LColor);
  h->SetLineWidth(1);
  h->SetFillStyle(FStyle);
  h->SetFillColor(FColor);
//  h->SetMinimum(0.8);

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(1.0);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(0.8);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(3);
}
/////////////////////////
void draw_ped(){
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadBottomMargin(0.15);

  TString ifname1=Form("hoge.root");//
  TFile *ifp1  = new TFile(ifname1  );
  cout<<"input filename1 : "<<ifname1<<endl;
  TH1F *h_Q_Trig1t, *h_Q_Trig1b, *h_Q_Trig2t, *h_Q_Trig2b;
  TH2F *h2_Q_Trig1t_ev, *h2_Q_Trig1b_ev, *h2_Q_Trig2t_ev, *h2_Q_Trig2b_ev;
  TH2F *h2_Q_Trig1t2t, *h2_Q_Trig1b2b;
  TH1F *h_TDC_AC1t, *h_TDC_AC1b, *h_Q_AC1t, *h_Q_AC1b;
  TH2F *h2_Q_AC1tb, *h2_Q_AC1tTr1t;
  h_Q_Trig1t      = (TH1F*)ifp1->Get("h_Q_Trig1t");
  h_Q_Trig1b      = (TH1F*)ifp1->Get("h_Q_Trig1b");
  h_Q_Trig2t      = (TH1F*)ifp1->Get("h_Q_Trig2t");
  h_Q_Trig2b      = (TH1F*)ifp1->Get("h_Q_Trig2b");
  h2_Q_Trig1t_ev  = (TH2F*)ifp1->Get("h2_Q_Trig1t_ev");
  h2_Q_Trig1b_ev  = (TH2F*)ifp1->Get("h2_Q_Trig1b_ev");
  h2_Q_Trig2t_ev  = (TH2F*)ifp1->Get("h2_Q_Trig2t_ev");
  h2_Q_Trig2b_ev  = (TH2F*)ifp1->Get("h2_Q_Trig2b_ev");
  h2_Q_Trig1t2t   = (TH2F*)ifp1->Get("h2_Q_Trig1t2t" );
  h2_Q_Trig1b2b   = (TH2F*)ifp1->Get("h2_Q_Trig1b2b" );
  h2_Q_AC1tTr1t   = (TH2F*)ifp1->Get("h2_Q_AC1tTr1t" ); 
  h2_Q_AC1tb      = (TH2F*)ifp1->Get("h2_Q_AC1tb"    );
  h_Q_AC1t        = (TH1F*)ifp1->Get("h_Q_AC1t"      );
  h_Q_AC1b        = (TH1F*)ifp1->Get("h_Q_AC1b"      );

  h2_Q_Trig1t_ev      ->GetXaxis()->SetRangeUser(0,5000);
  h2_Q_Trig1b_ev      ->GetXaxis()->SetRangeUser(0,5000);
  h2_Q_Trig2t_ev      ->GetXaxis()->SetRangeUser(0,5000);
  h2_Q_Trig2b_ev      ->GetXaxis()->SetRangeUser(0,5000);
  h_Q_Trig1t          ->GetXaxis()->SetRangeUser(400,650);
  h_Q_Trig2t          ->GetXaxis()->SetRangeUser(500,700);
  h2_Q_Trig1t2t       ->GetXaxis()->SetRangeUser(400,650);
  h2_Q_Trig1t2t       ->GetYaxis()->SetRangeUser(500,700);
  h_Q_AC1t            ->GetXaxis()->SetRangeUser(550,750);
  h_Q_AC1b            ->GetXaxis()->SetRangeUser(500,700);
  h2_Q_AC1tTr1t       ->GetXaxis()->SetRangeUser(550,750);
  h2_Q_AC1tTr1t       ->GetYaxis()->SetRangeUser(400,650);


  SetTH1(h_Q_Trig1t,"Pedestal Trig1top","QDC [ch]","Counts",1,1000,0);
  SetTH1(h_Q_Trig2t,"Pedestal Trig2top","QDC [ch]","Counts",1,1000,0);
  SetTH1(h_Q_AC1t  ,"Pedestal AC1top"  ,"QDC [ch]","Counts",1,1000,0);
  SetTH2(h2_Q_Trig1t2t  ,"Pedestal Trig1top vs Trig2top","Trig1 top QDC [ch]","Trig2 top QDC [ch]");
  SetTH2(h2_Q_AC1tTr1t  ,"Pedestal AC top vs Trig1top"  ,"AC top QDC [ch]","Trig1 top QDC [ch]");

  TCanvas *c[5];
  for(int i=0;i<5;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1800/2,1800/3);}
  c[0]->Clear();
  c[0]->cd(1);
  gPad->SetLogy(1);h_Q_Trig1t->Draw("");

  c[1]->Clear();
  c[1]->cd(1);
  gPad->SetLogy(1);h_Q_Trig2t->Draw("");

  c[2]->Clear();
  c[2]->cd(1);
  h2_Q_Trig1t2t->Draw("colz");

  c[3]->Clear();
  c[3]->cd(1);
  gPad->SetLogy(1);h_Q_AC1t->Draw("");

  c[4]->Clear();
  c[4]->cd(1);
  h2_Q_AC1tTr1t->Draw("colz");

  c[0] ->Print("pdf/ped_1415.pdf[");
  c[0] ->Print("pdf/ped_1415.pdf");
  c[1] ->Print("pdf/ped_1415.pdf");
  c[2] ->Print("pdf/ped_1415.pdf");
  c[3] ->Print("pdf/ped_1415.pdf");
  c[4] ->Print("pdf/ped_1415.pdf");
  c[4] ->Print("pdf/ped_1415.pdf]");



}

