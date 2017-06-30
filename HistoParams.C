char filename[50];
TFile *fileinput;
char name[100],label[100], labelout[50];
float gain0[21][4][32][6], gain1[21][4][32][6]; // gain for peaks 0-1 and 1-2
Double_t ntotal[21][4][32][6];         // number of total events in histo
Double_t par[9];
float events0;      // number of events in pedestal
float Ncells;       // number of initial fired cells
float DCR;          // dark count rate
float Ndrk;         // dark events, except pedestal peak
float Ndrk_xt;      // dark events, except 0 and 1st peaks 
float Pxt;          // cross-talk probability
float ENFgain;      // ENF for gain
int seq;            // sequence step
Double_t lambda[21][4][32][6];    //lambda of Generalized Poisson
Double_t rms[21][4][32][6], mean[21][4][32][6], mean_over_rms2[21][4][32][6]; //RMS and Mean of histogram
int binwidth[21][4][32][6];
//Double_t mu[21][4][32][6];
Double_t mu[21][4][32][6], mu_hist[21][4][32][6];
int fiber2BV[32][6];
float npe_peaks[21][4][32][6];
int RM, fiber, fiberch, sipm;
TH1F *histo[21][4][32][6];
TCanvas *c[4];
TF1 *func[21][4][32][6];
TF1 *myf[21][4][32][6];
TGraph *enfvsvolt = new TGraph();
TGraph *enfnpevsvolt = new TGraph();
float enf[21][4][32][6], gain[21][4][32][6], gain_mus[21][4][32][6], gain_spe[21][4][32][6], gain0_1[21][4][32][6], enf_mus1[21][4][32][6];
float voltages[20] = {65.9926, 66.1489, 66.3052, 66.4615, 66.6178, 66.7741, 66.9304, 67.0866, 67.2429, 67.3992, 67.5555, 67.7118, 67.8681, 68.0244, 68.1807, 68.3369, 68.4932, 68.6495, 68.8058, 68.9621};
void ReadBVmap();

void ReadBVmap(){ // fiber-fiberch to BVchannel map
  FILE *fileintxt;
  int RM,card,fiber,fiberch,QIE,BV,sipma,sipmx,sipmy;
  sprintf(name,"fiber2BV.txt");
  fileintxt = fopen (name, "rt");
  rewind(fileintxt);
  for (int ifiber=0;ifiber<32;ifiber++){
    for(int ifiberch=0;ifiberch<6;ifiberch++){
      fgets(label,50,fileintxt);
      sscanf(label,"%d\t%d\t%d\t%d\t%d\t%d",&RM,&card,&fiber,&fiberch,&QIE,&BV);
      // printf("%d %d %d %d\n",RM,fiber,fiberch,BV);
      //fiber2BV[(RM-1)*8+fiber+((fiber)%2)-1*((fiber%2)==0)-1][fiberch]=BV-1;
      //fiber2BV[(RM-1)*8+fiber-1][3*(fiberch<3)+fiberch-3*(fiberch>=3)]=BV-1;
      fiber2BV[(RM-1)*8+fiber-1][fiberch]=BV-1;
    }
  }
  fclose(fileintxt);
}


Double_t generpoiss(Double_t *x, Double_t *param){ // definition of the fit function
  Float_t xx = x[0];
  Double_t fitval = 0.;
  for (int k=0; k<npe_peaks[seq][RM-1][fiber][fiberch]+1; k++) {
    fitval += ntotal[seq][RM-1][fiber][fiberch]*binwidth[seq][RM-1][fiber][fiberch]*exp(-0.5*TMath::Power((xx-(par[1]+k*param[1])),2)/(par[2]*par[2]+param[2]*param[2]*k))*(param[3]*TMath::Power((param[3]+k*param[0]),(k-1))*exp(-(param[3]+k*param[0]))/(TMath::Factorial(k)*sqrt(2*3.141593)*sqrt(par[2]*par[2]+param[2]*param[2]*k)));
  }
  return fitval;
}

void HistoParams(int runnumber){
  ReadBVmap();
  for (seq=0; seq<21; seq++){   
    sprintf (name, "analysis_run_%d_RBX_0_seq%d.root", runnumber, seq);
    fileinput = TFile::Open(name);
    printf("file %s opened\n", name);
    sprintf (label, "rbx0_sumTS4/pulse");
    fileinput -> cd(label);
    for (RM = 1; RM<5; RM++){
      sprintf(label, "RM%d", RM);
      gDirectory -> cd(label);
      for (fiber = RM*8-8; fiber<RM*8; fiber++){
	for (fiberch = 0; fiberch<6; fiberch++){
	  sprintf(label, "f%d_fch%d_led=0_RN=%d_sipm=%d", fiber, fiberch, runnumber, fiber2BV[fiber][fiberch]);
	  histo[seq][RM-1][fiber][fiberch] = (TH1F*)gDirectory->Get(label);
	  if (histo[seq][RM-1][fiber][fiberch]->GetEntries() != 0){
	    binwidth[seq][RM-1][fiber][fiberch] = histo[seq][RM-1][fiber][fiberch] -> GetBinWidth(1);
	    int lastbin_x = histo[seq][RM-1][fiber][fiberch] -> GetBinCenter(histo[seq][RM-1][fiber][fiberch]->FindLastBinAbove(2,1));
	    int binmax = histo[seq][RM-1][fiber][fiberch] -> GetMaximumBin();
	    int maximumvalue = histo[seq][RM-1][fiber][fiberch] -> GetBinCenter(binmax);
	    histo[seq][RM-1][fiber][fiberch] -> GetXaxis() -> SetRangeUser(0,maximumvalue-4*(maximumvalue>4));
	    while( (histo[seq][RM-1][fiber][fiberch] -> Integral()) > ((histo[seq][RM-1][fiber][fiberch] -> GetEntries())*0.02))
	      {
		if ((histo[seq][RM-1][fiber][fiberch] -> Integral()) >= (histo[seq][RM-1][fiber][fiberch] -> GetEntries())) break;
		binmax = histo[seq][RM-1][fiber][fiberch] -> GetMaximumBin();
		maximumvalue = histo[seq][RM-1][fiber][fiberch] -> GetBinCenter(binmax);
		if (maximumvalue<5) break;
		histo[seq][RM-1][fiber][fiberch] -> GetXaxis() -> SetRangeUser(0,maximumvalue-4*(maximumvalue>4));
	      }
	    histo[seq][RM-1][fiber][fiberch] -> GetXaxis() -> SetRangeUser(maximumvalue-20,maximumvalue+20);
	    int n0 = histo[seq][RM-1][fiber][fiberch] -> Integral();
	    TF1 *g1 = new TF1("g1","gaus",maximumvalue-20,maximumvalue+20);
	    histo[seq][RM-1][fiber][fiberch] -> Fit(g1,"RQ");
	    g1 -> GetParameters(&par[0]);
	    histo[seq][RM-1][fiber][fiberch] -> GetXaxis() -> SetRangeUser(maximumvalue+20,maximumvalue+100);
	    binmax = histo[seq][RM-1][fiber][fiberch]->GetMaximumBin();
	    maximumvalue = histo[seq][RM-1][fiber][fiberch]->GetBinCenter(binmax);
	    TF1 *g2 = new TF1("g2","gaus",maximumvalue-15,maximumvalue+30);
	    histo[seq][RM-1][fiber][fiberch] -> Fit(g2,"RQ+");
	    g2 -> GetParameters(&par[3]);
	    //Double_t gain0_1 = par[4]-par[1];
	    //Double_t sigma1_pure = sqrt(par[5]*par[5] - par[2]*par[2]);
	    npe_peaks[seq][RM-1][fiber][fiberch] = TMath::Nint((lastbin_x-par[1])/(par[4]-par[1]));
	    histo[seq][RM-1][fiber][fiberch] -> GetXaxis() -> SetRangeUser(0,800);
	    ntotal[seq][RM-1][fiber][fiberch] = histo[seq][RM-1][fiber][fiberch]->GetEntries();
	    mu_hist[seq][RM-1][fiber][fiberch] = -TMath::Log(n0/ntotal[seq][RM-1][fiber][fiberch]);
	    sprintf(label, "generpoiss_fiber%d_ch%d_seq%d",fiber, fiberch, seq);
	    func[seq][RM-1][fiber][fiberch] = new TF1(label, generpoiss, par[1]-15, 800, 4);
	    func[seq][RM-1][fiber][fiberch] -> SetParameters( 0.1, par[4]-par[1], sqrt(abs(par[5]*par[5] - par[2]*par[2])), 2.); //initial parameters
	    func[seq][RM-1][fiber][fiberch] -> SetNpx(1000);
	    //func -> SetParLimits(0, 0.23, 2200.); //set limits (from 0.23 to 2200) for parameter 0
	    //func -> FixParameter(3, 0.23); //fix parameter 3 at 0.23 value
	    func[seq][RM-1][fiber][fiberch] -> SetParNames("lambda", "gain", "sigma1", "mu");
	    histo[seq][RM-1][fiber][fiberch] -> Fit(func[seq][RM-1][fiber][fiberch], "RQ");
	    myf[seq][RM-1][fiber][fiberch]=histo[seq][RM-1][fiber][fiberch] -> GetFunction(label);
	    lambda[seq][RM-1][fiber][fiberch] = myf[seq][RM-1][fiber][fiberch] -> GetParameter(0);
	    rms[seq][RM-1][fiber][fiberch] = histo[seq][RM-1][fiber][fiberch] -> GetRMS();
	    mean[seq][RM-1][fiber][fiberch] = histo[seq][RM-1][fiber][fiberch] -> GetMean();
	    mean_over_rms2[seq][RM-1][fiber][fiberch] = TMath::Power((mean[seq][RM-1][fiber][fiberch]-par[1])/(rms[seq][RM-1][fiber][fiberch]-par[2]),2);
	    mu[seq][RM-1][fiber][fiberch] = myf[seq][RM-1][fiber][fiberch] -> GetParameter(3);

	    //printf("lambda = %f\n", lambda[RM-1][fiber][fiberch]);
	    enf[seq][RM-1][fiber][fiberch] = 1/(1-lambda[seq][RM-1][fiber][fiberch]);
	    enf_mus1[seq][RM-1][fiber][fiberch] = mu[seq][RM-1][fiber][fiberch]*(TMath::Power(rms[seq][RM-1][fiber][fiberch],2)-TMath::Power(par[2],2))/TMath::Power((mean[seq][RM-1][fiber][fiberch]-par[1]),2);
	    //printf("ENF = %f\n", enf[RM-1][fiber][fiberch]);
	    gain[seq][RM-1][fiber][fiberch] = (rms[seq][RM-1][fiber][fiberch]*rms[seq][RM-1][fiber][fiberch]-par[2]*par[2])/((mean[seq][RM-1][fiber][fiberch]-par[1])*enf[seq][RM-1][fiber][fiberch]*enf[seq][RM-1][fiber][fiberch]);
	    gain_mus[seq][RM-1][fiber][fiberch] = (rms[seq][RM-1][fiber][fiberch]*rms[seq][RM-1][fiber][fiberch]-par[2]*par[2])/((mean[seq][RM-1][fiber][fiberch]-par[1])*enf_mus1[seq][RM-1][fiber][fiberch]*enf_mus1[seq][RM-1][fiber][fiberch]);

	    gain0_1[seq][RM-1][fiber][fiberch] = par[4] - par[1];
	    gain_spe[seq][RM-1][fiber][fiberch] = myf[seq][RM-1][fiber][fiberch] -> GetParameter(1);

	    //printf("Gain = %f\n", gain[RM-1][fiber][fiberch]);
	    /*
	    //below is a part concerning histogram saving 
	    sprintf(name, "run%d_generpoiss_fit_seq%d.root", runnumber, seq);
	    bool file_exists = gSystem->AccessPathName(name);
	    //if (file_exists == true){
	    TFile *f = new TFile(name, "update");
	    f -> cd();
	    TDirectory *RMfold;
	    sprintf(label, "RM%d", RM);
	    if (f -> GetDirectory(label) == 0){
	      //TDirectory *RM[4];
	      TDirectory *RMfold = f->mkdir(label);
	      RMfold->cd();
	    }
	    else{
	      TDirectory *RMfold = f->GetDirectory(label);
	      RMfold -> cd();
	    }	    
	    //histo[RM-1][fiber][fiberch] -> GetXaxis() -> SetRangeUser(1, 800);
	    //gStyle -> SetOptStat(1111111);
	    //gStyle -> SetOptFit(111);
	    sprintf(labelout, "f%d_fch%d_led=0_RN=%d_sipm=%d", fiber, fiberch, runnumber, fiber2BV[fiber][fiberch]);
	    bool object_exists = gDirectory->GetListOfKeys()->Contains(labelout);
	    if (object_exists == false){
	      histo[seq][RM-1][fiber][fiberch] -> SetName(labelout);
	      histo[seq][RM-1][fiber][fiberch] -> Write();
	      f -> Close();
	      sprintf (name, "analysis_run_%d_RBX_0_seq%d.root", runnumber, seq);
	      fileinput = TFile::Open(name,"READ");
	      sprintf (label, "rbx0_sumTS4/pulse/RM%d", RM);
	      fileinput -> cd(label);
	    }
	    else{
	      f -> Close();
	      sprintf (name, "analysis_run_%d_RBX_0_seq%d.root", runnumber, seq);
	      fileinput = TFile::Open(name,"READ");
	      sprintf (label, "rbx0_sumTS4/pulse/RM%d", RM);
	      fileinput -> cd(label);
	      }*/
	    //}
	  }
	  
	  else{
	    printf("fiber%d ch%d bv%d\n",fiber, fiberch, fiber2BV[fiber][fiberch]);
	  }
	  
	  if(fiber == RM*8-1 && fiberch == 5 && RM<5)
	    gDirectory->cd("..");
	  
	  if(fiber == 31 && fiberch == 5){
	    fileinput->Close();
	    delete fileinput;
	  }
	}
      }
      //}
    }
    /*
    fileinput -> Close("R");
    delete fileinput;
    printf("file %s finished\n", name);
    */
  }
    
  /*TH1F *gain_distr = new TH1F("gain_distr", "ENF from Npe", 200, 0., 2.);
  for (RM = 1; RM<5; RM++){
    for (fiber = RM*8-8; fiber<RM*8; fiber++){
      for (fiberch = 0; fiberch<6; fiberch++){
	//gain_distr ->Fill(abs(gain0_1[RM-1][fiber][fiberch]-gain[RM-1][fiber][fiberch]));
	if (gain[seq][RM-1][fiber][fiberch] != 0)
	  gain_distr ->Fill((mu[seq][RM-1][fiber][fiberch]*(1-lambda[seq][RM-1][fiber][fiberch]))/mean_over_rms2[seq][RM-1][fiber][fiberch]);
	  //gain_distr ->Fill(gain[RM-1][fiber][fiberch]);
	  //gain_distr ->Fill(enf[RM-1][fiber][fiberch]);///enf_mus1[RM-1][fiber][fiberch]);
	  //gain_distr ->Fill(enf[RM-1][fiber][fiberch]);
	  }
	  }
	  }
	  gain_distr->Draw();*/
  
  for (seq = 1; seq < 21; seq++){
    if ((64.236-voltages[seq-1]) < 0.){
      enfvsvolt -> SetPoint(seq, voltages[seq-1]-64.236, gain[seq][0][0][0]);
      enfnpevsvolt -> SetPoint(seq, voltages[seq-1]-64.236, gain_mus[seq][0][0][0]);
      //enfvsvolt -> SetPoint(seq, voltages[seq-1]-64.510, gain[seq][0][0][0]);
      enfvsvolt -> SetMarkerColor(2);
      enfvsvolt -> SetMarkerStyle(4);
      enfvsvolt -> SetName("ENF-vs-OverVoltage");
      enfvsvolt -> SetTitle("ENF-vs-OverVoltage");
      //enfvsvolt -> GetYaxis() -> SetRangeUser(0,70);
      enfvsvolt -> Draw("ap");
      enfnpevsvolt -> Draw("same");
    }
  }
}
