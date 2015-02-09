#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooAbsReal.h"
#include "RooDataSet.h"

#include "TROOT.h"
#include "TMinuit.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLine.h"
#include "TStopwatch.h"
#include "TRandom3.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;
using namespace RooFit;

Double_t GetLuminosity( Double_t time, Double_t current, Double_t T, TString Target );

void getLimits( Double_t Ebeam = 1.0, TString p0percent = "", Bool_t checkLumi = false, Bool_t fullAcc = false ){

  // Strings for access and output

  TString p0Str = "";
  if( p0percent != "" ){
    p0Str += "_p0pt"+p0percent;
  }

  TString mpPath   = "/Users/jbbeacham/APEX/ReachAnalyzer/APEX/";
  TString filePath = mpPath+"results/";
  //if( pathExt != "" ) filePath += pathExt+"/";
  TString cycle    = "APEX_Cycle.";

  TString EbStr;
  EbStr.Form("%i",int(Ebeam));

  TString outApp = "";
  if( fullAcc ) outApp = "_FullAcc";

  ////////////////////////////////////////////////////
  // Get TTrees:
  // Ten separate runs, between 4.5 and 5.5 degrees,
  // to be combined
  ////////////////////////////////////////////////////

  TFile* f_BH_HRS0 = TFile::Open(filePath+cycle+"MC.BH_Setting45to150mrad"+EbStr+"GeVInvMassGreq50MeV_No1_Rot2_FullRun_nom5degHRS0"+p0Str+outApp+".root","READ");
  TFile* f_Rad_HRS0 = TFile::Open(filePath+cycle+"MC.Rad_Setting45to150mrad"+EbStr+"GeVInvMassGreq50MeV_No1_Rot2_FullRun_nom5degHRS0"+p0Str+outApp+".root","READ");

  TFile* f_BH_HRS1 = TFile::Open(filePath+cycle+"MC.BH_Setting45to150mrad"+EbStr+"GeVInvMassGreq50MeV_No1_Rot2_FullRun_nom5degHRS1"+p0Str+outApp+".root","READ");
  TFile* f_Rad_HRS1 = TFile::Open(filePath+cycle+"MC.Rad_Setting45to150mrad"+EbStr+"GeVInvMassGreq50MeV_No1_Rot2_FullRun_nom5degHRS1"+p0Str+outApp+".root","READ");

  TFile* f_BH_HRS2 = TFile::Open(filePath+cycle+"MC.BH_Setting45to150mrad"+EbStr+"GeVInvMassGreq50MeV_No1_Rot2_FullRun_nom5degHRS2"+p0Str+outApp+".root","READ");
  TFile* f_Rad_HRS2 = TFile::Open(filePath+cycle+"MC.Rad_Setting45to150mrad"+EbStr+"GeVInvMassGreq50MeV_No1_Rot2_FullRun_nom5degHRS2"+p0Str+outApp+".root","READ");

  TFile* f_BH_HRS3 = TFile::Open(filePath+cycle+"MC.BH_Setting45to150mrad"+EbStr+"GeVInvMassGreq50MeV_No1_Rot2_FullRun_nom5degHRS3"+p0Str+outApp+".root","READ");
  TFile* f_Rad_HRS3 = TFile::Open(filePath+cycle+"MC.Rad_Setting45to150mrad"+EbStr+"GeVInvMassGreq50MeV_No1_Rot2_FullRun_nom5degHRS3"+p0Str+outApp+".root","READ");

  TFile* f_BH_HRS4 = TFile::Open(filePath+cycle+"MC.BH_Setting45to150mrad"+EbStr+"GeVInvMassGreq50MeV_No1_Rot2_FullRun_nom5degHRS4"+p0Str+outApp+".root","READ");
  TFile* f_Rad_HRS4 = TFile::Open(filePath+cycle+"MC.Rad_Setting45to150mrad"+EbStr+"GeVInvMassGreq50MeV_No1_Rot2_FullRun_nom5degHRS4"+p0Str+outApp+".root","READ");

  TFile* f_BH_HRS5 = TFile::Open(filePath+cycle+"MC.BH_Setting45to150mrad"+EbStr+"GeVInvMassGreq50MeV_No1_Rot2_FullRun_nom5degHRS5"+p0Str+outApp+".root","READ");
  TFile* f_Rad_HRS5 = TFile::Open(filePath+cycle+"MC.Rad_Setting45to150mrad"+EbStr+"GeVInvMassGreq50MeV_No1_Rot2_FullRun_nom5degHRS5"+p0Str+outApp+".root","READ");

  TFile* f_BH_HRS6 = TFile::Open(filePath+cycle+"MC.BH_Setting45to150mrad"+EbStr+"GeVInvMassGreq50MeV_No1_Rot2_FullRun_nom5degHRS6"+p0Str+outApp+".root","READ");
  TFile* f_Rad_HRS6 = TFile::Open(filePath+cycle+"MC.Rad_Setting45to150mrad"+EbStr+"GeVInvMassGreq50MeV_No1_Rot2_FullRun_nom5degHRS6"+p0Str+outApp+".root","READ");

  TFile* f_BH_HRS7 = TFile::Open(filePath+cycle+"MC.BH_Setting45to150mrad"+EbStr+"GeVInvMassGreq50MeV_No1_Rot2_FullRun_nom5degHRS7"+p0Str+outApp+".root","READ");
  TFile* f_Rad_HRS7 = TFile::Open(filePath+cycle+"MC.Rad_Setting45to150mrad"+EbStr+"GeVInvMassGreq50MeV_No1_Rot2_FullRun_nom5degHRS7"+p0Str+outApp+".root","READ");

  TFile* f_BH_HRS8 = TFile::Open(filePath+cycle+"MC.BH_Setting45to150mrad"+EbStr+"GeVInvMassGreq50MeV_No1_Rot2_FullRun_nom5degHRS8"+p0Str+outApp+".root","READ");
  TFile* f_Rad_HRS8 = TFile::Open(filePath+cycle+"MC.Rad_Setting45to150mrad"+EbStr+"GeVInvMassGreq50MeV_No1_Rot2_FullRun_nom5degHRS8"+p0Str+outApp+".root","READ");

  TFile* f_BH_HRS9 = TFile::Open(filePath+cycle+"MC.BH_Setting45to150mrad"+EbStr+"GeVInvMassGreq50MeV_No1_Rot2_FullRun_nom5degHRS9"+p0Str+outApp+".root","READ");
  TFile* f_Rad_HRS9 = TFile::Open(filePath+cycle+"MC.Rad_Setting45to150mrad"+EbStr+"GeVInvMassGreq50MeV_No1_Rot2_FullRun_nom5degHRS9"+p0Str+outApp+".root","READ");


  TTree* tree_BH_HRS0 = (TTree*) f_BH_HRS0->Get("invmass");
  TTree* tree_Rad_HRS0 = (TTree*) f_Rad_HRS0->Get("invmass");

  TTree* tree_BH_HRS1 = (TTree*) f_BH_HRS1->Get("invmass");
  TTree* tree_Rad_HRS1 = (TTree*) f_Rad_HRS1->Get("invmass");

  TTree* tree_BH_HRS2 = (TTree*) f_BH_HRS2->Get("invmass");
  TTree* tree_Rad_HRS2 = (TTree*) f_Rad_HRS2->Get("invmass");

  TTree* tree_BH_HRS3 = (TTree*) f_BH_HRS3->Get("invmass");
  TTree* tree_Rad_HRS3 = (TTree*) f_Rad_HRS3->Get("invmass");

  TTree* tree_BH_HRS4 = (TTree*) f_BH_HRS4->Get("invmass");
  TTree* tree_Rad_HRS4 = (TTree*) f_Rad_HRS4->Get("invmass");

  TTree* tree_BH_HRS5 = (TTree*) f_BH_HRS5->Get("invmass");
  TTree* tree_Rad_HRS5 = (TTree*) f_Rad_HRS5->Get("invmass");

  TTree* tree_BH_HRS6 = (TTree*) f_BH_HRS6->Get("invmass");
  TTree* tree_Rad_HRS6 = (TTree*) f_Rad_HRS6->Get("invmass");

  TTree* tree_BH_HRS7 = (TTree*) f_BH_HRS7->Get("invmass");
  TTree* tree_Rad_HRS7 = (TTree*) f_Rad_HRS7->Get("invmass");

  TTree* tree_BH_HRS8 = (TTree*) f_BH_HRS8->Get("invmass");
  TTree* tree_Rad_HRS8 = (TTree*) f_Rad_HRS8->Get("invmass");

  TTree* tree_BH_HRS9 = (TTree*) f_BH_HRS9->Get("invmass");
  TTree* tree_Rad_HRS9 = (TTree*) f_Rad_HRS9->Get("invmass");


  ///////////////////////////////////////
  // Cross-section and N_generated info
  ///////////////////////////////////////

  Double_t xsec_BH_gen  = 254272.0;
  Double_t xsec_Rad_gen = 95399.2;

  Int_t N_BH_gen  = 200000;
  Int_t N_Rad_gen = 188310;

  if( Ebeam > 1.056 ){
    xsec_BH_gen  = 58850.0;
    xsec_Rad_gen = 18815.1;
    N_Rad_gen = 142191;
  }
  if( Ebeam > 2.056 ){
    xsec_BH_gen  = 21006.8;
    xsec_Rad_gen = 6853.25;
    N_Rad_gen = 110093;
  }
  if( Ebeam > 3.056 ){
    xsec_BH_gen  = 9072.12;
    xsec_Rad_gen = 2759.19;
    N_BH_gen = 182577;
    N_Rad_gen = 132561;
  }

  // Calculate equivalent luminosity

  Double_t lumi_eq_BH  = double(N_BH_gen)/xsec_BH_gen;
  Double_t lumi_eq_Rad = double(N_Rad_gen)/xsec_Rad_gen;

  // Get total luminosity
  // Values used in proposal:

  // 1-pass: 50 uA, 0.7% X0, 6 days
  // ...or...
  // 1-pass: 65 uA, 0.7% X0, 6 days
  // 2-pass: 70 uA, 4% X0, 6 days
  // 3-pass: 80 uA, 8% X0, 6 days
  // 4-pass: 60 uA, 8% X0, 12 days

  Double_t time    = 6.0;   // Days
  Double_t current = 65.0;  // micro amp
  Double_t T       = 0.007;   // Target thickness, in radiation lengths
  TString Target   = "W";   // Target nucleus

  // Adjustment for 10-settings: 1/10 of target thickness each, combined
  // Ten degree settings

  Double_t thetaIn_4p5 = 0.0785;  // Opening angle setting, 4p5degHRS
  Double_t thetaIn_5p5 = 0.09599; // Opening angle setting, 5p5degHRS

  Double_t thetasIn[10];
  for( Int_t ii = 0; ii < 10; ++ii ){
    thetasIn[ii] = thetaIn_4p5+ii*(thetaIn_5p5-thetaIn_4p5)/9.0;
  }

  // Adjust for different beam energies
  if( Ebeam > 1.056 ){
    time = 6.0;
    current = 70.0;
    T = 0.04;
  }
  if( Ebeam > 2.056 ){
    time = 6.0;
    current = 80.0;
    T = 0.08;
  }
  if( Ebeam > 3.056 ){
    time = 12.0;
    current = 60.0;
    T = 0.08;
  }

  // Target thickness per foil, for 10 foils;
  // T/10 for 10-setting combination.

  Double_t T10 = T/10.0;
  Double_t lumi_tot = GetLuminosity(time,current,T10,Target);

  std::cout << "    lumi_tot, 1/10 target thickness  = " << lumi_tot << std::endl;
  std::cout << "                        weight, BH  = " << lumi_tot/lumi_eq_BH << std::endl;
  std::cout << "                        weight, Rad = " << lumi_tot/lumi_eq_Rad << std::endl;

  TString lumiTempFile = "LumiTempFile.txt";

  ofstream outFX;
  outFX.open(lumiTempFile,ios::app);
  outFX.setf(ios::fixed,ios::floatfield);

  if( checkLumi ){
    outFX << "    lumi_tot, 1/10 target thickness  = " << lumi_tot << "\n";
    outFX << "                        weight, BH  = " << lumi_tot/lumi_eq_BH << "\n";
    outFX << "                        weight, Rad = " << lumi_tot/lumi_eq_Rad << "\n";

    return;
  }

  /////////////////////////////////
  // Construct epsilon^2 function
  /////////////////////////////////

  // Get mA values from a histogram (for possible backwards compatibility)

  TString histName = "h_eplus_eminus_invmass_bin1MeV";
  //TString histName = "h_eplus_eminus_invmass";
  TH1F* invMassBH  = (TH1F*) f_BH_HRS0->Get( histName );

  Double_t inLo = invMassBH->GetXaxis()->GetBinLowEdge(1);
  Double_t inHi = invMassBH->GetXaxis()->GetBinUpEdge(invMassBH->GetNbinsX());
  Double_t mApt = invMassBH->GetXaxis()->GetBinCenter(1);

  RooRealVar* mass   = new RooRealVar("mass","mass",0.0,1000.0);

  RooDataSet* data_BH_HRS0  = new RooDataSet("data_BH_HRS0", "data_BH_HRS0", tree_BH_HRS0, RooArgSet(*mass));
  RooDataSet* data_Rad_HRS0 = new RooDataSet("data_Rad_HRS0", "data_Rad_HRS0", tree_Rad_HRS0, RooArgSet(*mass));

  RooDataSet* data_BH_HRS1  = new RooDataSet("data_BH_HRS1", "data_BH_HRS1", tree_BH_HRS1, RooArgSet(*mass));
  RooDataSet* data_Rad_HRS1 = new RooDataSet("data_Rad_HRS1", "data_Rad_HRS1", tree_Rad_HRS1, RooArgSet(*mass));

  RooDataSet* data_BH_HRS2  = new RooDataSet("data_BH_HRS2", "data_BH_HRS2", tree_BH_HRS2, RooArgSet(*mass));
  RooDataSet* data_Rad_HRS2 = new RooDataSet("data_Rad_HRS2", "data_Rad_HRS2", tree_Rad_HRS2, RooArgSet(*mass));

  RooDataSet* data_BH_HRS3  = new RooDataSet("data_BH_HRS3", "data_BH_HRS3", tree_BH_HRS3, RooArgSet(*mass));
  RooDataSet* data_Rad_HRS3 = new RooDataSet("data_Rad_HRS3", "data_Rad_HRS3", tree_Rad_HRS3, RooArgSet(*mass));

  RooDataSet* data_BH_HRS4  = new RooDataSet("data_BH_HRS4", "data_BH_HRS4", tree_BH_HRS4, RooArgSet(*mass));
  RooDataSet* data_Rad_HRS4 = new RooDataSet("data_Rad_HRS4", "data_Rad_HRS4", tree_Rad_HRS4, RooArgSet(*mass));

  RooDataSet* data_BH_HRS5  = new RooDataSet("data_BH_HRS5", "data_BH_HRS5", tree_BH_HRS5, RooArgSet(*mass));
  RooDataSet* data_Rad_HRS5 = new RooDataSet("data_Rad_HRS5", "data_Rad_HRS5", tree_Rad_HRS5, RooArgSet(*mass));

  RooDataSet* data_BH_HRS6  = new RooDataSet("data_BH_HRS6", "data_BH_HRS6", tree_BH_HRS6, RooArgSet(*mass));
  RooDataSet* data_Rad_HRS6 = new RooDataSet("data_Rad_HRS6", "data_Rad_HRS6", tree_Rad_HRS6, RooArgSet(*mass));

  RooDataSet* data_BH_HRS7  = new RooDataSet("data_BH_HRS7", "data_BH_HRS7", tree_BH_HRS7, RooArgSet(*mass));
  RooDataSet* data_Rad_HRS7 = new RooDataSet("data_Rad_HRS7", "data_Rad_HRS7", tree_Rad_HRS7, RooArgSet(*mass));

  RooDataSet* data_BH_HRS8  = new RooDataSet("data_BH_HRS8", "data_BH_HRS8", tree_BH_HRS8, RooArgSet(*mass));
  RooDataSet* data_Rad_HRS8 = new RooDataSet("data_Rad_HRS8", "data_Rad_HRS8", tree_Rad_HRS8, RooArgSet(*mass));

  RooDataSet* data_BH_HRS9  = new RooDataSet("data_BH_HRS9", "data_BH_HRS9", tree_BH_HRS9, RooArgSet(*mass));
  RooDataSet* data_Rad_HRS9 = new RooDataSet("data_Rad_HRS9", "data_Rad_HRS9", tree_Rad_HRS9, RooArgSet(*mass));

  vector< RooDataSet* > dataSets_BH;
  dataSets_BH.push_back(data_BH_HRS0); dataSets_BH.push_back(data_BH_HRS1);
  dataSets_BH.push_back(data_BH_HRS2); dataSets_BH.push_back(data_BH_HRS3);
  dataSets_BH.push_back(data_BH_HRS4); dataSets_BH.push_back(data_BH_HRS5);
  dataSets_BH.push_back(data_BH_HRS6); dataSets_BH.push_back(data_BH_HRS7);
  dataSets_BH.push_back(data_BH_HRS8); dataSets_BH.push_back(data_BH_HRS9);

  vector< RooDataSet* > dataSets_Rad;
  dataSets_Rad.push_back(data_Rad_HRS0); dataSets_Rad.push_back(data_Rad_HRS1);
  dataSets_Rad.push_back(data_Rad_HRS2); dataSets_Rad.push_back(data_Rad_HRS3);
  dataSets_Rad.push_back(data_Rad_HRS4); dataSets_Rad.push_back(data_Rad_HRS5);
  dataSets_Rad.push_back(data_Rad_HRS6); dataSets_Rad.push_back(data_Rad_HRS7);
  dataSets_Rad.push_back(data_Rad_HRS8); dataSets_Rad.push_back(data_Rad_HRS9);

  // Correction factors for full interference diagrams, from Natalia

  // Updated, from Natalia:
  // 1-pass: 3.6
  // 2-pass: 3.4
  // 3-pass: 2.95
  // 4-pass: 2.7

  Double_t CorrFact = 3.1;
  if( Ebeam > 1.056 ) CorrFact = 2.8;
  if( Ebeam > 2.056 ) CorrFact = 2.5;
  if( Ebeam > 3.056 ) CorrFact = 2.2;

  // Finally, build the epsilon^2 function

  RooRealVar* preF   = new RooRealVar("preF","preF",CorrFact);
  RooRealVar* B_BH   = new RooRealVar("B_BH","B_BH",100,0,10000000);
  RooRealVar* B_Rad  = new RooRealVar("B_Rad","B_Rad",100,0,10000000);
  RooAbsReal* Neff   = new RooRealVar("Neff","Neff",1.0);
  RooRealVar* alpha  = new RooRealVar("alpha","alpha",1.0/137.036);
  RooRealVar* mA     = new RooRealVar("mA","mA",mApt,inLo,inHi);

  RooRealVar* EbeamF  = new RooRealVar("EbeamF","EbeamF",Ebeam);

  RooRealVar* TfMS    = new RooRealVar("TfMS","TfMS",T10/2.0);
  RooRealVar* theta   = new RooRealVar("theta","theta",thetaIn_4p5); // Initialize to 4.5 deg

  RooAbsReal* dThetaMS = new RooFormulaVar("dThetaMS","(13.8/((@0/2.0)*1000.0))*sqrt(@1)*(1.0+0.038*log(@1))",
					   RooArgList(*EbeamF,*TfMS));

  /*
  RooAbsReal* dThetaMS = new RooFormulaVar("dThetaMS","(13.8/((@0/2.0)*1000.0))*sqrt(@1)*(1.0+0.038*log(@1))",
					   RooArgList(*EbeamF,*Teff));
  */

  RooAbsReal* deltam   = new RooFormulaVar("deltam","2.5*@0*sqrt( (3.0*10^(-4))^(2) + 0.5*((0.5*10^(-3)/@1)^2) + 0.5*((@2/@1)^2) )",
					   RooArgList(*mA,*theta,*dThetaMS));


  RooAbsReal* eps2 = new RooFormulaVar("eps2",
				       "( (2.0*(sqrt(@0*(@1+@2)))/(0.8*@2))*(@3/@4)*(2.0*@5*@6/(3.0*pi)) )",
				       RooArgList(*preF,*B_BH,*B_Rad,*deltam,*mA,*Neff,*alpha) );


  ///////////////
  // Get limits
  ///////////////

  Int_t nBins = invMassBH->GetNbinsX();

  Double_t eps2List[nBins]; // To hold the combined values of eps^2, for a given mA

  Double_t mAList[nBins];

  // For plots of dm(mA)

  Int_t dmMod = 8;
  Int_t nBins_dm = int(nBins/dmMod);

  Double_t mAList_dm[nBins_dm];

  Double_t dm_HRS0[nBins_dm];
  Double_t dm_HRS1[nBins_dm];
  Double_t dm_HRS2[nBins_dm];
  Double_t dm_HRS3[nBins_dm];
  Double_t dm_HRS4[nBins_dm];
  Double_t dm_HRS5[nBins_dm];
  Double_t dm_HRS6[nBins_dm];
  Double_t dm_HRS7[nBins_dm];
  Double_t dm_HRS8[nBins_dm];
  Double_t dm_HRS9[nBins_dm];

  vector< Double_t* > dm_All;
  dm_All.push_back(dm_HRS0); dm_All.push_back(dm_HRS1);
  dm_All.push_back(dm_HRS2); dm_All.push_back(dm_HRS3);
  dm_All.push_back(dm_HRS4); dm_All.push_back(dm_HRS5);
  dm_All.push_back(dm_HRS6); dm_All.push_back(dm_HRS7);
  dm_All.push_back(dm_HRS8); dm_All.push_back(dm_HRS9);

  TH1D* nExp_HRS0_BH = 0;
  TH1D* nExp_HRS1_BH = 0;
  TH1D* nExp_HRS2_BH = 0;
  TH1D* nExp_HRS3_BH = 0;
  TH1D* nExp_HRS4_BH = 0;
  TH1D* nExp_HRS5_BH = 0;
  TH1D* nExp_HRS6_BH = 0;
  TH1D* nExp_HRS7_BH = 0;
  TH1D* nExp_HRS8_BH = 0;
  TH1D* nExp_HRS9_BH = 0;

  TH1D* nExp_HRS0_Rad = 0;
  TH1D* nExp_HRS1_Rad = 0;
  TH1D* nExp_HRS2_Rad = 0;
  TH1D* nExp_HRS3_Rad = 0;
  TH1D* nExp_HRS4_Rad = 0;
  TH1D* nExp_HRS5_Rad = 0;
  TH1D* nExp_HRS6_Rad = 0;
  TH1D* nExp_HRS7_Rad = 0;
  TH1D* nExp_HRS8_Rad = 0;
  TH1D* nExp_HRS9_Rad = 0;

  vector< TH1D* > vec_BH;
  vec_BH.push_back(nExp_HRS0_BH); vec_BH.push_back(nExp_HRS1_BH);
  vec_BH.push_back(nExp_HRS2_BH); vec_BH.push_back(nExp_HRS3_BH);
  vec_BH.push_back(nExp_HRS4_BH); vec_BH.push_back(nExp_HRS5_BH);
  vec_BH.push_back(nExp_HRS6_BH); vec_BH.push_back(nExp_HRS7_BH);
  vec_BH.push_back(nExp_HRS8_BH); vec_BH.push_back(nExp_HRS9_BH);

  vector< TH1D* > vec_Rad;
  vec_Rad.push_back(nExp_HRS0_Rad); vec_Rad.push_back(nExp_HRS1_Rad);
  vec_Rad.push_back(nExp_HRS2_Rad); vec_Rad.push_back(nExp_HRS3_Rad);
  vec_Rad.push_back(nExp_HRS4_Rad); vec_Rad.push_back(nExp_HRS5_Rad);
  vec_Rad.push_back(nExp_HRS6_Rad); vec_Rad.push_back(nExp_HRS7_Rad);
  vec_Rad.push_back(nExp_HRS8_Rad); vec_Rad.push_back(nExp_HRS9_Rad);

  for( UInt_t ii = 0; ii < 10; ++ii ){
    TString HRSStr;
    HRSStr.Form("%i",ii);
    vec_BH[ii] = new TH1D("nExp_HRS"+HRSStr+"_BH","nExp_HRS"+HRSStr+"_BH",nBins,0,500);
    vec_Rad[ii] = new TH1D("nExp_HRS"+HRSStr+"_Rad","nExp_HRS"+HRSStr+"_Rad",nBins,0,500);
  }

  Int_t nSet = 10;
  Double_t eps2_vals[nSet];

  Double_t Nexp_win_BH = 0.0;
  Double_t Nexp_win_Rad = 0.0;

  Int_t dmInd = 0;

  // Loop over mass values

  for( Int_t ii = 0; ii < nBins; ++ii ){

    for( Int_t j = 0; j < nSet; ++j ) eps2_vals[j] = 0.0;
    Nexp_win_BH = 0.0;
    Nexp_win_Rad = 0.0;

    Double_t thismA = invMassBH->GetXaxis()->GetBinCenter(ii);

    mA->setVal( thismA );

    // Get eps2 contributions from each of the ten HRS settings:

    Double_t dm = 0.0;

    for( Int_t hh = 0; hh < nSet; ++hh ){

      theta->setVal(thetasIn[hh]);
      
      // Get mass resolution, to get range on data
      
      dm = deltam->getVal();

      TString whichOfTen;
      whichOfTen.Form("%i",hh);

      TString windowStr = "windowHRS"+whichOfTen;
      mass->setRange(windowStr,thismA-(dm/2.0),thismA+(dm/2.0));

      if( ii%dmMod == 0 ){
	dm_All[hh][dmInd] = dm;
	mAList_dm[dmInd] = thismA;
      }

      // Now get Nexp
      
      Nexp_win_BH = (lumi_tot/lumi_eq_BH)*(dataSets_BH[hh])->sumEntries(0,windowStr);
      Nexp_win_Rad = (lumi_tot/lumi_eq_Rad)*(dataSets_Rad[hh])->sumEntries(0,windowStr);
      
      vec_BH[hh]->Fill( thismA, Nexp_win_BH );
      vec_Rad[hh]->Fill( thismA, Nexp_win_Rad );
      
      // Now get eps2
      
      B_BH->setVal( Nexp_win_BH );
      B_Rad->setVal( Nexp_win_Rad );
      
      eps2_vals[hh] = eps2->getVal();
      
      dm = 0.0;
      Nexp_win_BH = 0.0;
      Nexp_win_Rad = 0.0;
      
    }

    ////////////////////////////
    // Now combine eps2 values
    ////////////////////////////

    Double_t oneOverEps2 = 0.0;

    for( Int_t hh = 0; hh < nSet; ++hh ){

      if( eps2_vals[hh] < 1e-10 ){
	continue;
      }else{
	oneOverEps2 += 1.0/(eps2_vals[hh]*eps2_vals[hh]);
      }

    }

    Double_t eps2_combined  = 1.0/(sqrt( oneOverEps2 ));

    ///////////////////////
    // Fill output arrays
    ///////////////////////

    mAList[ii] = thismA;
    eps2List[ii] = eps2_combined;

    std::cout << "mA = " << mAList[ii] << "     eps2 = " << eps2List[ii] << std::endl;
    //std::cout << "Nexp_win_BH  = " << Nexp_win_BH_tot << std::endl;
    //std::cout << "Nexp_win_Rad = " << Nexp_win_Rad_tot << std::endl;

    Nexp_win_BH  = 0.0;
    Nexp_win_Rad = 0.0;

    if( ii%dmMod == 0 ) ++dmInd;
  }

  /////////////////////
  // Make plots: eps2
  /////////////////////

  gROOT->SetStyle("Plain");

  Double_t XAlign2 = 0.14;
  Double_t YAlign2 = 0.84;

  TLatex* apexText = new TLatex(XAlign2,YAlign2,"#it{APEX} Internal");
  apexText->SetTextSize(0.038);
  apexText->SetNDC(kTRUE);


  TCanvas* epsExCanv = new TCanvas("epsEx","epsEx",2);
  TH1F* epsExFr = epsExCanv->DrawFrame(inLo,0.00000001,inHi,0.0001);

  epsExCanv->SetLogy();
  epsExFr->SetTitle("");
  epsExFr->SetXTitle("m_{A'} (MeV/c^{2})");
  epsExFr->SetYTitle("#varepsilon^{2}");
  epsExFr->SetTitleSize(0.04,"X");
  epsExFr->SetTitleSize(0.06,"Y");
  epsExFr->SetLabelSize(0.04,"X");
  epsExFr->SetLabelSize(0.035,"Y");
  epsExFr->SetTitleOffset(0.75,"Y");

  TGraph* greps2 = new TGraph(nBins,mAList,eps2List);
  greps2->Draw("PL");
  greps2->SetMarkerColor(kRed);
  greps2->SetLineColor(kBlue);
  greps2->SetMarkerStyle(20);
  greps2->SetMarkerSize(0.6);

  // Text on plot

  Double_t XAlign = 0.55;
  Double_t YAlign = 0.85;
  Double_t textSize = 0.03;

  TLatex* EbeamText = new TLatex(XAlign,YAlign,EbStr+" pass");
  EbeamText->SetTextSize(textSize);
  EbeamText->SetNDC(kTRUE);

  //TString dmStr;
  //dmStr.Form("%1.1f",massRes);

  TLatex* dmText = new TLatex(XAlign,YAlign-0.05,"Variable #delta m");
  dmText->SetTextSize(textSize);
  dmText->SetNDC(kTRUE);

  TString degStr = "nom5degHRS";

  TLatex* degText = new TLatex(XAlign,YAlign-0.1,degStr);
  degText->SetTextSize(textSize);
  degText->SetNDC(kTRUE);

  TString timeStr;
  timeStr.Form("%i",int(time));

  TString IStr;
  IStr.Form("%2.1f",current);

  TLatex* tIText = new TLatex(XAlign,YAlign-0.15,timeStr+" days at "+IStr+" #muA");
  tIText->SetTextSize(textSize);
  tIText->SetNDC(kTRUE);

  TString TStr;
  TStr.Form("%1.3f",T);

  TLatex* tTText = new TLatex(XAlign,YAlign-0.2,TStr+" rad lengths");
  tTText->SetTextSize(textSize);
  tTText->SetNDC(kTRUE);

  TLatex* tarText = new TLatex(XAlign,YAlign-0.25,Target+" target");
  tarText->SetTextSize(textSize);
  tarText->SetNDC(kTRUE);

  TLatex* p0Text = 0;
  if( p0percent != "" ){
    p0Text = new TLatex(XAlign-0.4,YAlign-0.05,"#frac{p_{0,in}}{p_{0,nom}} = 0."+p0percent);
  }else{
    p0Text = new TLatex(XAlign-0.4,YAlign-0.05,"#frac{p_{0,in}}{p_{0,nom}} = 1");
  }
  p0Text->SetNDC(kTRUE);
  p0Text->SetTextSize(textSize);

  EbeamText->Draw("same");
  dmText->Draw("same");
  degText->Draw("same");
  tIText->Draw("same");
  tTText->Draw("same");
  tarText->Draw("same");
  p0Text->Draw("same");
  apexText->Draw("same");

  greps2->SaveAs("EpsilonGraphs/eps2_"+EbStr+"GeV_"+degStr+"_ten_HRS_settings_combined_variable_dm_v2_"+timeStr+"days_"+IStr+"uA_target"+Target+outApp+".root");

  /*
  //////////////////////////////////////////////////////////
  // Indvidiual eps2 contributions for each degree setting
  //////////////////////////////////////////////////////////

  TCanvas* epsExCanvAlt = new TCanvas("epsExCanvAlt","epsExCanvAlt",2);
  TH1F* epsExFrAlt = epsExCanvAlt->DrawFrame(inLo,0.00000001,inHi,0.0001);

  epsExCanvAlt->SetLogy();
  epsExFrAlt->SetTitle("");
  epsExFrAlt->SetXTitle("m_{A'} (MeV/c^{2})");
  epsExFrAlt->SetYTitle("#varepsilon^{2}");
  epsExFrAlt->SetTitleSize(0.04,"X");
  epsExFrAlt->SetTitleSize(0.06,"Y");
  epsExFrAlt->SetLabelSize(0.04,"X");
  epsExFrAlt->SetLabelSize(0.035,"Y");
  epsExFrAlt->SetTitleOffset(0.75,"Y");

  TGraph* greps2_4p5 = new TGraph(nBins,mAList,eps2List_4p5);
  greps2_4p5->Draw("PL");
  greps2_4p5->SetMarkerColor(kCyan);
  greps2_4p5->SetLineColor(kBlue);
  greps2_4p5->SetMarkerStyle(20);
  greps2_4p5->SetMarkerSize(0.6);

  TGraph* greps2_5 = new TGraph(nBins,mAList,eps2List_5);
  greps2_5->Draw("PL");
  greps2_5->SetMarkerColor(kMagenta);
  greps2_5->SetLineColor(kBlue);
  greps2_5->SetMarkerStyle(20);
  greps2_5->SetMarkerSize(0.6);

  TGraph* greps2_5p5 = new TGraph(nBins,mAList,eps2List_5p5);
  greps2_5p5->Draw("PL");
  greps2_5p5->SetMarkerColor(kGreen);
  greps2_5p5->SetLineColor(kBlue);
  greps2_5p5->SetMarkerStyle(20);
  greps2_5p5->SetMarkerSize(0.6);

  greps2->Draw("PL");

  TLegend* legEps = new TLegend(XAlign2,0.62,XAlign2+0.2,0.8);
  legEps->SetFillColor(0);
  legEps->SetTextSize(0.022);
  legEps->SetBorderSize(0);
  legEps->AddEntry(greps2_4p5,EbStr+" pass; 4.5#circHRS","lp");
  legEps->AddEntry(greps2_5,EbStr+" pass; 5#circHRS","lp");
  legEps->AddEntry(greps2_5p5,EbStr+" pass; 5.5#circHRS","lp");
  legEps->AddEntry(greps2,EbStr+" pass; combined","lp");
  legEps->Draw();

  EbeamText->Draw("same");
  dmText->Draw("same");
  degText->Draw("same");
  tIText->Draw("same");
  tTText->Draw("same");
  tarText->Draw("same");
  //p0Text->Draw("same");
  apexText->Draw("same");

  greps2_4p5->SaveAs("EpsilonGraphs/eps2_individual_contributions_"+EbStr+"GeV_4p5degHRS_variable_dm_v2_"+timeStr+"days_"+IStr+"uA_target"+Target+".root");
  greps2_5->SaveAs("EpsilonGraphs/eps2_individual_contributions_"+EbStr+"GeV_5degHRS_variable_dm_v2_"+timeStr+"days_"+IStr+"uA_target"+Target+".root");
  greps2_5p5->SaveAs("EpsilonGraphs/eps2_individual_contributions_"+EbStr+"GeV_5p5degHRS_variable_dm_v2_"+timeStr+"days_"+IStr+"uA_target"+Target+".root");

  epsExCanvAlt->SaveAs("Plots_Epsilon/IndividualDegreeSettings/eps2_individual_contributions_"+EbStr+"GeV_variable_dm_v2_"+timeStr+"days_"+IStr+"uA_target"+Target+".pdf");
  */
  ///////////////////////////
  // Mass resolution graphs
  ///////////////////////////

  TGraph* grdm_HRS0 = 0;
  TGraph* grdm_HRS1 = 0;
  TGraph* grdm_HRS2 = 0;
  TGraph* grdm_HRS3 = 0;
  TGraph* grdm_HRS4 = 0;
  TGraph* grdm_HRS5 = 0;
  TGraph* grdm_HRS6 = 0;
  TGraph* grdm_HRS7 = 0;
  TGraph* grdm_HRS8 = 0;
  TGraph* grdm_HRS9 = 0;

  vector< TGraph* > graphs_dm;
  graphs_dm.push_back(grdm_HRS0); graphs_dm.push_back(grdm_HRS1);
  graphs_dm.push_back(grdm_HRS2); graphs_dm.push_back(grdm_HRS3);
  graphs_dm.push_back(grdm_HRS4); graphs_dm.push_back(grdm_HRS5);
  graphs_dm.push_back(grdm_HRS6); graphs_dm.push_back(grdm_HRS7);
  graphs_dm.push_back(grdm_HRS8); graphs_dm.push_back(grdm_HRS9);

  for( UInt_t ii = 0; ii < graphs_dm.size(); ++ii ){
    graphs_dm[ii] = new TGraph(nBins_dm,mAList_dm,dm_All[ii]);
    graphs_dm[ii]->SetLineColor(kBlue);
    graphs_dm[ii]->SetMarkerSize(0.4);
    graphs_dm[ii]->SetMarkerStyle(20);
  }

  graphs_dm[0]->SetMarkerColor(kBlue);
  graphs_dm[1]->SetMarkerColor(kMagenta);
  graphs_dm[2]->SetMarkerColor(kGreen);
  graphs_dm[3]->SetMarkerColor(kYellow+2);
  graphs_dm[4]->SetMarkerColor(kOrange);
  graphs_dm[5]->SetMarkerColor(kCyan);
  graphs_dm[6]->SetMarkerColor(kAzure);
  graphs_dm[7]->SetMarkerColor(kRed);
  graphs_dm[8]->SetMarkerColor(kPink);
  graphs_dm[9]->SetMarkerColor(kBlack);

  TCanvas* dmCanv = new TCanvas("dmCanv","dmCanv",2);
  TH1F* dmFr = dmCanv->DrawFrame(inLo,0.1,inHi,7);

  dmFr->SetTitle("");
  dmFr->SetXTitle("m_{A'} (MeV/c^{2})");
  dmFr->SetYTitle("#deltam");
  dmFr->SetTitleSize(0.04,"X");
  dmFr->SetTitleSize(0.04,"Y");
  dmFr->SetLabelSize(0.04,"X");
  dmFr->SetLabelSize(0.035,"Y");
  dmFr->SetTitleOffset(0.75,"Y");

  TLegend* legDM = new TLegend(0.15,0.52,0.45,0.82);
  legDM->SetFillColor(0);
  legDM->SetTextSize(0.02);
  legDM->SetBorderSize(0);

  for( UInt_t ii = 0; ii < graphs_dm.size(); ++ii ){
    TString HRSStr;
    HRSStr.Form("%i",ii);
    if( ii == 0 ) legDM->AddEntry(graphs_dm[ii],EbStr+" pass; HRS"+HRSStr+" (4.5#circ)","p");
    else if( ii == graphs_dm.size()-1 ) legDM->AddEntry(graphs_dm[ii],EbStr+" pass; HRS"+HRSStr+" (5.5#circ)","p");
    else legDM->AddEntry(graphs_dm[ii],EbStr+" pass; HRS"+HRSStr,"p");
    graphs_dm[ii]->Draw("PL");
    graphs_dm[ii]->SaveAs("MassResolutionGraphs/dm_HRS"+HRSStr+"_"+EbStr+"GeV_T"+TStr+"_target"+Target+outApp+".root");
  }

  legDM->Draw();

  apexText->Draw("same");

  dmCanv->SaveAs("Plots_MassResolution/dm_proposal_style_TenSettings_"+EbStr+"GeV_T"+TStr+"_target"+Target+outApp+".pdf");

  /////////////////
  // Plot of Nexp
  /////////////////

  vec_BH[0]->SetLineColor(kCyan);
  vec_Rad[0]->SetLineColor(kCyan+1);
  vec_BH[1]->SetLineColor(kYellow);
  vec_Rad[1]->SetLineColor(kYellow+1);
  vec_BH[2]->SetLineColor(kGreen);
  vec_Rad[2]->SetLineColor(kGreen+3);
  vec_BH[3]->SetLineColor(kBlue);
  vec_Rad[3]->SetLineColor(kBlue+1);
  vec_BH[4]->SetLineColor(kMagenta);
  vec_Rad[4]->SetLineColor(kMagenta+1);
  vec_BH[5]->SetLineColor(kRed);
  vec_Rad[5]->SetLineColor(kRed+1);
  vec_BH[6]->SetLineColor(kOrange);
  vec_Rad[6]->SetLineColor(kOrange+3);
  vec_BH[7]->SetLineColor(kCyan+2);
  vec_Rad[7]->SetLineColor(kCyan+3);
  vec_BH[8]->SetLineColor(kPink);
  vec_Rad[8]->SetLineColor(kPink+1);
  vec_BH[9]->SetLineColor(kPink+2);
  vec_Rad[9]->SetLineColor(kPink+3);

  vec_BH[0]->SetFillColor(kCyan);
  vec_Rad[0]->SetFillColor(kCyan+1);
  vec_BH[1]->SetFillColor(kYellow);
  vec_Rad[1]->SetFillColor(kYellow+1);
  vec_BH[2]->SetFillColor(kGreen);
  vec_Rad[2]->SetFillColor(kGreen+3);
  vec_BH[3]->SetFillColor(kBlue);
  vec_Rad[3]->SetFillColor(kBlue+1);
  vec_BH[4]->SetFillColor(kMagenta);
  vec_Rad[4]->SetFillColor(kMagenta+1);
  vec_BH[5]->SetFillColor(kRed);
  vec_Rad[5]->SetFillColor(kRed+1);
  vec_BH[6]->SetFillColor(kOrange);
  vec_Rad[6]->SetFillColor(kOrange+3);
  vec_BH[7]->SetFillColor(kCyan+2);
  vec_Rad[7]->SetFillColor(kCyan+3);
  vec_BH[8]->SetFillColor(kPink);
  vec_Rad[8]->SetFillColor(kPink+1);
  vec_BH[9]->SetFillColor(kPink+2);
  vec_Rad[9]->SetFillColor(kPink+3);

  TString nExpXTitle = "m_{A'} (MeV/c^{2})";
  TString nExpYTitle = "Expected event yield";

  THStack* stack = new THStack("Stack","");
  for( UInt_t ii = 0; ii < vec_BH.size(); ++ii ){
    vec_BH[ii]->SetStats(0);
    vec_Rad[ii]->SetStats(0);
    vec_BH[ii]->GetXaxis()->SetTitle(nExpXTitle);
    vec_BH[ii]->GetYaxis()->SetTitle(nExpYTitle);
    vec_BH[ii]->GetXaxis()->SetTitleSize(0.04);
    vec_BH[ii]->SetTitle("");
    vec_Rad[ii]->GetXaxis()->SetTitle(nExpXTitle);
    vec_Rad[ii]->GetYaxis()->SetTitle(nExpYTitle);
    vec_Rad[ii]->GetXaxis()->SetTitleSize(0.04);
    vec_Rad[ii]->SetTitle("");
    stack->Add(vec_BH[ii]);
    stack->Add(vec_Rad[ii]);
  }

  TCanvas* cExp = new TCanvas("cExp","cExp",600,600);
  stack->Draw("hist");
  stack->GetXaxis()->SetTitle(nExpXTitle);
  stack->GetXaxis()->SetTitleSize(0.04);
  stack->SetTitle("");
  stack->GetYaxis()->SetTitle(nExpYTitle);
  stack->GetYaxis()->SetLabelSize(0.02);
  stack->GetYaxis()->SetTitleOffset(1.25);

  if( EbStr == "1" ){
    stack->GetXaxis()->SetRangeUser(50,140);
  }
  if( EbStr == "2" ){
    stack->GetXaxis()->SetRangeUser(100,250);
  }
  if( EbStr == "3" ){
    stack->GetXaxis()->SetRangeUser(150,350);
  }
  if( EbStr == "4" ){
    stack->GetXaxis()->SetRangeUser(200,500);
  }

  TLegend* leg1 = new TLegend(0.59,0.48,0.88,0.885);
  leg1->SetFillColor(0);
  leg1->SetTextSize(0.02);
  leg1->SetBorderSize(0);

  for( UInt_t ii = 0; ii < vec_BH.size(); ++ii ){
    TString HRSStr;
    HRSStr.Form("%i",ii);

    if( ii == 0 ){
      leg1->AddEntry(vec_BH[ii],"BH, "+EbStr+" pass; HRS"+HRSStr+" (4.5#circ)","f");
      leg1->AddEntry(vec_Rad[ii],"Rad, "+EbStr+" pass; HRS"+HRSStr+" (4.5#circ)","f");
    }else if( ii == vec_BH.size()-1 ){
      leg1->AddEntry(vec_BH[ii],"BH, "+EbStr+" pass; HRS"+HRSStr+" (5.5#circ)","f");
      leg1->AddEntry(vec_Rad[ii],"Rad, "+EbStr+" pass; HRS"+HRSStr+" (5.5#circ)","f");
    }else{
      leg1->AddEntry(vec_BH[ii],"BH, "+EbStr+" pass; HRS"+HRSStr,"f");
      leg1->AddEntry(vec_Rad[ii],"Rad, "+EbStr+" pass; HRS"+HRSStr,"f");
    }

  }

  leg1->Draw();

  apexText->Draw("same");

  for( UInt_t ii = 0; ii < vec_BH.size(); ++ii ){
    TString HRSStr;
    HRSStr.Form("%i",ii);
    vec_BH[ii]->SaveAs("Hists_NExpected/nExp_HRS"+HRSStr+"_BH_"+EbStr+"GeV_T"+TStr+"_target"+Target+outApp+".root");
    vec_Rad[ii]->SaveAs("Hists_NExpected/nExp_HRS"+HRSStr+"_Rad_"+EbStr+"GeV_T"+TStr+"_target"+Target+outApp+".root");
  }

  cExp->SaveAs("Plots_NExpected/nExp_TenSettings_"+EbStr+"GeV_T"+TStr+"_target"+Target+outApp+".pdf");


  f_BH_HRS0->Close();
  f_Rad_HRS0->Close();
  f_BH_HRS1->Close();
  f_Rad_HRS1->Close();
  f_BH_HRS2->Close();
  f_Rad_HRS2->Close();
  f_BH_HRS3->Close();
  f_Rad_HRS3->Close();
  f_BH_HRS4->Close();
  f_Rad_HRS4->Close();
  f_BH_HRS5->Close();
  f_Rad_HRS5->Close();
  f_BH_HRS6->Close();
  f_Rad_HRS6->Close();
  f_BH_HRS7->Close();
  f_Rad_HRS7->Close();
  f_BH_HRS8->Close();
  f_Rad_HRS8->Close();
  f_BH_HRS9->Close();
  f_Rad_HRS9->Close();

  //std::cout << "Here 1" << std::endl;

  delete mass;
  delete preF;
  delete B_BH;
  delete B_Rad;
  delete Neff;
  delete alpha;
  delete mA;
  delete deltam;
  delete eps2;
  delete EbeamF;
  //delete Tf;
  //delete Teff;
  delete theta;
  delete dThetaMS;
  /*
  delete invMassBH;
  delete tree_BH_4p5;
  delete tree_Rad_4p5;
  delete tree_BH_5;
  delete tree_Rad_5;
  delete tree_BH_5p5;
  delete tree_Rad_5p5;
  delete data_BH_4p5;
  delete data_Rad_4p5;
  delete data_BH_5;
  delete data_Rad_5;
  delete data_BH_5p5;
  delete data_Rad_5p5;
  */
}


////////////////////////////////////////////////////////////////////////////////////////
Double_t GetLuminosity( Double_t time, Double_t current, Double_t T, TString Target ){

  Double32_t coulomb = 1.0/(1.6e-19);
  //Double32_t NAvo    = 6.02e23;
  Double32_t Xzero   = 0.0;
  Double32_t AtomicW = 0.0;
  //Double32_t pb = 1e-36;

  Double32_t NAvoXpb = 6.02e-13;

  if( Target == "W" ){

    Xzero = 6.76;
    AtomicW = 184.0;

  }

  Double32_t N_el = time*24.0*60.0*60.0*current*(1e-6)*coulomb;

  std::cout << "N_el = " << N_el << std::endl;

  Double_t lumi = (N_el*T*Xzero/AtomicW)*NAvoXpb;

  //std::cout << "Lumi = " << lumi << std::endl;

  return lumi;

}
