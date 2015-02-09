// $Id: CycleCreators.py 344 2012-12-13 13:10:53Z krasznaa $

// Local include(s):
#include "../include/APEX_Cycle.h"
#include <cmath>
#include <iostream>
#include <cstdio>
#include <vector>
#include <string.h>
#include <cstdlib>
#include "TMath.h"
#include "TMinuit.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TArrayD.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooFormulaVar.h"
#include "RooPoisson.h"
#include "RooProdPdf.h"
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

using namespace std;
ClassImp( APEX_Cycle );

APEX_Cycle::APEX_Cycle()
  : SCycleBase() {
  
  SetLogName( GetName() );
  DeclareProperty( "TreeName", m_treeName );
  DeclareProperty( "RunConfig", m_runConfig );
  DeclareProperty( "AccMapDir", m_accMapDir );

}

APEX_Cycle::~APEX_Cycle() {
  
}

void APEX_Cycle::BeginCycle() throw( SError ) {
  
  return;
  
}

void APEX_Cycle::EndCycle() throw( SError ) {
  
  return;
  
}

void APEX_Cycle::BeginInputData( const SInputData& ) throw( SError ) {

  m_h_nEvents = Book( TH1F( "h_nEvents", "No. events processed", 2, 0, 2) );

  m_h_cut_flow = Book( TH1F( "h_cut_flow", "Cut flow", 4, 0, 4) );
  m_h_cut_flow->GetXaxis()->SetBinLabel(1,"All events");
  m_h_cut_flow->GetXaxis()->SetBinLabel(2,"Two electrons");
  m_h_cut_flow->GetXaxis()->SetBinLabel(3,"Opposite sign");
  m_h_cut_flow->GetXaxis()->SetBinLabel(4,"Acceptance");

  m_h_eplus_eta = Book( TH1F( "h_eplus_eta", "e+ #eta", 100, 0, 4) );
  m_h_eminus_eta = Book( TH1F( "h_eminus_eta", "e- #eta", 100, 0, 4) );
  m_h_eplus_dP = Book( TH1F( "h_eplus_dP", "e+ dP", 500, -0.05, 0.05) );
  m_h_eminus_dP = Book( TH1F( "h_eminus_dP", "e- dP", 500, -0.05, 0.05) );
  m_h_eplus_P = Book( TH1F( "h_eplus_P", "e+ E", 1000, 0.45, 2.45) );
  m_h_eminus_P = Book( TH1F( "h_eminus_P", "e- E", 1000, 0.45, 2.45) );
  m_h_eplus_pt = Book( TH1F( "h_eplus_pt", "e+ p_{T}", 1000, -5, 5) );
  m_h_eminus_pt = Book( TH1F( "h_eminus_pt", "e- p_{T}", 1000, -5, 5) );
  m_h_eplus_px = Book( TH1F( "h_eplus_px", "e+ p_{x}", 1000, -5, 5) );
  m_h_eminus_px = Book( TH1F( "h_eminus_px", "e- p_{x}", 1000, -5, 5) );
  m_h_eplus_py = Book( TH1F( "h_eplus_py", "e+ p_{y}", 1000, -5, 5) );
  m_h_eminus_py = Book( TH1F( "h_eminus_py", "e- p_{y}", 1000, -5, 5) );
  m_h_eplus_pz = Book( TH1F( "h_eplus_pz", "e+ p_{z}", 1000, -5, 5) );
  m_h_eminus_pz = Book( TH1F( "h_eminus_pz", "e- p_{z}", 1000, -5, 5) );
  m_h_eplus_phi = Book( TH1F( "h_eplus_phi", "e+ #phi", 1000, -4, 4) );
  m_h_eminus_phi = Book( TH1F( "h_eminus_phi", "e- #phi", 1000, -4, 4) );
  m_h_eplus_eminus_invmass = Book( TH1F( "h_eplus_eminus_invmass", "e+e- invariant mass", 200, 0, 500) );
  m_h_eplus_eminus_invmass_bin1MeV = Book( TH1F( "h_eplus_eminus_invmass_bin1MeV", "e+e- invariant mass", 500, 0, 500) );
  m_h_eplus_eminus_invmass_bin0pt5MeV = Book( TH1F( "h_eplus_eminus_invmass_bin0pt5MeV", "e+e- invariant mass", 1000, 0, 500) );

  m_h_eplus_x_y = Book( TH2F( "h_eplus_x_y", "e+ x, y", 100, 0.05, 0.15, 100, -0.1, 0.1) );
  m_h_eminus_x_y = Book( TH2F( "h_eminus_x_y", "e- x, y", 100, -0.15, -0.05, 100, -0.1, 0.1) );

  m_h_P_eplus_vs_eminus = Book( TH2F( "h_P_eplus_vs_eminus", "eplus_P vs. eminus_P", 1000, 0.45, 2.45, 1000, 0.45, 2.45 ) );

  m_h_zBeam_vs_dP = Book( TH2F( "h_zBeam_vs_dP", "Zbeam vs. dP", 500, -0.05, 0.05, 520, -260.0, 260.0 ) );
  m_h_zBeam_vs_HRSSetting = Book( TH2F( "h_zBeam_vs_HRSSetting", "Zbeam vs. HRSSetting", 10, 0, 10, 520, -260.0, 260.0 ) );

  m_h_eplusPassed_vs_eminusPassed = Book( TH2F( "h_eplusPassed_vs_eminusPassed", "eplusPassed vs. eminusPassed", 2, 0, 2, 2, 0, 2 ) );

  m_h_eplus_xLo_vs_xHi = Book( TH2F( "h_eplus_xLo_vs_xHi", "e+ xHi, xLo", 100, 0.05, 0.15, 100, 0.05, 0.15) );
  m_h_eminus_xLo_vs_xHi = Book( TH2F( "h_eminus_xLo_vs_xHi", "e- xHi, xLo", 100, 0.05, 0.15, 100, 0.05, 0.15) );

  m_h_eplus_yLo_vs_yHi = Book( TH2F( "h_eplus_yLo_vs_yHi", "e+ yHi, yLo", 100, -0.1, 0.1, 100, -0.1, 0.1) );
  m_h_eminus_yLo_vs_yHi = Book( TH2F( "h_eminus_yLo_vs_yHi", "e- yHi, yLo", 100, -0.1, 0.1, 100, -0.1, 0.1) );

  m_h_Omega_eplus_vs_eminus = Book( TH2F( "h_Omega_eplus_vs_eminus", "h_Omega_eplus_vs_eminus", 30, 2, 5, 30, 2, 5 ) );

  // Variables for output TTree

  DeclareVariable( m_invMass, "mass" );
  
  return;

}

void APEX_Cycle::EndInputData( const SInputData& ) throw( SError ) {

  return;

}

void APEX_Cycle::BeginInputFile( const SInputData& ) throw( SError ) {

  ConnectVariable( m_treeName.c_str(), "pmin", m_pmin );
  ConnectVariable( m_treeName.c_str(), "pmax", m_pmax );
  ConnectVariable( m_treeName.c_str(), "elen", m_elen );
  ConnectVariable( m_treeName.c_str(), "eleeta", m_eleeta ); 
  ConnectVariable( m_treeName.c_str(), "elephi", m_elephi ); 
  ConnectVariable( m_treeName.c_str(), "elept", m_elept );
  ConnectVariable( m_treeName.c_str(), "elem", m_elem );
  ConnectVariable( m_treeName.c_str(), "elech", m_elech );
  ConnectVariable( m_treeName.c_str(), "nucln", m_nucln );
  ConnectVariable( m_treeName.c_str(), "nucleta", m_nucleta );
  ConnectVariable( m_treeName.c_str(), "nuclphi", m_nuclphi );
  ConnectVariable( m_treeName.c_str(), "nuclpt", m_nuclpt ); 
  ConnectVariable( m_treeName.c_str(), "nuclm", m_nuclm );
  ConnectVariable( m_treeName.c_str(), "nuclch", m_nuclch ); 
  ConnectVariable( m_treeName.c_str(), "An", m_An );
  ConnectVariable( m_treeName.c_str(), "Ae", m_Ae );
  ConnectVariable( m_treeName.c_str(), "Apx", m_Apx );
  ConnectVariable( m_treeName.c_str(), "Apy", m_Apy );
  ConnectVariable( m_treeName.c_str(), "Apz", m_Apz );
  ConnectVariable( m_treeName.c_str(), "Am", m_Am );

  return;

}

Double_t APEX_Cycle::GetSolidAngleAcceptance( Double_t dP, Double_t zBeam ){

  TString accDir(m_accMapDir);

  TFile* f1 = TFile::Open(accDir+"OmegaSolidAngle.root","READ");
  TGraph2D* tGr = (TGraph2D*) f1->Get("Graph2D");

  Double_t Omega = tGr->Interpolate(dP,zBeam);

  f1->Close();

  return Omega;

}

Double_t APEX_Cycle::GetHorizLo( Double_t zBeam ){

  TString accDir(m_accMapDir);

  TFile* f1 = TFile::Open(accDir+"ThetaHorizLo.root","READ");
  TGraph* tGr = (TGraph*) f1->Get("Graph");

  Double_t xLo = tGr->Eval(zBeam);

  f1->Close();

  return xLo;

}

Double_t APEX_Cycle::GetHorizHi( Double_t zBeam ){

  TString accDir(m_accMapDir);

  TFile* f1 = TFile::Open(accDir+"ThetaHorizHi.root","READ");
  TGraph* tGr = (TGraph*) f1->Get("Graph");

  Double_t xHi = tGr->Eval(zBeam);

  f1->Close();

  return xHi;

}


void APEX_Cycle::ExecuteEvent( const SInputData& sid, Double_t weight ) throw( SError ) {

  //std::cout << "Executing event" << std::endl;

  m_invMass = -999.0;

  m_h_nEvents->Fill( 1.5 );

  weight = 1.0;

  m_h_cut_flow->Fill("All events", weight);
  /*
  std::cout << "2" << std::endl;
  std::cout << "# electrons: " << m_elen << std::endl;
  std::cout << "  el0 pt container  ptr = " << m_elept << std::endl;
  //std::cout << "  el0 pt container size = " << m_elept.size() << std::endl;
  //std::cout << "  el0 pt container size = " << (*m_elept).size() << std::endl;
  //std::cout << "  el0 pt = " << (*m_elept)[0] << std::endl;
  std::cout << "  el0 pt = " << m_elept[0] << std::endl;
  //std::cout << "  el0 ch = " << (*m_elech)[0] << std::endl;
  */
  /*
  for( Int_t ii = 0; ii < m_elen; ++ii ){
    std::cout << "el # " << ii << "     ch = " << m_elech[ii] << std::endl;
  }
  */
  if( m_elen > 2){
    m_h_cut_flow->Fill("Two electrons", weight);
    
    if( m_elech[1]*m_elech[2] < 0 ){
      m_h_cut_flow->Fill("Opposite sign", weight);
    }
    
  }

  m_h_eplus_eta->Fill( m_eleeta[1], weight );
  m_h_eminus_eta->Fill( m_eleeta[2], weight );
  m_h_eplus_phi->Fill( m_elephi[1], weight );
  m_h_eminus_phi->Fill( m_elephi[2], weight );
  m_h_eplus_pt->Fill( m_elept[1], weight );
  m_h_eminus_pt->Fill( m_elept[2], weight );

  TLorentzVector eplus;
  TLorentzVector eminus;
  eplus.SetPtEtaPhiM( m_elept[1], m_eleeta[1], m_elephi[1], m_elem[1] );
  eminus.SetPtEtaPhiM( m_elept[2], m_eleeta[2], m_elephi[2], m_elem[2] );

  m_h_eplus_px->Fill( eplus.Px(), weight );
  m_h_eminus_px->Fill( eminus.Px(), weight );

  Double_t eplus_thx = atan( eplus.Px()/eplus.Pz() );
  Double_t eplus_thy = atan( eplus.Py()/eplus.Pz() );

  Double_t eminus_thx = atan( eminus.Px()/eminus.Pz() );
  Double_t eminus_thy = atan( eminus.Py()/eminus.Pz() );

  //std::cout << "e+ (x,y) = (" << eplus_thx << "," << eplus_thy << ")" << std::endl;
  //std::cout << "e- (x,y) = (" << eminus_thx << "," << eminus_thy << ")" << std::endl;

  m_h_eplus_x_y->Fill( eplus_thx, eplus_thy, weight );
  m_h_eminus_x_y->Fill( eminus_thx, eminus_thy, weight );

  //std::cout << "e+ E = " << eplus.E() << std::endl;
  //std::cout << "e- E = " << eminus.E() << std::endl;

  m_h_eplus_P->Fill( eplus.P(), weight );
  m_h_eminus_P->Fill( eminus.P(), weight );

  m_h_P_eplus_vs_eminus->Fill( eplus.P(), eminus.P(), weight );

  Double_t pcent = (m_pmax+m_pmin)/2.0;
  Double_t eplus_dP  = (eplus.P()-pcent)/pcent;
  Double_t eminus_dP = (eminus.P()-pcent)/pcent;

  m_h_eplus_dP->Fill( eplus_dP, weight );
  m_h_eminus_dP->Fill( eminus_dP, weight );

  // First approximation: Uniform acceptance across the horizontal range 3.5 to 7.1 degrees,
  // a.k.a, 0.06106 to 0.12386 radians

  ///////////////////////////////////
  // John Lerose acceptance numbers
  ///////////////////////////////////

  Double_t zBeamLo = -250.0; // In mm
  Double_t zBeamHi = 250.0; // In mm

  Int_t thisHRSSetting = -999;

  if( (sid.GetVersion()).Contains("nom5degHRS0") ) thisHRSSetting = 0;
  else if( (sid.GetVersion()).Contains("nom5degHRS1") ) thisHRSSetting = 1;
  else if( (sid.GetVersion()).Contains("nom5degHRS2") ) thisHRSSetting = 2;
  else if( (sid.GetVersion()).Contains("nom5degHRS3") ) thisHRSSetting = 3;
  else if( (sid.GetVersion()).Contains("nom5degHRS4") ) thisHRSSetting = 4;
  else if( (sid.GetVersion()).Contains("nom5degHRS5") ) thisHRSSetting = 5;
  else if( (sid.GetVersion()).Contains("nom5degHRS6") ) thisHRSSetting = 6;
  else if( (sid.GetVersion()).Contains("nom5degHRS7") ) thisHRSSetting = 7;
  else if( (sid.GetVersion()).Contains("nom5degHRS8") ) thisHRSSetting = 8;
  else if( (sid.GetVersion()).Contains("nom5degHRS9") ) thisHRSSetting = 9;
  else{
    std::cout << "Not a valid input file for the ten-setting setup." << std::endl;
    return;
  }

  // 1) Get z-beam, from the input filename (e.g., for zBeam = 0.015)

  Double_t zBeam = zBeamLo + (thisHRSSetting)*(zBeamHi-zBeamLo)/9.0;

  m_h_zBeam_vs_dP->Fill(eplus_dP,zBeam,1.0);
  m_h_zBeam_vs_dP->Fill(eminus_dP,zBeam,1.0);

  m_h_zBeam_vs_HRSSetting->Fill(thisHRSSetting,zBeam,1.0);

  // 2) Get the central angle, from the input filename (e.g., for 5degHRS)

  //Double_t angleCent = 0.087245*180.0/TMath::Pi();

  // 3) Using the P and dP of each particle, get the particle's solid angle acceptance,
  //    Omega, from the lookup TGraph

  Double_t eplus_Omega  = GetSolidAngleAcceptance(eplus_dP,zBeam)*1e-3;
  Double_t eminus_Omega = GetSolidAngleAcceptance(eminus_dP,zBeam)*1e-3;

  m_h_Omega_eplus_vs_eminus->Fill( eminus_Omega, eplus_Omega, 1.0 );

  // 4) Using zBeam, get horizontal acceptance

  Double_t xLo = GetHorizLo(zBeam);
  Double_t xHi = GetHorizHi(zBeam);

  Double_t xLoRad = xLo*TMath::Pi()/180.0;
  Double_t xHiRad = xHi*TMath::Pi()/180.0;

  // 5) Using each particle's Omega and the horizontal acceptance, get the +- vertical acceptance,
  //    and translate into yLo and yHi

  Double_t eplus_VertPM  = eplus_Omega/(xHiRad-xLoRad);
  Double_t eminus_VertPM = eminus_Omega/(xHiRad-xLoRad);

  //m_h_eplus_vs_eminus->Fill( xHi, xLo, 1.0 );

  Double_t eplus_yLo  = -1.0*eplus_VertPM;
  Double_t eminus_yLo = -1.0*eminus_VertPM;

  Double_t eplus_yHi  = eplus_VertPM;
  Double_t eminus_yHi = eminus_VertPM;

  m_h_eplus_xLo_vs_xHi->Fill( xHi, xLo, 1.0 );

  m_h_eplus_yLo_vs_yHi->Fill( eplus_yHi, eplus_yLo, 1.0 );
  m_h_eminus_yLo_vs_yHi->Fill( eminus_yHi, eminus_yLo, 1.0 );

  // 6) See if both eplus and eminus are accepted

  Bool_t eplusYes  = false;
  Bool_t eminusYes = false;

  if( eplus_thx > xLoRad && eplus_thx < xHiRad &&
      eplus_thy > eplus_yLo && eplus_thy < eplus_yHi ){
    eplusYes = true;
  }

  if( eminus_thx > -1.0*xHiRad && eminus_thx < -1.0*xLoRad &&
      eminus_thy > eminus_yLo && eminus_thy < eminus_yHi ){
    eminusYes = true;
  }

  m_h_eplusPassed_vs_eminusPassed->Fill( eminusYes, eplusYes, 1.0 );

  if( !eplusYes || !eminusYes ) return;

  /*
  if( pcent > eplus.P() ){
    std::cout << pcent << "   " << eplus.P() << std::endl;
  }
  if( pcent > eminus.P() ){
    std::cout << pcent << "   " << eminus.P() << std::endl;
  }
  */

  TLorentzVector ee;
  ee = eplus+eminus;

  Double_t eeM = ee.M();

  m_invMass = eeM*1e3;

  //std::cout << "eeM = " << eeM << std::endl;

  m_h_eplus_eminus_invmass->Fill( eeM*1e3, weight );
  m_h_eplus_eminus_invmass_bin1MeV->Fill( eeM*1e3, weight );
  m_h_eplus_eminus_invmass_bin0pt5MeV->Fill( eeM*1e3, weight );

  return;

}

