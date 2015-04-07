// Dear emacs, this is -*- c++ -*-
// $Id: CycleCreators.py 344 2012-12-13 13:10:53Z krasznaa $
#ifndef APEX_Cycle_H
#define APEX_Cycle_H

// SFrame include(s):
#include "core/include/SCycleBase.h"
#include "RooGlobalFunc.h"
#include <vector>
#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TROOT.h"

using namespace std;
/*
#ifdef _CINT_
#pragma link C++ class std::vector<std::vector<int> >;
#else
template class std::vector<std::vector<int> >;
#endif
*/
/**
 *   @short Put short description of class here
 *
 *          Put a longer description over here...
 *
 *  @author Put your name here
 * @version $Revision: 344 $
 */
class APEX_Cycle : public SCycleBase {

public:
  /// Default constructor
  APEX_Cycle();
  /// Default destructor
  ~APEX_Cycle();
  
  /// Function called at the beginning of the cycle
  virtual void BeginCycle() throw( SError );
  /// Function called at the end of the cycle
  virtual void EndCycle() throw( SError );
  
  /// Function called at the beginning of a new input data
  virtual void BeginInputData( const SInputData& ) throw( SError );
  /// Function called after finishing to process an input data
  virtual void EndInputData  ( const SInputData& ) throw( SError );
  
  /// Function called after opening each new input file
  virtual void BeginInputFile( const SInputData& ) throw( SError );
  
  /// Function called for every event
  virtual void ExecuteEvent( const SInputData&, Double_t ) throw( SError );

  Double_t GetSolidAngleAcceptance( Double_t dP, Double_t zBeam );
  Double_t GetHorizLo( Double_t zBeam );
  Double_t GetHorizHi( Double_t zBeam );
  
private:
  //
  // Put all your private variables here
  //

  std::string m_treeName;
  std::string m_runConfig;
  std::string m_accMapDir;

  //Int_t max_size = 20;

  Float_t m_pmin;
  Float_t m_pmax;

  Double_t m_invMass;

  Int_t m_elen;
  Float_t m_eleeta[20];
  Float_t m_elephi[20];
  Float_t m_elept[20];
  Float_t m_elem[20];
  Float_t m_elech[20];
  Int_t m_nucln;
  Float_t m_nucleta[20];
  Float_t m_nuclphi[20];
  Float_t m_nuclpt[20];
  Float_t m_nuclm[20];
  Float_t m_nuclch[20];
  Int_t m_An;
  Float_t m_Ae[20];
  Float_t m_Apx[20];
  Float_t m_Apy[20];
  Float_t m_Apz[20];
  Float_t m_Am[20];

  TH1F* m_h_nEvents;

  TH1F* m_h_cut_flow;

  TH1F* m_h_eplus_eta;
  TH1F* m_h_eplus_phi;
  TH1F* m_h_eplus_dP;
  TH1F* m_h_eplus_P;
  TH1F* m_h_eplus_pt;
  TH1F* m_h_eplus_px;
  TH1F* m_h_eplus_py;
  TH1F* m_h_eplus_pz;
  TH1F* m_h_eminus_eta;
  TH1F* m_h_eminus_phi;
  TH1F* m_h_eminus_dP;
  TH1F* m_h_eminus_P;
  TH1F* m_h_eminus_pt;
  TH1F* m_h_eminus_px;
  TH1F* m_h_eminus_py;
  TH1F* m_h_eminus_pz;

  TH1F* m_h_eplus_eminus_invmass;
  TH1F* m_h_eplus_eminus_invmass_bin1MeV;
  TH1F* m_h_eplus_eminus_invmass_bin0pt5MeV;

  TH2F* m_h_eplus_x_y;
  TH2F* m_h_eminus_x_y;

  TH2F* m_h_P_eplus_vs_eminus;

  TH2F* m_h_zBeam_vs_dP;
  TH2F* m_h_zBeam_vs_HRSSetting;
  TH2F* m_h_eplusPassed_vs_eminusPassed;

  TH2F* m_h_eplus_xLo_vs_xHi;
  TH2F* m_h_eminus_xLo_vs_xHi;

  TH2F* m_h_eplus_yLo_vs_yHi;
  TH2F* m_h_eminus_yLo_vs_yHi;

  TH2F* m_h_Omega_eplus_vs_eminus;
  
  // Macro adding the functions for dictionary generation
  ClassDef( APEX_Cycle, 0 );
  
}; // class APEX_Cycle

#endif // APEX_Cycle_H

