#include "TFile.h"
#include "TMath.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include "TProfile.h"
#include "math.h"
#include <fstream>
#include <string>
#include <iostream>
#include <TStyle.h>
using namespace std;
int main(int argc,char **argv){
  gROOT->SetStyle("Plain");

  TChain muon("prom/Muon/muonHB");
  //  /castor/cern.ch/user/e/efe/CRAFT/TrackerPoint_JC
  // muon.Add("/data1/efe/prompt_out.root");//TrkPointing224
  //  muon.Add("/data1/efe/Commissioning08_CRAFT_ALL_V9_225_ReReco_FromTrackerPointing_v1.root");
  //muon.Add("/data1/efe/CosmicMCBon_10GeV.root");
  //muon.Add("/data1/efe/CosmicMCBoff_10GeV.root");
  //muon.Add("/data1/efe/Boff_70195_70674_70675.root");

  //muon.Add("prompt_out.root");
  //  muon.Add("rfio:/castor/cern.ch/user/j/jdamgov/CRAFT/DATA/Bon_Commissioning08_CRAFT_ALL_V12_229_Tosca090322_ReReco_FromTrackerPointing_v1.root");
 

  //muon.Add("/data1/efe/Bon_Commissioning08_CRAFT_ALL_V12_229_Tosca090322_ReReco_FromTrackerPointing_v1.root");//DATA BON
  //muon.Add("/data1/efe/BON_10GeV_AllCMS_COSMMC_22X_V6_v1.root");//MC BON
  muon.Add("prompt_out.root");

  //muon.Add(".root");

  int MagnetEffectStudy = 0;
  int Bon = 1;
  /*
  ifstream pdg;
  float pdg_mom[1000];
  float rad_loss[1000];
  pdg_E[1000];
  pdg.open("pdg_dedx.txt");
  int pp = 0;
  while(1){
    if (!pdg.good()) break;
    pdg >> a1 >> a2 >> a3 >> a4 >> a5 >> a6 >> a7 >> a8 >> a9 >> a10 >> a11;
    pdg_mom[pp] = a2;
    pdg_E[pp]   = a8;
    rad_loss[pp] = a7;
  }
  */
  //-------------------corrections---------------------------------
  ifstream in;
  in.open("CRAFT08respMu_v2.26.txt");
  int phi_channel[80];
  float corr_HBp[80]={};
  float corr_HBm[80]={};
  float av_corr = 1.86876;
  int kkk = 0;
  while(1){
    if (!in.good() || kkk == 72) break;
    ++kkk;
    in >> phi_channel[kkk] >> corr_HBp[kkk] >> corr_HBm[kkk];
    //   cout<<kkk<<"   "<<phi_channel[kkk]<<"   "<<corr_HBp[kkk]<<endl;
    
  }
  for (int i=0;i<kkk;i++){
          if (corr_HBp[i] == 0 || i == 13){ 
	    corr_HBp[i] = av_corr;
	    corr_HBp[i] = 1;
	  }
          if (corr_HBp[i] != 0 && i != 13){ 
	    corr_HBp[i] = av_corr/corr_HBp[i];
	    corr_HBp[i] = 1;
	  }
          if (corr_HBm[i] == 0 || i == 64){
	    corr_HBm[i] = av_corr;
	    corr_HBm[i] = 1;
	  }
          if (corr_HBm[i] != 0 && i != 64){
	    corr_HBm[i] = av_corr/corr_HBm[i];
	    corr_HBm[i] = 1;
	  }
  }
  float lEHBtopP1_corr = -99999.;
  float lEHBbotP1_corr = -99999.;
  float Esum_top_minus_corr = -99999.;
  float Esum_top_plus_corr = -99999.;
  float Esum_bot_minus_corr = -99999.;
  float Esum_bot_plus_corr = -99999.;
  
  //---------------------------------------------------------------



  /////////////////////////////////////////////////////////
  // Parameters
  // dist from the phi tower edge in tow widths.
  double dTphi=0.1;
  /////////////////////////////////////////////////////////
  // Ntpl vars
  /////////////////////////////////////////////////////////
  Int_t hbheC,run, event, TriggerBit[4];
  // TK
  Float_t HBpedP[14][72];
  Float_t HBpedM[14][72];
  Float_t tk_mom;
  Float_t tk_momt;
  Int_t tk_ndof;
  Float_t tk_chi2;
  Int_t tk_lost;
  Int_t tk_charge;
  Float_t tk_d0;

  // DT
  Int_t NumDTtracks;
  Int_t DT_valid[20];
  Float_t DT_mom[20];
  Float_t DT_momt[20];
  Int_t DT_ndof[20];
  Int_t DT_lost[20];
  Int_t DT_charge[20];
  Float_t DT_chi2[20];
  Float_t DT_d0[20];
  Float_t DT_XinPos[20];
  Float_t DT_YinPos[20];
  Float_t DT_ZinPos[20];
  Float_t DT_XoutPos[20];
  Float_t DT_YoutPos[20];
  Float_t DT_ZoutPos[20];
  Float_t DT_XoutPosHB[20];
  Float_t DT_YoutPosHB[20];
  Float_t DT_ZoutPosHB[20];
  Float_t DT_etaHB[20];
  Float_t DT_phiHB[20];
  Int_t DT_ietaHB[20];
  Int_t DT_iphiHB[20];
  Int_t DT_isTop[20];

  // HB top
  Int_t NTowHBtop;
  Int_t isValidHBtop;
  Float_t ETowHBtop[30][3];
  Float_t ETowHBtopCr[30][3];
  Float_t TTowHBtop[30][3];
  Float_t EHBtop;
  Float_t EHBtopCr;
  Float_t THBtop;
  Int_t iPhiInTowHBtop;
  Int_t iPhiOutTowHBtop;
  Int_t iEtaInTowHBtop;
  Int_t iEtaOutTowHBtop;
  Float_t PhiInTowHBtop;
  Float_t PhiOutTowHBtop;
  Float_t EtaInTowHBtop;
  Float_t EtaOutTowHBtop;
  Float_t TK_XinPosHBtop;
  Float_t TK_YinPosHBtop;
  Float_t TK_ZinPosHBtop;
  Float_t TK_XoutPosHBtop;
  Float_t TK_YoutPosHBtop;
  Float_t TK_ZoutPosHBtop;
  Float_t TK_LengthHBtop;
  //HB Bottom
  Int_t NTowHBbot;
  Int_t isValidHBbot;
  Float_t ETowHBbot[30][3];
  Float_t ETowHBbotCr[30][3];
  Float_t TTowHBbot[30][3];
  Float_t EHBbot;
  Float_t EHBbotCr;
  Float_t THBbot;
  Int_t iPhiInTowHBbot;
  Int_t iPhiOutTowHBbot;
  Int_t iEtaInTowHBbot;
  Int_t iEtaOutTowHBbot;
  Float_t PhiInTowHBbot;
  Float_t PhiOutTowHBbot;
  Float_t EtaInTowHBbot;
  Float_t EtaOutTowHBbot;
  Float_t TK_XinPosHBbot;
  Float_t TK_YinPosHBbot;
  Float_t TK_ZinPosHBbot;
  Float_t TK_XoutPosHBbot;
  Float_t TK_YoutPosHBbot;
  Float_t TK_ZoutPosHBbot;
  Float_t TK_LengthHBbot;

  // Set branch addresses.
  muon.SetBranchAddress("run",&run);
  muon.SetBranchAddress("event",&event);

  muon.SetBranchAddress("run",  &run);
  muon.SetBranchAddress("event",  &event);
  muon.SetBranchAddress("hbheC",  &hbheC);
  muon.SetBranchAddress("TriggerBit",  TriggerBit);
  muon.SetBranchAddress("HBpedP",  HBpedP);
  muon.SetBranchAddress("HBpedM",  HBpedM);
  muon.SetBranchAddress("tk_mom", &tk_mom);
  muon.SetBranchAddress("tk_momt", &tk_momt);
  muon.SetBranchAddress("tk_ndof", &tk_ndof);
  muon.SetBranchAddress("tk_chi2", &tk_chi2);
  muon.SetBranchAddress("tk_lost", &tk_lost);
  muon.SetBranchAddress("tk_charge", &tk_charge);
  muon.SetBranchAddress("tk_d0", &tk_d0);
  muon.SetBranchAddress("NumDTtracks", &NumDTtracks);
  muon.SetBranchAddress("DT_valid",  DT_valid);
  muon.SetBranchAddress("DT_mom",  DT_mom);
  muon.SetBranchAddress("DT_momt",  DT_momt);
  muon.SetBranchAddress("DT_ndof",  DT_ndof);
  muon.SetBranchAddress("DT_lost",  DT_lost);
  muon.SetBranchAddress("DT_charge",  DT_charge);
  muon.SetBranchAddress("DT_chi2",  DT_chi2);
  muon.SetBranchAddress("DT_d0",  DT_d0);
  muon.SetBranchAddress("DT_XinPos", DT_XinPos);
  muon.SetBranchAddress("DT_YinPos", DT_YinPos);
  muon.SetBranchAddress("DT_ZinPos", DT_ZinPos);
  muon.SetBranchAddress("DT_XoutPos", DT_XoutPos);
  muon.SetBranchAddress("DT_YoutPos", DT_YoutPos);
  muon.SetBranchAddress("DT_ZoutPos", DT_ZoutPos);
  muon.SetBranchAddress("DT_XoutPosHB", DT_XoutPosHB);
  muon.SetBranchAddress("DT_YoutPosHB", DT_YoutPosHB);
  muon.SetBranchAddress("DT_ZoutPosHB", DT_ZoutPosHB);
  muon.SetBranchAddress("DT_etaHB", DT_etaHB);
  muon.SetBranchAddress("DT_phiHB", DT_phiHB);
  muon.SetBranchAddress("DT_ietaHB", DT_ietaHB);
  muon.SetBranchAddress("DT_iphiHB", DT_iphiHB);
  muon.SetBranchAddress("DT_isTop", DT_isTop);
  muon.SetBranchAddress("NTowHBtop",&NTowHBtop);
  muon.SetBranchAddress("isValidHBtop",&isValidHBtop);
  muon.SetBranchAddress("ETowHBtop",  ETowHBtop);
  muon.SetBranchAddress("ETowHBtopCr",  ETowHBtopCr);
  muon.SetBranchAddress("TTowHBtop",  TTowHBtop);
  muon.SetBranchAddress("EHBtop",  &EHBtop);
  muon.SetBranchAddress("EHBtopCr",  &EHBtopCr);
  muon.SetBranchAddress("THBtop",  &THBtop);
  muon.SetBranchAddress("iPhiInTowHBtop",  &iPhiInTowHBtop);
  muon.SetBranchAddress("iPhiOutTowHBtop",  &iPhiOutTowHBtop);
  muon.SetBranchAddress("iEtaInTowHBtop",  &iEtaInTowHBtop);
  muon.SetBranchAddress("iEtaOutTowHBtop",  &iEtaOutTowHBtop);
  muon.SetBranchAddress("PhiInTowHBtop",  &PhiInTowHBtop);
  muon.SetBranchAddress("PhiOutTowHBtop",  &PhiOutTowHBtop);
  muon.SetBranchAddress("EtaInTowHBtop",  &EtaInTowHBtop);
  muon.SetBranchAddress("EtaOutTowHBtop",  &EtaOutTowHBtop);
  muon.SetBranchAddress("TK_XinPosHBtop", &TK_XinPosHBtop);
  muon.SetBranchAddress("TK_YinPosHBtop", &TK_YinPosHBtop);
  muon.SetBranchAddress("TK_ZinPosHBtop", &TK_ZinPosHBtop);
  muon.SetBranchAddress("TK_XoutPosHBtop", &TK_XoutPosHBtop);
  muon.SetBranchAddress("TK_YoutPosHBtop", &TK_YoutPosHBtop);
  muon.SetBranchAddress("TK_ZoutPosHBtop", &TK_ZoutPosHBtop);
  muon.SetBranchAddress("TK_LengthHBtop", &TK_LengthHBtop);
  muon.SetBranchAddress("NTowHBbot",&NTowHBbot);
  muon.SetBranchAddress("isValidHBbot",&isValidHBbot);
  muon.SetBranchAddress("ETowHBbot",  ETowHBbot);
  muon.SetBranchAddress("ETowHBbotCr",  ETowHBbotCr);
  muon.SetBranchAddress("TTowHBbot",  TTowHBbot);
  muon.SetBranchAddress("EHBbot",  &EHBbot);
  muon.SetBranchAddress("EHBbotCr",  &EHBbotCr);
  muon.SetBranchAddress("THBbot",  &THBbot);
  muon.SetBranchAddress("iPhiInTowHBbot",  &iPhiInTowHBbot);
  muon.SetBranchAddress("iPhiOutTowHBbot",  &iPhiOutTowHBbot);
  muon.SetBranchAddress("iEtaInTowHBbot",  &iEtaInTowHBbot);
  muon.SetBranchAddress("iEtaOutTowHBbot",  &iEtaOutTowHBbot);
  muon.SetBranchAddress("PhiInTowHBbot",  &PhiInTowHBbot);
  muon.SetBranchAddress("PhiOutTowHBbot",  &PhiOutTowHBbot);
  muon.SetBranchAddress("EtaInTowHBbot",  &EtaInTowHBbot);
  muon.SetBranchAddress("EtaOutTowHBbot",  &EtaOutTowHBbot);
  muon.SetBranchAddress("TK_XinPosHBbot", &TK_XinPosHBbot);
  muon.SetBranchAddress("TK_YinPosHBbot", &TK_YinPosHBbot);
  muon.SetBranchAddress("TK_ZinPosHBbot", &TK_ZinPosHBbot);
  muon.SetBranchAddress("TK_XoutPosHBbot", &TK_XoutPosHBbot);
  muon.SetBranchAddress("TK_YoutPosHBbot", &TK_YoutPosHBbot);
  muon.SetBranchAddress("TK_ZoutPosHBbot", &TK_ZoutPosHBbot);
  muon.SetBranchAddress("TK_LengthHBbot", &TK_LengthHBbot);

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Out File
  TH1::AddDirectory(true);  
  //TFile *theFile    =new TFile("hist.root", "RECREATE");
  //TFile *theFile    =new TFile("boff_for_mag_eff_on_resp.root", "RECREATE");
  //TFile *theFile     =new TFile("boff_v9.root", "RECREATE");
  //TFile *theFile = new TFile("MC_Bon_out.root","RECREATE");

  //TFile *theFile = new TFile("outtest_MC.root","RECREATE");
  TFile *theFile = new TFile("outtest.root","RECREATE");
  theFile->cd();

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Out hist in out root file
  TProfile *pEminusPhi = new TProfile ("pEminusPhi","pEminusPhi",72,0.5,72.5);
  TProfile *pEplusPhi =  new TProfile ("pEplusPhi","pEplusPhi",72,0.5,72.5);

  TProfile *pEminusPhi_corr = new TProfile ("pEminusPhi_corr","HB- corrected",72,0.5,72.5);
  TProfile *pEplusPhi_corr =  new TProfile ("pEplusPhi_corr","HB+ corrected",72,0.5,72.5);

  TH2F *hIdPhiPlusVsE = new TH2F("hIdPhiPlusVsE","IdPhiTower vs Emuon, HB+",  60, -2, 10, 74, -0.5, 73.5); // x, y coordiants
  TH2F *hIdPhiMinusVsE = new TH2F("hIdPhiMinusVsE","IdPhiTower vs Emuon, HB-",  60, -2, 10, 74, -0.5, 73.5); // x, y coordiants
  TH2F *hIdPhiPlusVsE_corrected = new TH2F("hIdPhiPlusVsE_corrected","IdPhiTower vs Emuon, HB+",  60, -2, 10, 74, -0.5, 73.5); // x, y coordiants
  TH2F *hIdPhiMinusVsE_corrected = new TH2F("hIdPhiMinusVsE_corrected","IdPhiTower vs Emuon, HB-",  60, -2, 10, 74, -0.5, 73.5); // x, y coordiants

  TH1F *iphi_plus_lt_10GeV = new TH1F("iphi_plus_lt_10GeV","",72,-0.5,71.5);
  TH1F *iphi_minus_lt_10GeV = new TH1F("iphi_minus_lt_10GeV","",72,-0.5,71.5);
  TH1F *iphi_plus_gt_10GeV = new TH1F("iphi_plus_gt_10GeV","",72,-0.5,71.5);
  TH1F *iphi_minus_gt_10GeV = new TH1F("iphi_minus_gt_10GeV","",72,-0.5,71.5);

  pEminusPhi->SetLineColor(2);
  pEplusPhi->SetLineColor(4);

  Float_t muX[100];
  int nmuX=0;
  for(float ix=log(3);ix<log(1000);ix+=0.2){
      muX[nmuX]=exp(ix);
      nmuX++;
  }

  TProfile *hdEdXmuon = new TProfile("hdEdXmuon","P muon vs dE/dX",nmuX-1,muX);
  TProfile *hdEdXmuon_corr = new TProfile("hdEdXmuon_corr","P muon vs dE/dX (corrected)",nmuX-1,muX);

  TProfile *hdEdXmuon_top_minus = new TProfile("hdEdXmuon_top_minus","P muon vs dE/dX (top-)",nmuX-1,muX);
  TProfile *hdEdXmuon_top_plus = new TProfile("hdEdXmuon_top_plus","P muon vs dE/dX (top+)",nmuX-1,muX);
  TProfile *hdEdXmuon_bot_minus = new TProfile("hdEdXmuon_bot_minus","P muon vs dE/dX (bot-)",nmuX-1,muX);
  TProfile *hdEdXmuon_bot_plus = new TProfile("hdEdXmuon_bot_plus","P muon vs dE/dX (bot+)",nmuX-1,muX);

  TProfile *hdEdXmuon_top_minus_corr = new TProfile("hdEdXmuon_top_minus_corr","P muon vs dE/dX (top-) corrected",nmuX-1,muX);
  TProfile *hdEdXmuon_top_plus_corr = new TProfile("hdEdXmuon_top_plus_corr","P muon vs dE/dX (top+) corrected",nmuX-1,muX);
  TProfile *hdEdXmuon_bot_minus_corr = new TProfile("hdEdXmuon_bot_minus_corr","P muon vs dE/dX (bot-) corrected",nmuX-1,muX);
  TProfile *hdEdXmuon_bot_plus_corr = new TProfile("hdEdXmuon_bot_plus_corr","P muon vs dE/dX (bot+) corrected",nmuX-1,muX);

  TProfile *dedx_over_p_vs_p_top_minus_corr = new TProfile("dedx_over_p_vs_p_top_minus_corr","P muon vs dE/dX/P (top-) corrected",nmuX-1,muX);
  TProfile *dedx_over_p_vs_p_top_plus_corr  = new TProfile("dedx_over_p_vs_p_top_plus_corr","P muon vs dE/dX/P (top+) corrected",nmuX-1,muX);
  TProfile *dedx_over_p_vs_p_bot_minus_corr = new TProfile("dedx_over_p_vs_p_bot_minus_corr","P muon vs dE/dX/P (bot-) corrected",nmuX-1,muX);
  TProfile *dedx_over_p_vs_p_bot_plus_corr  = new TProfile("dedx_over_p_vs_p_bot_plus_corr","P muon vs dE/dX/P (bot+) corrected",nmuX-1,muX);

  TProfile *hdEdXmuon_top_minus_w = new TProfile("hdEdXmuon_top_minus_w","P muon vs dE/dX (top-)",nmuX-1,muX);
  TProfile *hdEdXmuon_top_plus_w = new TProfile("hdEdXmuon_top_plus_w","P muon vs dE/dX (top+)",nmuX-1,muX);
  TProfile *hdEdXmuon_bot_minus_w = new TProfile("hdEdXmuon_bot_minus_w","P muon vs dE/dX (bot-)",nmuX-1,muX);
  TProfile *hdEdXmuon_bot_plus_w = new TProfile("hdEdXmuon_bot_plus_w","P muon vs dE/dX (bot+)",nmuX-1,muX);

  TH1F *hp_log_top_minus = new TH1F("hp_log_top_minus","Top-",nmuX-1,muX);
  TH1F *hp_log_top_plus = new TH1F("hp_log_top_plus","Top+",nmuX-1,muX);
  TH1F *hp_log_bot_minus = new TH1F("hp_log_bot_minus","Bottom-",nmuX-1,muX);
  TH1F *hp_log_bot_plus = new TH1F("hp_log_bot_plus","Bottom+",nmuX-1,muX);
  TH1F *hp_log = new TH1F("hp_log","P",nmuX-1,muX);
  
  TH1F *hdEdXmuon_Integ_top_minus = new TH1F("hdEdXmuon_Integ_top_minus","HB Top-",100,-5.,20.);
  TH1F *hdEdXmuon_Integ_top_plus = new TH1F("hdEdXmuon_Integ_top_plus","HB Top+",100,-5.,20.);
  TH1F *hdEdXmuon_Integ_bot_minus = new TH1F("hdEdXmuon_Integ_bot_minus","HB Bottom-",100,-5.,20.);
  TH1F *hdEdXmuon_Integ_bot_plus = new TH1F("hdEdXmuon_Integ_bot_plus","HB Bottom+",100,-5.,20.);

  TH1F *hdEdXmuon_Integ_top_minus_corr = new TH1F("hdEdXmuon_Integ_top_minus_corr","HB Top- (corr)",100,-5.,20.);
  TH1F *hdEdXmuon_Integ_top_plus_corr = new TH1F("hdEdXmuon_Integ_top_plus_corr","HB Top+ (corr)",100,-5.,20.);
  TH1F *hdEdXmuon_Integ_bot_minus_corr = new TH1F("hdEdXmuon_Integ_bot_minus_corr","HB Bottom- (corr)",100,-5.,20.);
  TH1F *hdEdXmuon_Integ_bot_plus_corr = new TH1F("hdEdXmuon_Integ_bot_plus_corr","HB Bottom+ (corr)",100,-5.,20.);

  TH1F *hdEdXmuon_Integ = new TH1F("hdEdXmuon_Integ","",100,-5.,20.);

  TH2F *hphi_p_E_top_minus = new TH2F("hphi_p_E_top_minus","HB Top- ", 100, -0., 100.,72,0.5,72.5);
  TH2F *hphi_p_E_top_plus = new TH2F("hphi_p_E_top_plus","HB Top+ ", 100, -0., 100.,72,0.5,72.5);
  TH2F *hphi_p_E_bot_minus = new TH2F("hphi_p_E_bot_minus","HB Bottom- ", 100, -0., 100.,72,0.5,72.5);
  TH2F *hphi_p_E_bot_plus = new TH2F("hphi_p_E_bot_plus","HB Bottom+ ", 100, -0., 100.,72,0.5,72.5);

  TH1F *hp_top_minus = new TH1F("hp_top_minus","Top-",1000,0.,1000.);
  TH1F *hp_top_plus = new TH1F("hp_top_plus","Top+",1000,0.,1000.);
  TH1F *hp_bot_minus = new TH1F("hp_bot_minus","Bottom-",1000,0.,1000.);
  TH1F *hp_bot_plus = new TH1F("hp_bot_plus","Bottom+",1000,0.,1000.);

  TH2F *correlation_top_minus_vs_bot_minus = new TH2F("correlation_top_minus_vs_bot_minus","correlation_top_minus_vs_bot_minus",200,-20,20,200,-20,20);

  ////////////////////////////////////////////////////////////////////////////////////////////

  Int_t nevent = muon.GetEntries();

  double towW = TMath::Pi()/36.;

  int topminus = 0;
  int topplus = 0;
  int botminus = 0;
  int botplus = 0;


  ////////////////////////////////////////////////////////////////////////////////////////////
  // Evt loop
  int kk = 0;
  for (Int_t i=kk;i<nevent;i++) {

    double muon_energy[4]={0,0,0,0};
    //for tkpoint224
    /*
    if (i == 1463725) continue;
    if (i == 4395626) continue;
    if (i == 4558215) continue;
    */

    //for tk point 225 v9 
    
    /* 
    if (i == 37051) continue;
    if (i == 954254) continue;
    if (i == 1664985) continue;
    */

    if(i%10000==0) cout<<"event:     "<<i<<"/"<<nevent<<"        run:    "<<run<<endl;
    muon.GetEntry(i);

    if(hbheC<4800) continue; // hmm, wedges off?
    if(run==67126)continue; // crappy ped in first 150k evt        

    //    if (run!=67147 && run!=67173 && run !=67219 && run !=67141) continue;//bon
    //if (run!=70195) continue;//boff
    //if (run!=67141) continue;

    //if(i%1000==0) cout<<"event:     "<<i<<"/"<<nevent<<"        run:    "<<run<<endl;

    int usetop=1;
    int usebot=1;

    int fill_top_minus = 0;
    int fill_top_plus = 0;
    int fill_bot_minus = 0;
    int fill_bot_plus = 0;

    //    if(tk_mom>100)continue;

    
    //    if(MagnetEffectStudy == 0 && Bon == 1 && tk_mom < 10) continue;//!!!!!!!!!!
    if(tk_ndof < 7) continue;
    

    if(isValidHBtop<1) usetop=0;
    if(isValidHBbot<1) usebot=0;

    if(NTowHBtop>7)usetop=0;
    if(NTowHBbot>7)usebot=0;

    if(iPhiInTowHBtop!=iPhiOutTowHBtop)usetop=0;
    if(iPhiInTowHBbot!=iPhiOutTowHBbot)usebot=0;
       
    int lEtaMinTop=iEtaInTowHBtop;
    if(iEtaInTowHBtop>iEtaOutTowHBtop) lEtaMinTop=iEtaOutTowHBtop;
    lEtaMinTop-=1;
    if(lEtaMinTop==0) usetop=0;
    if(abs(lEtaMinTop)>14)usetop=0;
    if(abs(lEtaMinTop+NTowHBtop)>14)usetop=0;
    if(lEtaMinTop*(lEtaMinTop+NTowHBtop)<=0)usetop=0;

    int lEtaMinBot=iEtaInTowHBbot;
    if(iEtaInTowHBbot>iEtaOutTowHBbot)lEtaMinBot=iEtaOutTowHBbot;
    lEtaMinBot-=1;
    if(lEtaMinBot==0) usebot=0;
    if(abs(lEtaMinBot)>14)usebot=0;
    if(abs(lEtaMinBot+NTowHBbot)>14)usebot=0;
    if(lEtaMinBot*(lEtaMinBot+NTowHBbot)<0)usebot=0;

    // Some TK-DT matching here ...
    double RtkdtTop=99999;
    double RtkdtBot=99999;
    double PdtTop=-1;
    double PdtBot=-1;
    for(int im=0;im<NumDTtracks;im++)
      {
	if(DT_valid[im]!=1||DT_ndof[im]<25)continue;
	double tRtkdtTop=999999;
	double tRtkdtBot=999999;
	tRtkdtTop=sqrt(
		       pow((DT_XoutPosHB[im]-TK_XoutPosHBtop),2)+
		       pow((DT_YoutPosHB[im]-TK_YoutPosHBtop),2)+
		       pow((DT_ZoutPosHB[im]-TK_ZoutPosHBtop),2));
	tRtkdtBot=sqrt(
		       pow((DT_XoutPosHB[im]-TK_XoutPosHBbot),2)+
		       pow((DT_YoutPosHB[im]-TK_YoutPosHBbot),2)+
		       pow((DT_ZoutPosHB[im]-TK_ZoutPosHBbot),2));
    
	if(tRtkdtTop<RtkdtTop){RtkdtTop=tRtkdtTop;PdtTop=DT_mom[im];}
	if(tRtkdtBot<RtkdtBot){RtkdtBot=tRtkdtBot;PdtBot=DT_mom[im];}
      }
    if(RtkdtTop>15)usetop=0;
    if(RtkdtBot>15)usebot=0;
    //  Matching done

    //  Fix this h-whedges

    if(iPhiInTowHBtop<3||iPhiInTowHBtop>34)usetop=0;
    if(iPhiInTowHBbot<39||iPhiInTowHBbot>70)usebot=0;

    //  Well placed in phi in tower 
    double dPhiInTowHBtop = fabs(PhiInTowHBtop/towW - double(int(PhiInTowHBtop/towW)));
    double dPhiOutTowHBtop = fabs(PhiOutTowHBtop/towW - double(int(PhiOutTowHBtop/towW)));
    double dPhiInTowHBbot = fabs(PhiInTowHBbot/towW - double(int(PhiInTowHBbot/towW)));
    double dPhiOutTowHBbot = fabs(PhiOutTowHBbot/towW - double(int(PhiOutTowHBbot/towW)));

    if(dPhiInTowHBtop>0.5) dPhiInTowHBtop=1.-dPhiInTowHBtop;
    if(dPhiOutTowHBtop>0.5) dPhiOutTowHBtop=1.-dPhiOutTowHBtop;
    if(dPhiInTowHBbot>0.5) dPhiInTowHBbot=1.-dPhiInTowHBbot;
    if(dPhiOutTowHBbot>0.5) dPhiOutTowHBbot=1.-dPhiOutTowHBbot;

    if(dPhiInTowHBtop<dTphi||dPhiOutTowHBtop<dTphi)usetop=0;
    if(dPhiInTowHBbot<dTphi||dPhiOutTowHBbot<dTphi)usebot=0;

    float Esum = 0;

    if (usetop && lEtaMinTop < 0 && usebot && lEtaMinBot <0) {
      fill_top_minus = 1;
      fill_bot_minus = 1;
    }
    if (usetop && lEtaMinTop < 0 && usebot && lEtaMinBot >0) {
      fill_top_minus = 1;
      fill_bot_plus = 1;
    }
    if (usetop && lEtaMinTop > 0 && usebot && lEtaMinBot <0) {
      fill_top_plus = 1;
      fill_bot_minus = 1;
    }
    if (usetop && lEtaMinTop > 0 && usebot && lEtaMinBot >0) {
      fill_top_plus = 1;
      fill_bot_plus = 1;
    }

    double Esum_top_minus = 0;
    double Esum_top_plus = 0;

    //  Top is clean - doing  analysis 
    if(usetop)
      {
	double lEHBtopP1 = 0;
	Esum_top_minus = 0;
	Esum_top_plus = 0;
	double Etopminus = 0;
	double Etopplus  = 0;

        if (lEtaMinTop < 0 && tk_mom < 10.) iphi_minus_lt_10GeV->Fill(iPhiInTowHBtop);
	if (lEtaMinTop > 0 && tk_mom < 10.) iphi_plus_lt_10GeV->Fill(iPhiInTowHBtop);
	if (lEtaMinTop < 0 && tk_mom > 10.) iphi_minus_gt_10GeV->Fill(iPhiInTowHBtop);
	if (lEtaMinTop > 0 && tk_mom > 10.) iphi_plus_gt_10GeV->Fill(iPhiInTowHBtop);

	if ((MagnetEffectStudy == 0 && Bon == 1 && tk_mom < 100.) || (MagnetEffectStudy == 1 && Bon == 1) || Bon == 0){
	  // 1phi
	  for(int i=lEtaMinTop+1;i<(lEtaMinTop+NTowHBtop-1);i++) lEHBtopP1+=ETowHBtop[i-lEtaMinTop][1];

	  //-----------apply corrections for 1-phi---------------------------------------
	  if (lEtaMinTop < 0) lEHBtopP1_corr = lEHBtopP1*corr_HBm[iPhiInTowHBtop];
	  if (lEtaMinTop > 0) lEHBtopP1_corr = lEHBtopP1*corr_HBp[iPhiInTowHBtop];
	  //-----------------------------------------------------------------------------

	  

	  if(lEtaMinTop < 0){ 
	    pEminusPhi->Fill(iPhiInTowHBtop,lEHBtopP1/TK_LengthHBtop);
	    pEminusPhi_corr->Fill(iPhiInTowHBtop,lEHBtopP1_corr/TK_LengthHBtop);
	    hIdPhiMinusVsE->Fill(lEHBtopP1/TK_LengthHBtop,iPhiInTowHBtop);
	    hIdPhiMinusVsE_corrected->Fill(lEHBtopP1_corr/TK_LengthHBtop,iPhiInTowHBtop);
	  }
	  if(lEtaMinTop > 0){
	    pEplusPhi->Fill(iPhiInTowHBtop,lEHBtopP1/TK_LengthHBtop);
	    pEplusPhi_corr->Fill(iPhiInTowHBtop,lEHBtopP1_corr/TK_LengthHBtop);
	    hIdPhiPlusVsE->Fill(lEHBtopP1/TK_LengthHBtop,iPhiInTowHBtop);
	    hIdPhiPlusVsE_corrected->Fill(lEHBtopP1_corr/TK_LengthHBtop,iPhiInTowHBtop);
	  }
	}

	for(int i=lEtaMinTop;i<(lEtaMinTop+NTowHBtop);i++){ 
	  if (lEtaMinTop < 0){ 
	    Esum_top_minus += ETowHBtop[i-lEtaMinTop][0] + ETowHBtop[i-lEtaMinTop][1] + ETowHBtop[i-lEtaMinTop][2];
	    Etopminus = EHBtop;
	  }
	  if (lEtaMinTop > 0){
	    Esum_top_plus += ETowHBtop[i-lEtaMinTop][0] + ETowHBtop[i-lEtaMinTop][1] + ETowHBtop[i-lEtaMinTop][2];
	    Etopplus = EHBtop;
	  }
	}

	//-------------apply corrections--------------------------------------
	Esum_top_minus_corr = Esum_top_minus*corr_HBm[iPhiInTowHBtop];
	Esum_top_plus_corr  = Esum_top_plus*corr_HBp[iPhiInTowHBtop];
	//--------------------------------------------------------------------

	
	if (lEtaMinTop < 0){
	  if (iPhiInTowHBtop!= 13 && iPhiInTowHBtop!=64){
	    muon_energy[0] = Esum_top_minus/TK_LengthHBtop;
	    hdEdXmuon_top_minus->Fill(tk_mom,Esum_top_minus/TK_LengthHBtop);
	    hdEdXmuon_top_minus_corr->Fill(tk_mom,Esum_top_minus_corr/TK_LengthHBtop);
	    hdEdXmuon->Fill(tk_mom,Esum_top_minus/TK_LengthHBtop);
	    hdEdXmuon_corr->Fill(tk_mom,Esum_top_minus_corr/TK_LengthHBtop);
	    dedx_over_p_vs_p_top_minus_corr->Fill(tk_mom,Esum_top_minus_corr/(TK_LengthHBtop*(tk_mom)));
	    hdEdXmuon_Integ_top_minus->Fill(Esum_top_minus/TK_LengthHBtop);
	    hdEdXmuon_Integ_top_minus_corr->Fill(Esum_top_minus_corr/TK_LengthHBtop);
	    if (tk_mom > 7. && tk_mom < 1000.) hdEdXmuon_Integ->Fill(Esum_top_minus/TK_LengthHBtop);
	    hphi_p_E_top_minus->Fill(tk_mom,iPhiInTowHBtop,Esum_top_minus/TK_LengthHBtop);
	    hp_top_minus->Fill(tk_mom);
	    hp_log_top_minus->Fill(tk_mom);
	    hp_log->Fill(tk_mom);
	    ++topminus;
	  }
	}
	if (lEtaMinTop > 0){
	  if (iPhiInTowHBtop!= 13){
	    muon_energy[1] = Esum_top_plus/TK_LengthHBtop;
	    hdEdXmuon_top_plus->Fill(tk_mom,Esum_top_plus/TK_LengthHBtop);
	    hdEdXmuon_top_plus_corr->Fill(tk_mom,Esum_top_plus_corr/TK_LengthHBtop);
	    hdEdXmuon->Fill(tk_mom,Esum_top_plus/TK_LengthHBtop);
	    hdEdXmuon_corr->Fill(tk_mom,Esum_top_plus_corr/TK_LengthHBtop);
	    dedx_over_p_vs_p_top_plus_corr->Fill(tk_mom,Esum_top_plus_corr/(TK_LengthHBtop*(tk_mom)));
	    hdEdXmuon_Integ_top_plus->Fill(Esum_top_plus/TK_LengthHBtop);
	    hdEdXmuon_Integ_top_plus_corr->Fill(Esum_top_plus_corr/TK_LengthHBtop);
	    if (tk_mom > 7. && tk_mom < 1000.) hdEdXmuon_Integ->Fill(Esum_top_plus/TK_LengthHBtop);
	    hphi_p_E_top_plus->Fill(tk_mom,iPhiInTowHBtop,Esum_top_plus/TK_LengthHBtop);
	    hp_top_plus->Fill(tk_mom);
	    hp_log_top_plus->Fill(tk_mom);
	    hp_log->Fill(tk_mom);
	    ++topplus;
	  }
	}
	if (lEtaMinTop < 0) hdEdXmuon_top_minus_w->Fill(tk_mom,Etopminus/TK_LengthHBtop);
	if (lEtaMinTop > 0) hdEdXmuon_top_plus_w->Fill(tk_mom,Etopplus/TK_LengthHBtop);
      }

    double Esum_bot_minus = 0;
    double Esum_bot_plus = 0;

    //  Bottom is clean - doing  analysis 
    if(usebot)
      {
	double lEHBbotP1 = 0;
	Esum_bot_minus = 0;
	Esum_bot_plus = 0;
	double Ebotminus = 0;
	double Ebotplus = 0;

        if (lEtaMinBot < 0 && tk_mom < 10. && iPhiInTowHBbot!=13 && iPhiInTowHBbot!=64) iphi_minus_lt_10GeV->Fill(iPhiInTowHBbot);
	if (lEtaMinBot > 0 && tk_mom < 10. && iPhiInTowHBbot!=13) iphi_plus_lt_10GeV->Fill(iPhiInTowHBbot);
	if (lEtaMinBot < 0 && tk_mom > 10. && iPhiInTowHBbot!=13 && iPhiInTowHBbot!=64) iphi_minus_gt_10GeV->Fill(iPhiInTowHBbot);
	if (lEtaMinBot > 0 && tk_mom > 10. && iPhiInTowHBbot!=13) iphi_plus_gt_10GeV->Fill(iPhiInTowHBbot);

	if ((MagnetEffectStudy == 0 && Bon == 1 && tk_mom < 100.) || (MagnetEffectStudy == 1 && Bon == 1) || Bon == 0){
	  // 1phi
	  for(int i=lEtaMinBot+1;i<(lEtaMinBot+NTowHBbot-1);i++) lEHBbotP1+=ETowHBbot[i-lEtaMinBot][1];

	  //-----------apply corrections for 1-phi---------------------------------------
	  if (lEtaMinBot < 0) lEHBbotP1_corr = lEHBbotP1*corr_HBm[iPhiInTowHBbot];
	  if (lEtaMinBot > 0) lEHBbotP1_corr = lEHBbotP1*corr_HBp[iPhiInTowHBbot];
	  //-----------------------------------------------------------------------------

	  if(lEtaMinBot<0){
	    pEminusPhi->Fill(iPhiInTowHBbot,lEHBbotP1/TK_LengthHBbot);
	    pEminusPhi_corr->Fill(iPhiInTowHBbot,lEHBbotP1_corr/TK_LengthHBbot);
	    hIdPhiMinusVsE->Fill(lEHBbotP1/TK_LengthHBbot,iPhiInTowHBbot);
	    hIdPhiMinusVsE_corrected->Fill(lEHBbotP1_corr/TK_LengthHBbot,iPhiInTowHBbot);
	  }
	  if(lEtaMinBot>0){
	    pEplusPhi->Fill(iPhiInTowHBbot,lEHBbotP1/TK_LengthHBbot);
	    pEplusPhi_corr->Fill(iPhiInTowHBbot,lEHBbotP1_corr/TK_LengthHBbot);
	    hIdPhiPlusVsE->Fill(lEHBbotP1/TK_LengthHBbot,iPhiInTowHBbot);
	    hIdPhiPlusVsE_corrected->Fill(lEHBbotP1_corr/TK_LengthHBbot,iPhiInTowHBbot);
	  }
	}

	for(int i=lEtaMinBot;i<(lEtaMinBot+NTowHBbot);i++){ 
	  if (lEtaMinBot < 0){ 
	    Esum_bot_minus += ETowHBbot[i-lEtaMinBot][0] + ETowHBbot[i-lEtaMinBot][1] + ETowHBbot[i-lEtaMinBot][2];
	    Ebotminus = EHBbot;
	  }
	  if (lEtaMinBot > 0){
	    Esum_bot_plus  += ETowHBbot[i-lEtaMinBot][0] + ETowHBbot[i-lEtaMinBot][1] + ETowHBbot[i-lEtaMinBot][2];
	    Ebotplus = EHBbot;
	  }
	}


	//-------------apply corrections--------------------------------------
	Esum_bot_minus_corr = Esum_bot_minus*corr_HBm[iPhiInTowHBbot];
	Esum_bot_plus_corr  = Esum_bot_plus*corr_HBp[iPhiInTowHBbot];
	//--------------------------------------------------------------------


	if (lEtaMinBot < 0){ 
	  if (iPhiInTowHBbot!=13 && iPhiInTowHBbot!=64){
	    muon_energy[2] = Esum_bot_minus/TK_LengthHBbot;
	    hdEdXmuon_bot_minus->Fill(tk_mom,Esum_bot_minus/TK_LengthHBbot);
	    hdEdXmuon_bot_minus_corr->Fill(tk_mom,Esum_bot_minus_corr/TK_LengthHBbot);
	    hdEdXmuon->Fill(tk_mom,Esum_bot_minus/TK_LengthHBbot);
	    hdEdXmuon_corr->Fill(tk_mom,Esum_bot_minus_corr/TK_LengthHBbot);
	    dedx_over_p_vs_p_bot_minus_corr->Fill(tk_mom,Esum_bot_minus_corr/(TK_LengthHBbot*(tk_mom)));
	    hdEdXmuon_Integ_bot_minus->Fill(Esum_bot_minus/TK_LengthHBbot);
	    hdEdXmuon_Integ_bot_minus_corr->Fill(Esum_bot_minus_corr/TK_LengthHBbot);
	    if (tk_mom > 7. && tk_mom < 1000.) hdEdXmuon_Integ->Fill(Esum_bot_minus/TK_LengthHBbot);
	    hphi_p_E_bot_minus->Fill(tk_mom,iPhiInTowHBbot,Esum_bot_minus/TK_LengthHBbot);
	    hp_bot_minus->Fill(tk_mom);
	    hp_log_bot_minus->Fill(tk_mom);
	    hp_log->Fill(tk_mom);
	    ++botminus;
	  }
	}
	if (lEtaMinBot > 0){
	  if (iPhiInTowHBbot!=13){
	    muon_energy[3] = Esum_bot_plus/TK_LengthHBbot;
	    hdEdXmuon_bot_plus->Fill(tk_mom,Esum_bot_plus/TK_LengthHBbot);
	    hdEdXmuon_bot_plus_corr->Fill(tk_mom,Esum_bot_plus_corr/TK_LengthHBbot);
	    hdEdXmuon->Fill(tk_mom,Esum_bot_plus/TK_LengthHBbot);
	    hdEdXmuon_corr->Fill(tk_mom,Esum_bot_plus_corr/TK_LengthHBbot);
	    dedx_over_p_vs_p_bot_plus_corr->Fill(tk_mom,Esum_bot_plus_corr/(TK_LengthHBbot*(tk_mom)));
	    hdEdXmuon_Integ_bot_plus->Fill(Esum_bot_plus/TK_LengthHBbot);
	    hdEdXmuon_Integ_bot_plus_corr->Fill(Esum_bot_plus_corr/TK_LengthHBbot);
	    if (tk_mom > 7. && tk_mom < 1000.) hdEdXmuon_Integ->Fill(Esum_bot_plus/TK_LengthHBbot);
	    hphi_p_E_bot_plus->Fill(tk_mom,iPhiInTowHBbot,Esum_bot_plus/TK_LengthHBbot);
	    hp_bot_plus->Fill(tk_mom);
	    hp_log_bot_plus->Fill(tk_mom);
	    hp_log->Fill(tk_mom);
	    ++botplus;
	  }
	}

	if (lEtaMinBot < 0){
	  hdEdXmuon_bot_minus_w->Fill(tk_mom,Ebotminus/TK_LengthHBtop);
	}
	if (lEtaMinBot > 0){
	  hdEdXmuon_bot_plus_w->Fill(tk_mom,Ebotplus/TK_LengthHBtop);
	}

       

      }



    if (fill_top_minus && fill_bot_minus) correlation_top_minus_vs_bot_minus->Fill(muon_energy[0],muon_energy[2]);

    /*
    if (fill_top_minus && fill_bot_minus){
      cout<<fill_top_minus<<fill_top_plus<<fill_bot_minus<<fill_bot_plus<<endl;
      cout<<muon_energy[0]<<"  "<<muon_energy[1]<<"   "<<muon_energy[2]<<"    "<<muon_energy[3]<<endl;
      cout<<Esum_top_minus/TK_LengthHBtop<<"  "<<Esum_top_plus/TK_LengthHBtop<<"    "<<Esum_bot_minus/TK_LengthHBbot<<"    "<<Esum_bot_plus/TK_LengthHBbot<<endl;
    }
    */

    // Analysis done
  } //end event loop

  // Dump out root
  theFile->Write();
  theFile->Close();
  cout<<"Top- ="<<topminus<<endl;
  cout<<"Top+ ="<<topplus<<endl;
  cout<<"Bot- ="<<botminus<<endl;
  cout<<"Bot+ ="<<botplus<<endl;
}
