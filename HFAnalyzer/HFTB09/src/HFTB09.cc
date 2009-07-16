// -*- C++ -*-
//
// Package:    HFTB09
// Class:      HFTB09
// 
/**\class HFTB09 HFTB09.cc HFAnalyzer/HFTB09/src/HFTB09.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Efe Yazgan
//         Created:  Thu Jul 16 14:04:07 CEST 2009
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TBDataFormats/HcalTBObjects/interface/HcalTBTriggerData.h"
#include "TBDataFormats/HcalTBObjects/interface/HcalTBEventPosition.h"
#include "TBDataFormats/HcalTBObjects/interface/HcalTBBeamCounters.h"
#include "TBDataFormats/HcalTBObjects/interface/HcalTBTiming.h"

//TFile Service
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

//
#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include <TROOT.h>
#include <TSystem.h>
#include "TFile.h"
#include <TCanvas.h>
#include <cmath>

//
// class decleration
//

class HFTB09 : public edm::EDAnalyzer {
   public:
      explicit HFTB09(const edm::ParameterSet&);
      ~HFTB09();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  edm::Service<TFileService> fs;
  TH1F* h_S1;
  TH1F* h_S2;
  TH1F* h_S3;
  TH1F* h_S4;

  TH1F* h_BH1;
  TH1F* h_BH2;
  TH1F* h_BH3;
  TH1F* h_BH4;

  TH2F* h_WCA;
  TH2F* h_WCB;
  TH2F* h_WCC;
  TH2F* h_WCD;
  TH2F* h_WCE; 

  TH2F*  h_WCA_cut;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HFTB09::HFTB09(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


HFTB09::~HFTB09()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HFTB09::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
 
  edm::Handle<HcalTBTriggerData> trigger;
  iEvent.getByLabel("tbUnpacker", trigger);
  edm::Handle<HcalTBBeamCounters> tb_qadc;
  iEvent.getByLabel("tbUnpacker", tb_qadc);
  edm::Handle<HcalTBEventPosition> tb_epos;
  iEvent.getByLabel("tbUnpacker", tb_epos);
  edm::Handle<HcalTBTiming> tb_time;
  iEvent.getByLabel("tbUnpacker", tb_time);

  h_S1->Fill(tb_qadc->S1adc());
  h_S2->Fill(tb_qadc->S2adc());
  h_S3->Fill(tb_qadc->S3adc());
  h_S4->Fill(tb_qadc->S4adc());

  h_BH1->Fill(tb_qadc->BH1adc());
  h_BH2->Fill(tb_qadc->BH2adc());
  h_BH3->Fill(tb_qadc->BH3adc());
  h_BH4->Fill(tb_qadc->BH4adc());

  std::vector<double> Axvec,Ayvec,Bxvec,Byvec,Cxvec,Cyvec,Dxvec,Dyvec;     
  std::vector<double> Exvec,Eyvec,Fxvec,Fyvec,Gxvec,Gyvec,Hxvec,Hyvec;
  tb_epos->getChamberHits('A',Axvec,Ayvec);     
  tb_epos->getChamberHits('B',Bxvec,Byvec);     
  tb_epos->getChamberHits('C',Cxvec,Cyvec);     
  tb_epos->getChamberHits('D',Dxvec,Dyvec);     
  tb_epos->getChamberHits('E',Exvec,Eyvec);     
  tb_epos->getChamberHits('F',Fxvec,Fyvec);     
  tb_epos->getChamberHits('G',Gxvec,Gyvec);     
  tb_epos->getChamberHits('H',Hxvec,Hyvec);

  if(Axvec.size()==1&&Ayvec.size()==1){ 
    h_WCA->Fill(Axvec[0],Ayvec[0]);     
    if (tb_qadc->BH1adc() > 150.) h_WCA_cut->Fill(Axvec[0],Ayvec[0]);
  }
  if(Bxvec.size()==1&&Byvec.size()==1) h_WCB->Fill(Bxvec[0],Byvec[0]);  
  if(Cxvec.size()==1&&Cyvec.size()==1) h_WCC->Fill(Cxvec[0],Cyvec[0]);  
  if(Dxvec.size()==1&&Dyvec.size()==1) h_WCD->Fill(Dxvec[0],Dyvec[0]);  
  if(Exvec.size()==1&&Eyvec.size()==1) h_WCE->Fill(Exvec[0],Eyvec[0]);  

  


}


// ------------ method called once each job just before starting event loop  ------------
void 
HFTB09::beginJob(const edm::EventSetup&)
{
  TFileDirectory BeamCntrs = fs->mkdir("Beam Counters");
  h_S1 = BeamCntrs.make<TH1F>("h_S1","S1",200,-50.,700.);
  h_S2 = BeamCntrs.make<TH1F>("h_S2","S2",200,-50.,700.);
  h_S3 = BeamCntrs.make<TH1F>("h_S3","S3",200,-50.,700.);
  h_S4 = BeamCntrs.make<TH1F>("h_S4","S4",200,-50.,1200.);
  h_BH1 = BeamCntrs.make<TH1F>("h_BH1","BH1",200,-5.,1000.);
  h_BH2 = BeamCntrs.make<TH1F>("h_BH2","BH2",200,-5.,1000.);
  h_BH3 = BeamCntrs.make<TH1F>("h_BH3","BH3",200,-5.,1000.);
  h_BH4 = BeamCntrs.make<TH1F>("h_BH4","BH4",200,-5.,1000.);
  h_WCA = BeamCntrs.make<TH2F>("WCApos","Wire Chamber A hits position",500,-55,75,500,-55,75);h_WCA->SetMarkerStyle(6);  
  h_WCB = BeamCntrs.make<TH2F>("WCBpos","Wire Chamber B hits position",500,-55,75,500,-55,75);h_WCB->SetMarkerStyle(6);  
  h_WCC = BeamCntrs.make<TH2F>("WCCpos","Wire Chamber C hits position",500,-55,75,500,-55,75);h_WCC->SetMarkerStyle(6);
  h_WCD = BeamCntrs.make<TH2F>("WCDpos","Wire Chamber D hits position",500,-55,75,500,-55,75);h_WCD->SetMarkerStyle(6);  
  h_WCE = BeamCntrs.make<TH2F>("WCEpos","Wire Chamber E hits position",500,-55,75,500,-55,75);h_WCE->SetMarkerStyle(6);  

  h_WCA_cut = BeamCntrs.make<TH2F>("WCApos_cut","Wire Chamber A hits position (cut)",500,-55,75,500,-55,75);h_WCA_cut->SetMarkerStyle(6);  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HFTB09::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFTB09);
