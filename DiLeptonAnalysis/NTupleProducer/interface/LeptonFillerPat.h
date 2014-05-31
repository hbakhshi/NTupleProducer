#ifndef __DiLeptonAnalysis_NTupleProducer_LeptonFillerPat_H__
#define __DiLeptonAnalysis_NTupleProducer_LeptonFillerPat_H__
// 
// Package: NTupleProducer
// Class:   LeptonFillerPat
//
/* class LeptonFillerPat
   LeptonFillerPat.h
   Description:  generic class for basic jet dumper

*/
//
// $Id: LeptonFillerPat.h,v 1.5.2.8 2012/05/24 12:21:46 paktinat Exp $
//
//

#include <string>
#include <vector>

#include "TTree.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "DiLeptonAnalysis/NTupleProducer/interface/FillerBase.h"

template <class LeptonType>
class LeptonFillerPat : public FillerBase {
public:
  /// Constructor: set pointer to tree
  LeptonFillerPat<LeptonType>( const edm::ParameterSet&, const bool& isRealData );
  virtual ~LeptonFillerPat(void) {}

  /// Define all branches
  virtual const std::vector<filler::PPair> declareProducts(void);
  /// Reset all branch containers
  virtual void resetProducts(void);
  /// Fill all branches
  virtual void fillProducts(edm::Event&, const edm::EventSetup& );
  /// Put products in the event data
  virtual void putProducts( edm::Event& );

  enum Type {
    El, Mu, Tau, unknown
  };
  Type fType;

private:


  /// retrieve specific lepton information
  virtual void getSpecific(const LeptonType& lepton) {}

  /// Set and get jet type
  void setType( const Type& type ) { fType = type; }
  const Type getType(void) const { return fType; }
  
  //- Configuration parameters
  edm::InputTag fTag; 
  
  // Pre-selection
  double fMinpt;
  double fMaxeta;

  size_t gMaxnobjs;
  
  // Tree variables
  std::auto_ptr<int>     fTMaxLepExc;
  std::auto_ptr<int>     fTNObjsTot;
  std::auto_ptr<int>     fTNObjs;
  std::auto_ptr<std::vector<float> >  fTPx;
  std::auto_ptr<std::vector<float> >  fTPy;
  std::auto_ptr<std::vector<float> >  fTPz;
  std::auto_ptr<std::vector<float> >  fTPt;
  std::auto_ptr<std::vector<float> >  fTPtErr;
  std::auto_ptr<std::vector<float> >  fTEta;
  std::auto_ptr<std::vector<float> >  fTPhi;
  std::auto_ptr<std::vector<float> >  fTE;
  std::auto_ptr<std::vector<float> >  fTEt;
  std::auto_ptr<std::vector<int> >    fTCharge;
  std::auto_ptr<std::vector<int> >    fTIsPFTau;
  std::auto_ptr<std::vector<int> >    fTDecayMode;
  std::auto_ptr<std::vector<float> >  fTVz; 
  std::auto_ptr<std::vector<float> >  fTEmFraction; 
  std::auto_ptr<std::vector<float> >  fTJetPt;
  std::auto_ptr<std::vector<float> >  fTJetEta;
  std::auto_ptr<std::vector<float> >  fTJetPhi;
  std::auto_ptr<std::vector<float> >  fTJetMass;
  std::auto_ptr<std::vector<float> >  fTLeadingTkPt;
  std::auto_ptr<std::vector<float> >  fTLeadingNeuPt;
  std::auto_ptr<std::vector<float> >  fTLeadingTkHcalenergy;
  std::auto_ptr<std::vector<float> >  fTLeadingTkEcalenergy;
  std::auto_ptr<std::vector<int> >    fTNumChargedHadronsSignalCone;
  std::auto_ptr<std::vector<int> >    fTNumNeutralHadronsSignalCone;
  std::auto_ptr<std::vector<int> >    fTNumPhotonsSignalCone;
  std::auto_ptr<std::vector<int> >    fTNumParticlesSignalCone;
  std::auto_ptr<std::vector<int> >    fTNumChargedHadronsIsoCone;
  std::auto_ptr<std::vector<int> >    fTNumNeutralHadronsIsoCone;
  std::auto_ptr<std::vector<int> >    fTNumPhotonsIsolationCone;
  std::auto_ptr<std::vector<int> >    fTNumParticlesIsolationCone;
  std::auto_ptr<std::vector<float> >  fTPtSumChargedParticlesIsoCone;
  std::auto_ptr<std::vector<float> >  fTPtSumPhotonsIsoCone;

//defining the tauid variables, by Hamed
std::auto_ptr<std::vector<float> >  fTtauIDdecayModeFindingNewDMs;
std::auto_ptr<std::vector<float> >  fTtauIDdecayModeFindingOldDMs;
std::auto_ptr<std::vector<float> >  fTtauIDdecayModeFinding;
std::auto_ptr<std::vector<float> >  fTtauIDbyLooseIsolation;
std::auto_ptr<std::vector<float> >  fTtauIDbyVLooseCombinedIsolationDeltaBetaCorr;
std::auto_ptr<std::vector<float> >  fTtauIDbyLooseCombinedIsolationDeltaBetaCorr;
std::auto_ptr<std::vector<float> >  fTtauIDbyMediumCombinedIsolationDeltaBetaCorr;
std::auto_ptr<std::vector<float> >  fTtauIDbyTightCombinedIsolationDeltaBetaCorr;
std::auto_ptr<std::vector<float> >  fTtauIDbyCombinedIsolationDeltaBetaCorrRaw;
std::auto_ptr<std::vector<float> >  fTtauIDbyLooseCombinedIsolationDeltaBetaCorr3Hits;
std::auto_ptr<std::vector<float> >  fTtauIDbyMediumCombinedIsolationDeltaBetaCorr3Hits;
std::auto_ptr<std::vector<float> >  fTtauIDbyTightCombinedIsolationDeltaBetaCorr3Hits;
std::auto_ptr<std::vector<float> >  fTtauIDbyCombinedIsolationDeltaBetaCorrRaw3Hits;
std::auto_ptr<std::vector<float> >  fTtauIDchargedIsoPtSum;
std::auto_ptr<std::vector<float> >  fTtauIDneutralIsoPtSum;
std::auto_ptr<std::vector<float> >  fTtauIDpuCorrPtSum;
std::auto_ptr<std::vector<float> >  fTtauIDbyIsolationMVA3oldDMwoLTraw;
std::auto_ptr<std::vector<float> >  fTtauIDbyVLooseIsolationMVA3oldDMwoLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyLooseIsolationMVA3oldDMwoLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyMediumIsolationMVA3oldDMwoLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyTightIsolationMVA3oldDMwoLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyVTightIsolationMVA3oldDMwoLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyVVTightIsolationMVA3oldDMwoLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyIsolationMVA3oldDMwLTraw;
std::auto_ptr<std::vector<float> >  fTtauIDbyVLooseIsolationMVA3oldDMwLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyLooseIsolationMVA3oldDMwLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyMediumIsolationMVA3oldDMwLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyTightIsolationMVA3oldDMwLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyVTightIsolationMVA3oldDMwLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyVVTightIsolationMVA3oldDMwLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyIsolationMVA3newDMwoLTraw;
std::auto_ptr<std::vector<float> >  fTtauIDbyVLooseIsolationMVA3newDMwoLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyLooseIsolationMVA3newDMwoLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyMediumIsolationMVA3newDMwoLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyTightIsolationMVA3newDMwoLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyVTightIsolationMVA3newDMwoLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyVVTightIsolationMVA3newDMwoLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyIsolationMVA3newDMwLTraw;
std::auto_ptr<std::vector<float> >  fTtauIDbyVLooseIsolationMVA3newDMwLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyLooseIsolationMVA3newDMwLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyMediumIsolationMVA3newDMwLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyTightIsolationMVA3newDMwLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyVTightIsolationMVA3newDMwLT;
std::auto_ptr<std::vector<float> >  fTtauIDbyVVTightIsolationMVA3newDMwLT;
std::auto_ptr<std::vector<float> >  fTtauIDagainstElectronLoose;
std::auto_ptr<std::vector<float> >  fTtauIDagainstElectronMedium;
std::auto_ptr<std::vector<float> >  fTtauIDagainstElectronTight;
std::auto_ptr<std::vector<float> >  fTtauIDagainstElectronMVA5raw;
std::auto_ptr<std::vector<float> >  fTtauIDagainstElectronMVA5category;
std::auto_ptr<std::vector<float> >  fTtauIDagainstElectronVLooseMVA5;
std::auto_ptr<std::vector<float> >  fTtauIDagainstElectronLooseMVA5;
std::auto_ptr<std::vector<float> >  fTtauIDagainstElectronMediumMVA5;
std::auto_ptr<std::vector<float> >  fTtauIDagainstElectronTightMVA5;
std::auto_ptr<std::vector<float> >  fTtauIDagainstElectronVTightMVA5;
std::auto_ptr<std::vector<float> >  fTtauIDagainstElectronDeadECAL;
std::auto_ptr<std::vector<float> >  fTtauIDagainstMuonLoose;
std::auto_ptr<std::vector<float> >  fTtauIDagainstMuonMedium;
std::auto_ptr<std::vector<float> >  fTtauIDagainstMuonTight;
std::auto_ptr<std::vector<float> >  fTtauIDagainstMuonLoose2;
std::auto_ptr<std::vector<float> >  fTtauIDagainstMuonMedium2;
std::auto_ptr<std::vector<float> >  fTtauIDagainstMuonTight2;
std::auto_ptr<std::vector<float> >  fTtauIDagainstMuonLoose3;
std::auto_ptr<std::vector<float> >  fTtauIDagainstMuonTight3;
std::auto_ptr<std::vector<float> >  fTtauIDagainstMuonMVAraw;
std::auto_ptr<std::vector<float> >  fTtauIDagainstMuonLooseMVA;
std::auto_ptr<std::vector<float> >  fTtauIDagainstMuonMediumMVA;
std::auto_ptr<std::vector<float> >  fTtauIDagainstMuonTightMVA;
//END



  std::auto_ptr<std::vector<int> >    fTID80;
  std::auto_ptr<std::vector<int> >    fTID85;
  std::auto_ptr<std::vector<int> >    fTID90;
  std::auto_ptr<std::vector<int> >    fTID95;
  std::auto_ptr<std::vector<int> >    fTMuNMatches;

  std::auto_ptr<std::vector<float> > fTParticleIso;
  std::auto_ptr<std::vector<float> > fTChargedHadronIso;
  std::auto_ptr<std::vector<float> > fTNeutralHadronIso;
  std::auto_ptr<std::vector<float> > fTPhotonIso;

};

typedef LeptonFillerPat<pat::Muon>     PatMuonFiller;
typedef LeptonFillerPat<pat::Electron> PatElectronFiller;
typedef LeptonFillerPat<pat::Tau>      PatTauFiller;


//________________________________________________________________________________________
template <class LeptonType>
LeptonFillerPat<LeptonType>::LeptonFillerPat( const edm::ParameterSet& config, const bool& isRealData )
  : FillerBase( config, isRealData )
{

  // Retrieve configuration parameters
  std::string leptontype    = config.getParameter<std::string>("type");
  fMinpt                    = config.getParameter<double>("sel_minpt");
  fMaxeta                   = config.getParameter<double>("sel_maxeta");
  gMaxnobjs                 = config.getParameter<uint>("maxnobjs");
  
  fTag                      = config.getParameter<edm::InputTag>("tag");

  edm::LogVerbatim("NTP") << " ==> LeptonFillerPat Constructor - " << fPrefix;
  edm::LogVerbatim("NTP") << "  Input Tag:      " << fTag.label();
  edm::LogVerbatim("NTP") << "  Max n(objs):    " << gMaxnobjs;
  edm::LogVerbatim("NTP") << "  Min pt:         " << fMinpt;
  edm::LogVerbatim("NTP") << "  Max eta:        " << fMaxeta;
  edm::LogVerbatim("NTP") << "---------------------------------";

  if      (leptontype == "electron")   setType(El);
  else if (leptontype == "muon")       setType(Mu);
  else if (leptontype == "tau")        setType(Tau);
	  else {
    setType(unknown);
    edm::LogWarning("NTP") << "!! Don't know Lepton Type !!" << leptontype;
  }
  
}


//________________________________________________________________________________________
template <class LeptonType>
void LeptonFillerPat<LeptonType>::fillProducts(edm::Event& iEvent,
                                               const edm::EventSetup& iSetup ) {

  // Retrieve collection
  edm::Handle<edm::View<LeptonType> > collection;
  iEvent.getByLabel(fTag,collection);
    
  size_t pfqi(0);  // Index of qualified leptons
  for (typename edm::View<LeptonType>::const_iterator it = collection->begin(); 
       it != collection->end(); ++it ) {
    // Check if maximum number of leptons is exceeded already:
    if(pfqi >= gMaxnobjs){
      edm::LogWarning("NTP") << "@SUB=analyze()"
                             << "Maximum number of " << fPrefix << " leptons exceeded";
      *fTMaxLepExc = 1;
      break;
    }
    (*fTNObjsTot)++;

    // PfMuon preselection:
    if (it->pt() < fMinpt) continue;
    if (fabs(it->eta()) > fMaxeta) continue;
        
    const LeptonType& lepton = (*it);
    fTPx    ->push_back( lepton.px() );
    fTPy    ->push_back( lepton.py() );
    fTPz    ->push_back( lepton.pz() );
    fTPt    ->push_back( lepton.pt() );
    fTEta   ->push_back( lepton.eta() );
    fTPhi   ->push_back( lepton.phi() );
    fTE     ->push_back( lepton.energy() );
    fTEt    ->push_back( lepton.et() );
    fTCharge->push_back( lepton.charge() );
          
    fTParticleIso     ->push_back( (lepton.chargedHadronIso()+lepton.neutralHadronIso()+lepton.photonIso())/lepton.pt() );
    fTChargedHadronIso->push_back( lepton.chargedHadronIso() );
    fTNeutralHadronIso->push_back( lepton.neutralHadronIso() );
    fTPhotonIso       ->push_back( lepton.photonIso() );
  
    getSpecific(lepton);
          
    ++pfqi;

  }
  *fTNObjs = pfqi;

}

//________________________________________________________________________________________
template <class LeptonType>
const std::vector<filler::PPair> LeptonFillerPat<LeptonType>::declareProducts(void) {

  addProduct("MaxLepExc",typeid(*fTMaxLepExc));
  addProduct("NObjsTot", typeid(*fTNObjsTot));
  addProduct("NObjs",    typeid(*fTNObjs));

  addProduct("Px",       typeid(*fTPx));
  addProduct("Py",       typeid(*fTPy));
  addProduct("Pz",       typeid(*fTPz));
  addProduct("Pt",       typeid(*fTPt));
  addProduct("E",        typeid(*fTE));
  addProduct("Et",       typeid(*fTEt));
  addProduct("Eta",      typeid(*fTEta));
  addProduct("Phi",      typeid(*fTPhi));
  addProduct("Charge",   typeid(*fTCharge));

  addProduct("ParticleIso",      typeid(*fTParticleIso));
  addProduct("ChargedHadronIso", typeid(*fTChargedHadronIso));
  addProduct("NeutralHadronIso", typeid(*fTNeutralHadronIso));
  addProduct("PhotonIso",        typeid(*fTPhotonIso));
  
  if(fType == Tau){
    addProduct("IsPFTau", typeid(*fTIsPFTau));
    addProduct("DecayMode", typeid(*fTDecayMode));
    addProduct("Vz",        typeid(*fTVz)); 
    addProduct("EmFraction",typeid(*fTEmFraction)); 
    addProduct("JetPt",     typeid(*fTJetPt));
    addProduct("JetEta",    typeid(*fTJetEta));
    addProduct("JetPhi",    typeid(*fTJetPhi));
    addProduct("JetMass",   typeid(*fTJetMass));
    addProduct("LeadingTkPt", typeid(*fTLeadingTkPt));
    addProduct("LeadingNeuPt",typeid(*fTLeadingNeuPt));
    addProduct("LeadingTkHcalenergy", typeid(*fTLeadingTkHcalenergy));
    addProduct("LeadingTkEcalenergy", typeid(*fTLeadingTkEcalenergy));
    addProduct("NumChargedHadronsSignalCone", typeid(*fTNumChargedHadronsSignalCone));
    addProduct("NumNeutralHadronsSignalCone", typeid(*fTNumNeutralHadronsSignalCone));
    addProduct("NumPhotonsSignalCone",     typeid(*fTNumPhotonsSignalCone));
    addProduct("NumParticlesSignalCone",   typeid(*fTNumParticlesSignalCone));
    addProduct("NumChargedHadronsIsoCone", typeid(*fTNumChargedHadronsIsoCone));
    addProduct("NumNeutralHadronsIsoCone", typeid(*fTNumNeutralHadronsIsoCone));
    addProduct("NumPhotonsIsolationCone",  typeid(*fTNumPhotonsIsolationCone));
    addProduct("NumParticlesIsolationCone",typeid(*fTNumParticlesIsolationCone));
    addProduct("PtSumChargedParticlesIsoCone", typeid(*fTPtSumChargedParticlesIsoCone));
    addProduct("PtSumPhotonsIsoCone",      typeid(*fTPtSumPhotonsIsoCone));


//tauid addProducts, by Hamed
addProduct("decayModeFindingNewDMs", typeid(*fTtauIDdecayModeFindingNewDMs));
addProduct("decayModeFindingOldDMs", typeid(*fTtauIDdecayModeFindingOldDMs));
addProduct("decayModeFinding", typeid(*fTtauIDdecayModeFinding));
addProduct("byLooseIsolation", typeid(*fTtauIDbyLooseIsolation));
addProduct("byVLooseCombinedIsolationDeltaBetaCorr", typeid(*fTtauIDbyVLooseCombinedIsolationDeltaBetaCorr));
addProduct("byLooseCombinedIsolationDeltaBetaCorr", typeid(*fTtauIDbyLooseCombinedIsolationDeltaBetaCorr));
addProduct("byMediumCombinedIsolationDeltaBetaCorr", typeid(*fTtauIDbyMediumCombinedIsolationDeltaBetaCorr));
addProduct("byTightCombinedIsolationDeltaBetaCorr", typeid(*fTtauIDbyTightCombinedIsolationDeltaBetaCorr));
addProduct("byCombinedIsolationDeltaBetaCorrRaw", typeid(*fTtauIDbyCombinedIsolationDeltaBetaCorrRaw));
addProduct("byLooseCombinedIsolationDeltaBetaCorr3Hits", typeid(*fTtauIDbyLooseCombinedIsolationDeltaBetaCorr3Hits));
addProduct("byMediumCombinedIsolationDeltaBetaCorr3Hits", typeid(*fTtauIDbyMediumCombinedIsolationDeltaBetaCorr3Hits));
addProduct("byTightCombinedIsolationDeltaBetaCorr3Hits", typeid(*fTtauIDbyTightCombinedIsolationDeltaBetaCorr3Hits));
addProduct("byCombinedIsolationDeltaBetaCorrRaw3Hits", typeid(*fTtauIDbyCombinedIsolationDeltaBetaCorrRaw3Hits));
addProduct("chargedIsoPtSum", typeid(*fTtauIDchargedIsoPtSum));
addProduct("neutralIsoPtSum", typeid(*fTtauIDneutralIsoPtSum));
addProduct("puCorrPtSum", typeid(*fTtauIDpuCorrPtSum));
addProduct("byIsolationMVA3oldDMwoLTraw", typeid(*fTtauIDbyIsolationMVA3oldDMwoLTraw));
addProduct("byVLooseIsolationMVA3oldDMwoLT", typeid(*fTtauIDbyVLooseIsolationMVA3oldDMwoLT));
addProduct("byLooseIsolationMVA3oldDMwoLT", typeid(*fTtauIDbyLooseIsolationMVA3oldDMwoLT));
addProduct("byMediumIsolationMVA3oldDMwoLT", typeid(*fTtauIDbyMediumIsolationMVA3oldDMwoLT));
addProduct("byTightIsolationMVA3oldDMwoLT", typeid(*fTtauIDbyTightIsolationMVA3oldDMwoLT));
addProduct("byVTightIsolationMVA3oldDMwoLT", typeid(*fTtauIDbyVTightIsolationMVA3oldDMwoLT));
addProduct("byVVTightIsolationMVA3oldDMwoLT", typeid(*fTtauIDbyVVTightIsolationMVA3oldDMwoLT));
addProduct("byIsolationMVA3oldDMwLTraw", typeid(*fTtauIDbyIsolationMVA3oldDMwLTraw));
addProduct("byVLooseIsolationMVA3oldDMwLT", typeid(*fTtauIDbyVLooseIsolationMVA3oldDMwLT));
addProduct("byLooseIsolationMVA3oldDMwLT", typeid(*fTtauIDbyLooseIsolationMVA3oldDMwLT));
addProduct("byMediumIsolationMVA3oldDMwLT", typeid(*fTtauIDbyMediumIsolationMVA3oldDMwLT));
addProduct("byTightIsolationMVA3oldDMwLT", typeid(*fTtauIDbyTightIsolationMVA3oldDMwLT));
addProduct("byVTightIsolationMVA3oldDMwLT", typeid(*fTtauIDbyVTightIsolationMVA3oldDMwLT));
addProduct("byVVTightIsolationMVA3oldDMwLT", typeid(*fTtauIDbyVVTightIsolationMVA3oldDMwLT));
addProduct("byIsolationMVA3newDMwoLTraw", typeid(*fTtauIDbyIsolationMVA3newDMwoLTraw));
addProduct("byVLooseIsolationMVA3newDMwoLT", typeid(*fTtauIDbyVLooseIsolationMVA3newDMwoLT));
addProduct("byLooseIsolationMVA3newDMwoLT", typeid(*fTtauIDbyLooseIsolationMVA3newDMwoLT));
addProduct("byMediumIsolationMVA3newDMwoLT", typeid(*fTtauIDbyMediumIsolationMVA3newDMwoLT));
addProduct("byTightIsolationMVA3newDMwoLT", typeid(*fTtauIDbyTightIsolationMVA3newDMwoLT));
addProduct("byVTightIsolationMVA3newDMwoLT", typeid(*fTtauIDbyVTightIsolationMVA3newDMwoLT));
addProduct("byVVTightIsolationMVA3newDMwoLT", typeid(*fTtauIDbyVVTightIsolationMVA3newDMwoLT));
addProduct("byIsolationMVA3newDMwLTraw", typeid(*fTtauIDbyIsolationMVA3newDMwLTraw));
addProduct("byVLooseIsolationMVA3newDMwLT", typeid(*fTtauIDbyVLooseIsolationMVA3newDMwLT));
addProduct("byLooseIsolationMVA3newDMwLT", typeid(*fTtauIDbyLooseIsolationMVA3newDMwLT));
addProduct("byMediumIsolationMVA3newDMwLT", typeid(*fTtauIDbyMediumIsolationMVA3newDMwLT));
addProduct("byTightIsolationMVA3newDMwLT", typeid(*fTtauIDbyTightIsolationMVA3newDMwLT));
addProduct("byVTightIsolationMVA3newDMwLT", typeid(*fTtauIDbyVTightIsolationMVA3newDMwLT));
addProduct("byVVTightIsolationMVA3newDMwLT", typeid(*fTtauIDbyVVTightIsolationMVA3newDMwLT));
addProduct("againstElectronLoose", typeid(*fTtauIDagainstElectronLoose));
addProduct("againstElectronMedium", typeid(*fTtauIDagainstElectronMedium));
addProduct("againstElectronTight", typeid(*fTtauIDagainstElectronTight));
addProduct("againstElectronMVA5raw", typeid(*fTtauIDagainstElectronMVA5raw));
addProduct("againstElectronMVA5category", typeid(*fTtauIDagainstElectronMVA5category));
addProduct("againstElectronVLooseMVA5", typeid(*fTtauIDagainstElectronVLooseMVA5));
addProduct("againstElectronLooseMVA5", typeid(*fTtauIDagainstElectronLooseMVA5));
addProduct("againstElectronMediumMVA5", typeid(*fTtauIDagainstElectronMediumMVA5));
addProduct("againstElectronTightMVA5", typeid(*fTtauIDagainstElectronTightMVA5));
addProduct("againstElectronVTightMVA5", typeid(*fTtauIDagainstElectronVTightMVA5));
addProduct("againstElectronDeadECAL", typeid(*fTtauIDagainstElectronDeadECAL));
addProduct("againstMuonLoose", typeid(*fTtauIDagainstMuonLoose));
addProduct("againstMuonMedium", typeid(*fTtauIDagainstMuonMedium));
addProduct("againstMuonTight", typeid(*fTtauIDagainstMuonTight));
addProduct("againstMuonLoose2", typeid(*fTtauIDagainstMuonLoose2));
addProduct("againstMuonMedium2", typeid(*fTtauIDagainstMuonMedium2));
addProduct("againstMuonTight2", typeid(*fTtauIDagainstMuonTight2));
addProduct("againstMuonLoose3", typeid(*fTtauIDagainstMuonLoose3));
addProduct("againstMuonTight3", typeid(*fTtauIDagainstMuonTight3));
addProduct("againstMuonMVAraw", typeid(*fTtauIDagainstMuonMVAraw));
addProduct("againstMuonLooseMVA", typeid(*fTtauIDagainstMuonLooseMVA));
addProduct("againstMuonMediumMVA", typeid(*fTtauIDagainstMuonMediumMVA));
addProduct("againstMuonTightMVA", typeid(*fTtauIDagainstMuonTightMVA));
//END

  }else if(fType == El){
    addProduct("ID95", typeid(*fTID95));
    addProduct("ID90", typeid(*fTID90));
    addProduct("ID85", typeid(*fTID85));
    addProduct("ID80", typeid(*fTID80));
  }else if(fType == Mu){
    addProduct("PtErr"   , typeid(*fTPtErr));
    addProduct("NMatches", typeid(*fTMuNMatches));
  }

  return typeList;

}

//________________________________________________________________________________________
template <class LeptonType>
void LeptonFillerPat<LeptonType>::putProducts(edm::Event& e) {

  e.put(fTMaxLepExc,fullName("MaxLepExc"));
  e.put(fTNObjsTot,fullName("NObjsTot"));
  e.put(fTNObjs,fullName("NObjs"));

  e.put(fTPx,fullName("Px"));
  e.put(fTPy,fullName("Py"));
  e.put(fTPz,fullName("Pz"));
  e.put(fTPt,fullName("Pt"));
  e.put(fTE,fullName("E"));
  e.put(fTEt,fullName("Et"));
  e.put(fTEta,fullName("Eta"));
  e.put(fTPhi,fullName("Phi"));
  e.put(fTCharge,fullName("Charge"));

  e.put(fTParticleIso,fullName("ParticleIso"));
  e.put(fTChargedHadronIso,fullName("ChargedHadronIso"));
  e.put(fTNeutralHadronIso,fullName("NeutralHadronIso"));
  e.put(fTPhotonIso,fullName("PhotonIso"));
  
  if(fType == Tau){
    e.put(fTDecayMode,fullName("DecayMode"));
    e.put(fTIsPFTau,fullName("IsPFTau"));
    e.put(fTVz,fullName("Vz")); 
    e.put(fTEmFraction,fullName("EmFraction")); 
    e.put(fTJetPt,fullName("JetPt"));
    e.put(fTJetEta,fullName("JetEta"));
    e.put(fTJetPhi,fullName("JetPhi"));
    e.put(fTJetMass,fullName("JetMass"));
    e.put(fTLeadingTkPt,fullName("LeadingTkPt"));
    e.put(fTLeadingNeuPt,fullName("LeadingNeuPt"));
    e.put(fTLeadingTkHcalenergy,fullName("LeadingTkHcalenergy"));
    e.put(fTLeadingTkEcalenergy,fullName("LeadingTkEcalenergy"));
    e.put(fTNumChargedHadronsSignalCone,fullName("NumChargedHadronsSignalCone"));
    e.put(fTNumNeutralHadronsSignalCone,fullName("NumNeutralHadronsSignalCone"));
    e.put(fTNumPhotonsSignalCone,fullName("NumPhotonsSignalCone"));
    e.put(fTNumParticlesSignalCone,fullName("NumParticlesSignalCone"));
    e.put(fTNumChargedHadronsIsoCone,fullName("NumChargedHadronsIsoCone"));
    e.put(fTNumNeutralHadronsIsoCone,fullName("NumNeutralHadronsIsoCone"));
    e.put(fTNumPhotonsIsolationCone,fullName("NumPhotonsIsolationCone"));
    e.put(fTNumParticlesIsolationCone,fullName("NumParticlesIsolationCone"));
    e.put(fTPtSumChargedParticlesIsoCone,fullName("PtSumChargedParticlesIsoCone"));
    e.put(fTPtSumPhotonsIsoCone,fullName("PtSumPhotonsIsoCone"));

//tauid puts, by Hamed
e.put(fTtauIDdecayModeFindingNewDMs,fullName("decayModeFindingNewDMs"));
e.put(fTtauIDdecayModeFindingOldDMs,fullName("decayModeFindingOldDMs"));
e.put(fTtauIDdecayModeFinding,fullName("decayModeFinding"));
e.put(fTtauIDbyLooseIsolation,fullName("byLooseIsolation"));
e.put(fTtauIDbyVLooseCombinedIsolationDeltaBetaCorr,fullName("byVLooseCombinedIsolationDeltaBetaCorr"));
e.put(fTtauIDbyLooseCombinedIsolationDeltaBetaCorr,fullName("byLooseCombinedIsolationDeltaBetaCorr"));
e.put(fTtauIDbyMediumCombinedIsolationDeltaBetaCorr,fullName("byMediumCombinedIsolationDeltaBetaCorr"));
e.put(fTtauIDbyTightCombinedIsolationDeltaBetaCorr,fullName("byTightCombinedIsolationDeltaBetaCorr"));
e.put(fTtauIDbyCombinedIsolationDeltaBetaCorrRaw,fullName("byCombinedIsolationDeltaBetaCorrRaw"));
e.put(fTtauIDbyLooseCombinedIsolationDeltaBetaCorr3Hits,fullName("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
e.put(fTtauIDbyMediumCombinedIsolationDeltaBetaCorr3Hits,fullName("byMediumCombinedIsolationDeltaBetaCorr3Hits"));
e.put(fTtauIDbyTightCombinedIsolationDeltaBetaCorr3Hits,fullName("byTightCombinedIsolationDeltaBetaCorr3Hits"));
e.put(fTtauIDbyCombinedIsolationDeltaBetaCorrRaw3Hits,fullName("byCombinedIsolationDeltaBetaCorrRaw3Hits"));
e.put(fTtauIDchargedIsoPtSum,fullName("chargedIsoPtSum"));
e.put(fTtauIDneutralIsoPtSum,fullName("neutralIsoPtSum"));
e.put(fTtauIDpuCorrPtSum,fullName("puCorrPtSum"));
e.put(fTtauIDbyIsolationMVA3oldDMwoLTraw,fullName("byIsolationMVA3oldDMwoLTraw"));
e.put(fTtauIDbyVLooseIsolationMVA3oldDMwoLT,fullName("byVLooseIsolationMVA3oldDMwoLT"));
e.put(fTtauIDbyLooseIsolationMVA3oldDMwoLT,fullName("byLooseIsolationMVA3oldDMwoLT"));
e.put(fTtauIDbyMediumIsolationMVA3oldDMwoLT,fullName("byMediumIsolationMVA3oldDMwoLT"));
e.put(fTtauIDbyTightIsolationMVA3oldDMwoLT,fullName("byTightIsolationMVA3oldDMwoLT"));
e.put(fTtauIDbyVTightIsolationMVA3oldDMwoLT,fullName("byVTightIsolationMVA3oldDMwoLT"));
e.put(fTtauIDbyVVTightIsolationMVA3oldDMwoLT,fullName("byVVTightIsolationMVA3oldDMwoLT"));
e.put(fTtauIDbyIsolationMVA3oldDMwLTraw,fullName("byIsolationMVA3oldDMwLTraw"));
e.put(fTtauIDbyVLooseIsolationMVA3oldDMwLT,fullName("byVLooseIsolationMVA3oldDMwLT"));
e.put(fTtauIDbyLooseIsolationMVA3oldDMwLT,fullName("byLooseIsolationMVA3oldDMwLT"));
e.put(fTtauIDbyMediumIsolationMVA3oldDMwLT,fullName("byMediumIsolationMVA3oldDMwLT"));
e.put(fTtauIDbyTightIsolationMVA3oldDMwLT,fullName("byTightIsolationMVA3oldDMwLT"));
e.put(fTtauIDbyVTightIsolationMVA3oldDMwLT,fullName("byVTightIsolationMVA3oldDMwLT"));
e.put(fTtauIDbyVVTightIsolationMVA3oldDMwLT,fullName("byVVTightIsolationMVA3oldDMwLT"));
e.put(fTtauIDbyIsolationMVA3newDMwoLTraw,fullName("byIsolationMVA3newDMwoLTraw"));
e.put(fTtauIDbyVLooseIsolationMVA3newDMwoLT,fullName("byVLooseIsolationMVA3newDMwoLT"));
e.put(fTtauIDbyLooseIsolationMVA3newDMwoLT,fullName("byLooseIsolationMVA3newDMwoLT"));
e.put(fTtauIDbyMediumIsolationMVA3newDMwoLT,fullName("byMediumIsolationMVA3newDMwoLT"));
e.put(fTtauIDbyTightIsolationMVA3newDMwoLT,fullName("byTightIsolationMVA3newDMwoLT"));
e.put(fTtauIDbyVTightIsolationMVA3newDMwoLT,fullName("byVTightIsolationMVA3newDMwoLT"));
e.put(fTtauIDbyVVTightIsolationMVA3newDMwoLT,fullName("byVVTightIsolationMVA3newDMwoLT"));
e.put(fTtauIDbyIsolationMVA3newDMwLTraw,fullName("byIsolationMVA3newDMwLTraw"));
e.put(fTtauIDbyVLooseIsolationMVA3newDMwLT,fullName("byVLooseIsolationMVA3newDMwLT"));
e.put(fTtauIDbyLooseIsolationMVA3newDMwLT,fullName("byLooseIsolationMVA3newDMwLT"));
e.put(fTtauIDbyMediumIsolationMVA3newDMwLT,fullName("byMediumIsolationMVA3newDMwLT"));
e.put(fTtauIDbyTightIsolationMVA3newDMwLT,fullName("byTightIsolationMVA3newDMwLT"));
e.put(fTtauIDbyVTightIsolationMVA3newDMwLT,fullName("byVTightIsolationMVA3newDMwLT"));
e.put(fTtauIDbyVVTightIsolationMVA3newDMwLT,fullName("byVVTightIsolationMVA3newDMwLT"));
e.put(fTtauIDagainstElectronLoose,fullName("againstElectronLoose"));
e.put(fTtauIDagainstElectronMedium,fullName("againstElectronMedium"));
e.put(fTtauIDagainstElectronTight,fullName("againstElectronTight"));
e.put(fTtauIDagainstElectronMVA5raw,fullName("againstElectronMVA5raw"));
e.put(fTtauIDagainstElectronMVA5category,fullName("againstElectronMVA5category"));
e.put(fTtauIDagainstElectronVLooseMVA5,fullName("againstElectronVLooseMVA5"));
e.put(fTtauIDagainstElectronLooseMVA5,fullName("againstElectronLooseMVA5"));
e.put(fTtauIDagainstElectronMediumMVA5,fullName("againstElectronMediumMVA5"));
e.put(fTtauIDagainstElectronTightMVA5,fullName("againstElectronTightMVA5"));
e.put(fTtauIDagainstElectronVTightMVA5,fullName("againstElectronVTightMVA5"));
e.put(fTtauIDagainstElectronDeadECAL,fullName("againstElectronDeadECAL"));
e.put(fTtauIDagainstMuonLoose,fullName("againstMuonLoose"));
e.put(fTtauIDagainstMuonMedium,fullName("againstMuonMedium"));
e.put(fTtauIDagainstMuonTight,fullName("againstMuonTight"));
e.put(fTtauIDagainstMuonLoose2,fullName("againstMuonLoose2"));
e.put(fTtauIDagainstMuonMedium2,fullName("againstMuonMedium2"));
e.put(fTtauIDagainstMuonTight2,fullName("againstMuonTight2"));
e.put(fTtauIDagainstMuonLoose3,fullName("againstMuonLoose3"));
e.put(fTtauIDagainstMuonTight3,fullName("againstMuonTight3"));
e.put(fTtauIDagainstMuonMVAraw,fullName("againstMuonMVAraw"));
e.put(fTtauIDagainstMuonLooseMVA,fullName("againstMuonLooseMVA"));
e.put(fTtauIDagainstMuonMediumMVA,fullName("againstMuonMediumMVA"));
e.put(fTtauIDagainstMuonTightMVA,fullName("againstMuonTightMVA"));
//END


  }else if(fType == El){
    e.put(fTID95,fullName("ID95"));
    e.put(fTID90,fullName("ID90"));
    e.put(fTID85,fullName("ID85"));
    e.put(fTID80,fullName("ID80"));
  }else if(fType == Mu){
    e.put(fTPtErr,     fullName("PtErr"));
    e.put(fTMuNMatches,fullName("NMatches"));
  }

}

//______________________________________________________________________________
template <class LeptonType>
void LeptonFillerPat<LeptonType>::resetProducts(void) {

  fTMaxLepExc.reset(new int(0));
  fTNObjsTot .reset(new int(0));
  fTNObjs    .reset(new int(0));

  // Reset all arrays
  fTPx.reset(new std::vector<float>);
  fTPy.reset(new std::vector<float>);
  fTPz.reset(new std::vector<float>);
  fTPt.reset(new std::vector<float>);
  fTEta.reset(new std::vector<float>);
  fTPhi.reset(new std::vector<float>);
  fTE.reset(new std::vector<float>);
  fTEt.reset(new std::vector<float>);
  fTCharge.reset(new std::vector<int>);

  fTParticleIso.reset(new std::vector<float>);
  fTChargedHadronIso.reset(new std::vector<float>);
  fTNeutralHadronIso.reset(new std::vector<float>);
  fTPhotonIso.reset(new std::vector<float>);

  if(fType == Tau){
    fTDecayMode.reset(new std::vector<int>);
    fTIsPFTau.reset(new std::vector<int>);
    fTVz.reset(new std::vector<float>); 
    fTEmFraction.reset(new std::vector<float>); 
    fTJetPt.reset(new std::vector<float>);
    fTJetEta.reset(new std::vector<float>);
    fTJetPhi.reset(new std::vector<float>);
    fTJetMass.reset(new std::vector<float>);
    fTLeadingTkPt.reset(new std::vector<float>);
    fTLeadingNeuPt.reset(new std::vector<float>);
    fTLeadingTkHcalenergy.reset(new std::vector<float>);
    fTLeadingTkEcalenergy.reset(new std::vector<float>);
    fTNumChargedHadronsSignalCone.reset(new std::vector<int>);
    fTNumNeutralHadronsSignalCone.reset(new std::vector<int>);
    fTNumPhotonsSignalCone.reset(new std::vector<int>);
    fTNumParticlesSignalCone.reset(new std::vector<int>);
    fTNumChargedHadronsIsoCone.reset(new std::vector<int>);
    fTNumNeutralHadronsIsoCone.reset(new std::vector<int>);
    fTNumPhotonsIsolationCone.reset(new std::vector<int>);
    fTNumParticlesIsolationCone.reset(new std::vector<int>);
    fTPtSumChargedParticlesIsoCone.reset(new std::vector<float>);
    fTPtSumPhotonsIsoCone.reset(new std::vector<float>);

//tauid resets, by Hamed
fTtauIDdecayModeFindingNewDMs.reset(new std::vector<float>);
fTtauIDdecayModeFindingOldDMs.reset(new std::vector<float>);
fTtauIDdecayModeFinding.reset(new std::vector<float>);
fTtauIDbyLooseIsolation.reset(new std::vector<float>);
fTtauIDbyVLooseCombinedIsolationDeltaBetaCorr.reset(new std::vector<float>);
fTtauIDbyLooseCombinedIsolationDeltaBetaCorr.reset(new std::vector<float>);
fTtauIDbyMediumCombinedIsolationDeltaBetaCorr.reset(new std::vector<float>);
fTtauIDbyTightCombinedIsolationDeltaBetaCorr.reset(new std::vector<float>);
fTtauIDbyCombinedIsolationDeltaBetaCorrRaw.reset(new std::vector<float>);
fTtauIDbyLooseCombinedIsolationDeltaBetaCorr3Hits.reset(new std::vector<float>);
fTtauIDbyMediumCombinedIsolationDeltaBetaCorr3Hits.reset(new std::vector<float>);
fTtauIDbyTightCombinedIsolationDeltaBetaCorr3Hits.reset(new std::vector<float>);
fTtauIDbyCombinedIsolationDeltaBetaCorrRaw3Hits.reset(new std::vector<float>);
fTtauIDchargedIsoPtSum.reset(new std::vector<float>);
fTtauIDneutralIsoPtSum.reset(new std::vector<float>);
fTtauIDpuCorrPtSum.reset(new std::vector<float>);
fTtauIDbyIsolationMVA3oldDMwoLTraw.reset(new std::vector<float>);
fTtauIDbyVLooseIsolationMVA3oldDMwoLT.reset(new std::vector<float>);
fTtauIDbyLooseIsolationMVA3oldDMwoLT.reset(new std::vector<float>);
fTtauIDbyMediumIsolationMVA3oldDMwoLT.reset(new std::vector<float>);
fTtauIDbyTightIsolationMVA3oldDMwoLT.reset(new std::vector<float>);
fTtauIDbyVTightIsolationMVA3oldDMwoLT.reset(new std::vector<float>);
fTtauIDbyVVTightIsolationMVA3oldDMwoLT.reset(new std::vector<float>);
fTtauIDbyIsolationMVA3oldDMwLTraw.reset(new std::vector<float>);
fTtauIDbyVLooseIsolationMVA3oldDMwLT.reset(new std::vector<float>);
fTtauIDbyLooseIsolationMVA3oldDMwLT.reset(new std::vector<float>);
fTtauIDbyMediumIsolationMVA3oldDMwLT.reset(new std::vector<float>);
fTtauIDbyTightIsolationMVA3oldDMwLT.reset(new std::vector<float>);
fTtauIDbyVTightIsolationMVA3oldDMwLT.reset(new std::vector<float>);
fTtauIDbyVVTightIsolationMVA3oldDMwLT.reset(new std::vector<float>);
fTtauIDbyIsolationMVA3newDMwoLTraw.reset(new std::vector<float>);
fTtauIDbyVLooseIsolationMVA3newDMwoLT.reset(new std::vector<float>);
fTtauIDbyLooseIsolationMVA3newDMwoLT.reset(new std::vector<float>);
fTtauIDbyMediumIsolationMVA3newDMwoLT.reset(new std::vector<float>);
fTtauIDbyTightIsolationMVA3newDMwoLT.reset(new std::vector<float>);
fTtauIDbyVTightIsolationMVA3newDMwoLT.reset(new std::vector<float>);
fTtauIDbyVVTightIsolationMVA3newDMwoLT.reset(new std::vector<float>);
fTtauIDbyIsolationMVA3newDMwLTraw.reset(new std::vector<float>);
fTtauIDbyVLooseIsolationMVA3newDMwLT.reset(new std::vector<float>);
fTtauIDbyLooseIsolationMVA3newDMwLT.reset(new std::vector<float>);
fTtauIDbyMediumIsolationMVA3newDMwLT.reset(new std::vector<float>);
fTtauIDbyTightIsolationMVA3newDMwLT.reset(new std::vector<float>);
fTtauIDbyVTightIsolationMVA3newDMwLT.reset(new std::vector<float>);
fTtauIDbyVVTightIsolationMVA3newDMwLT.reset(new std::vector<float>);
fTtauIDagainstElectronLoose.reset(new std::vector<float>);
fTtauIDagainstElectronMedium.reset(new std::vector<float>);
fTtauIDagainstElectronTight.reset(new std::vector<float>);
fTtauIDagainstElectronMVA5raw.reset(new std::vector<float>);
fTtauIDagainstElectronMVA5category.reset(new std::vector<float>);
fTtauIDagainstElectronVLooseMVA5.reset(new std::vector<float>);
fTtauIDagainstElectronLooseMVA5.reset(new std::vector<float>);
fTtauIDagainstElectronMediumMVA5.reset(new std::vector<float>);
fTtauIDagainstElectronTightMVA5.reset(new std::vector<float>);
fTtauIDagainstElectronVTightMVA5.reset(new std::vector<float>);
fTtauIDagainstElectronDeadECAL.reset(new std::vector<float>);
fTtauIDagainstMuonLoose.reset(new std::vector<float>);
fTtauIDagainstMuonMedium.reset(new std::vector<float>);
fTtauIDagainstMuonTight.reset(new std::vector<float>);
fTtauIDagainstMuonLoose2.reset(new std::vector<float>);
fTtauIDagainstMuonMedium2.reset(new std::vector<float>);
fTtauIDagainstMuonTight2.reset(new std::vector<float>);
fTtauIDagainstMuonLoose3.reset(new std::vector<float>);
fTtauIDagainstMuonTight3.reset(new std::vector<float>);
fTtauIDagainstMuonMVAraw.reset(new std::vector<float>);
fTtauIDagainstMuonLooseMVA.reset(new std::vector<float>);
fTtauIDagainstMuonMediumMVA.reset(new std::vector<float>);
fTtauIDagainstMuonTightMVA.reset(new std::vector<float>);
//END
  }else if(fType == El){
    fTID95.reset(new std::vector<int>);
    fTID90.reset(new std::vector<int>);
    fTID85.reset(new std::vector<int>);
    fTID80.reset(new std::vector<int>);
  }else if(fType == Mu){
    fTMuNMatches.reset(new std::vector<int>);
    fTPtErr.reset(new std::vector<float>);
  }
  

}

//________________________________________________________________________________________
template <>
void LeptonFillerPat<pat::Tau>::getSpecific(const pat::Tau& lepton){
  // speficic for PFTaus
  fTIsPFTau  ->push_back( lepton.isPFTau() );	
  fTDecayMode  ->push_back( lepton.decayMode() );	
  fTVz         ->push_back( lepton.vz() );  
  fTEmFraction ->push_back( lepton.emFraction() ); 
  fTJetPt      ->push_back( lepton.pfJetRef().get()->pt() );
  fTJetEta     ->push_back( lepton.pfJetRef().get()->eta() );
  fTJetPhi     ->push_back( lepton.pfJetRef().get()->phi() );
  fTJetMass    ->push_back( lepton.pfJetRef().get()->mass() );

  fTLeadingTkPt        ->push_back( (lepton.leadPFChargedHadrCand())->pt() );
  fTLeadingNeuPt       ->push_back( (lepton.leadPFNeutralCand().isNonnull() ? lepton.leadPFNeutralCand()->pt() : 0.) );
  fTLeadingTkHcalenergy->push_back( lepton.leadPFChargedHadrCand()->hcalEnergy() );
  fTLeadingTkEcalenergy->push_back( lepton.leadPFChargedHadrCand()->ecalEnergy() );

  fTNumChargedHadronsSignalCone->push_back( lepton.signalPFChargedHadrCands().size() );
  fTNumNeutralHadronsSignalCone->push_back( lepton.signalPFNeutrHadrCands().size() );
  fTNumPhotonsSignalCone->push_back( lepton.signalPFGammaCands().size() );
  fTNumParticlesSignalCone->push_back( lepton.signalPFCands().size() );

  fTNumChargedHadronsIsoCone->push_back( lepton.isolationPFChargedHadrCands().size() );
  fTNumNeutralHadronsIsoCone->push_back( lepton.isolationPFNeutrHadrCands().size() );
  fTNumPhotonsIsolationCone->push_back( lepton.isolationPFGammaCands().size() );
  fTNumParticlesIsolationCone->push_back( lepton.isolationPFCands().size() );
  fTPtSumChargedParticlesIsoCone->push_back( lepton.isolationPFChargedHadrCandsPtSum() );
  fTPtSumPhotonsIsoCone->push_back( lepton.isolationPFGammaCandsEtSum() );


//tauid push_backs, by Hamed
fTtauIDdecayModeFindingNewDMs->push_back( lepton.tauID("decayModeFindingNewDMs") );
fTtauIDdecayModeFindingOldDMs->push_back( lepton.tauID("decayModeFindingOldDMs") );
fTtauIDdecayModeFinding->push_back( lepton.tauID("decayModeFinding") );
fTtauIDbyLooseIsolation->push_back( lepton.tauID("byLooseIsolation") );
fTtauIDbyVLooseCombinedIsolationDeltaBetaCorr->push_back( lepton.tauID("byVLooseCombinedIsolationDeltaBetaCorr") );
fTtauIDbyLooseCombinedIsolationDeltaBetaCorr->push_back( lepton.tauID("byLooseCombinedIsolationDeltaBetaCorr") );
fTtauIDbyMediumCombinedIsolationDeltaBetaCorr->push_back( lepton.tauID("byMediumCombinedIsolationDeltaBetaCorr") );
fTtauIDbyTightCombinedIsolationDeltaBetaCorr->push_back( lepton.tauID("byTightCombinedIsolationDeltaBetaCorr") );
fTtauIDbyCombinedIsolationDeltaBetaCorrRaw->push_back( lepton.tauID("byCombinedIsolationDeltaBetaCorrRaw") );
fTtauIDbyLooseCombinedIsolationDeltaBetaCorr3Hits->push_back( lepton.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") );
fTtauIDbyMediumCombinedIsolationDeltaBetaCorr3Hits->push_back( lepton.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits") );
fTtauIDbyTightCombinedIsolationDeltaBetaCorr3Hits->push_back( lepton.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits") );
fTtauIDbyCombinedIsolationDeltaBetaCorrRaw3Hits->push_back( lepton.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") );
fTtauIDchargedIsoPtSum->push_back( lepton.tauID("chargedIsoPtSum") );
fTtauIDneutralIsoPtSum->push_back( lepton.tauID("neutralIsoPtSum") );
fTtauIDpuCorrPtSum->push_back( lepton.tauID("puCorrPtSum") );
fTtauIDbyIsolationMVA3oldDMwoLTraw->push_back( lepton.tauID("byIsolationMVA3oldDMwoLTraw") );
fTtauIDbyVLooseIsolationMVA3oldDMwoLT->push_back( lepton.tauID("byVLooseIsolationMVA3oldDMwoLT") );
fTtauIDbyLooseIsolationMVA3oldDMwoLT->push_back( lepton.tauID("byLooseIsolationMVA3oldDMwoLT") );
fTtauIDbyMediumIsolationMVA3oldDMwoLT->push_back( lepton.tauID("byMediumIsolationMVA3oldDMwoLT") );
fTtauIDbyTightIsolationMVA3oldDMwoLT->push_back( lepton.tauID("byTightIsolationMVA3oldDMwoLT") );
fTtauIDbyVTightIsolationMVA3oldDMwoLT->push_back( lepton.tauID("byVTightIsolationMVA3oldDMwoLT") );
fTtauIDbyVVTightIsolationMVA3oldDMwoLT->push_back( lepton.tauID("byVVTightIsolationMVA3oldDMwoLT") );
fTtauIDbyIsolationMVA3oldDMwLTraw->push_back( lepton.tauID("byIsolationMVA3oldDMwLTraw") );
fTtauIDbyVLooseIsolationMVA3oldDMwLT->push_back( lepton.tauID("byVLooseIsolationMVA3oldDMwLT") );
fTtauIDbyLooseIsolationMVA3oldDMwLT->push_back( lepton.tauID("byLooseIsolationMVA3oldDMwLT") );
fTtauIDbyMediumIsolationMVA3oldDMwLT->push_back( lepton.tauID("byMediumIsolationMVA3oldDMwLT") );
fTtauIDbyTightIsolationMVA3oldDMwLT->push_back( lepton.tauID("byTightIsolationMVA3oldDMwLT") );
fTtauIDbyVTightIsolationMVA3oldDMwLT->push_back( lepton.tauID("byVTightIsolationMVA3oldDMwLT") );
fTtauIDbyVVTightIsolationMVA3oldDMwLT->push_back( lepton.tauID("byVVTightIsolationMVA3oldDMwLT") );
fTtauIDbyIsolationMVA3newDMwoLTraw->push_back( lepton.tauID("byIsolationMVA3newDMwoLTraw") );
fTtauIDbyVLooseIsolationMVA3newDMwoLT->push_back( lepton.tauID("byVLooseIsolationMVA3newDMwoLT") );
fTtauIDbyLooseIsolationMVA3newDMwoLT->push_back( lepton.tauID("byLooseIsolationMVA3newDMwoLT") );
fTtauIDbyMediumIsolationMVA3newDMwoLT->push_back( lepton.tauID("byMediumIsolationMVA3newDMwoLT") );
fTtauIDbyTightIsolationMVA3newDMwoLT->push_back( lepton.tauID("byTightIsolationMVA3newDMwoLT") );
fTtauIDbyVTightIsolationMVA3newDMwoLT->push_back( lepton.tauID("byVTightIsolationMVA3newDMwoLT") );
fTtauIDbyVVTightIsolationMVA3newDMwoLT->push_back( lepton.tauID("byVVTightIsolationMVA3newDMwoLT") );
fTtauIDbyIsolationMVA3newDMwLTraw->push_back( lepton.tauID("byIsolationMVA3newDMwLTraw") );
fTtauIDbyVLooseIsolationMVA3newDMwLT->push_back( lepton.tauID("byVLooseIsolationMVA3newDMwLT") );
fTtauIDbyLooseIsolationMVA3newDMwLT->push_back( lepton.tauID("byLooseIsolationMVA3newDMwLT") );
fTtauIDbyMediumIsolationMVA3newDMwLT->push_back( lepton.tauID("byMediumIsolationMVA3newDMwLT") );
fTtauIDbyTightIsolationMVA3newDMwLT->push_back( lepton.tauID("byTightIsolationMVA3newDMwLT") );
fTtauIDbyVTightIsolationMVA3newDMwLT->push_back( lepton.tauID("byVTightIsolationMVA3newDMwLT") );
fTtauIDbyVVTightIsolationMVA3newDMwLT->push_back( lepton.tauID("byVVTightIsolationMVA3newDMwLT") );
fTtauIDagainstElectronLoose->push_back( lepton.tauID("againstElectronLoose") );
fTtauIDagainstElectronMedium->push_back( lepton.tauID("againstElectronMedium") );
fTtauIDagainstElectronTight->push_back( lepton.tauID("againstElectronTight") );
fTtauIDagainstElectronMVA5raw->push_back( lepton.tauID("againstElectronMVA5raw") );
fTtauIDagainstElectronMVA5category->push_back( lepton.tauID("againstElectronMVA5category") );
fTtauIDagainstElectronVLooseMVA5->push_back( lepton.tauID("againstElectronVLooseMVA5") );
fTtauIDagainstElectronLooseMVA5->push_back( lepton.tauID("againstElectronLooseMVA5") );
fTtauIDagainstElectronMediumMVA5->push_back( lepton.tauID("againstElectronMediumMVA5") );
fTtauIDagainstElectronTightMVA5->push_back( lepton.tauID("againstElectronTightMVA5") );
fTtauIDagainstElectronVTightMVA5->push_back( lepton.tauID("againstElectronVTightMVA5") );
fTtauIDagainstElectronDeadECAL->push_back( lepton.tauID("againstElectronDeadECAL") );
fTtauIDagainstMuonLoose->push_back( lepton.tauID("againstMuonLoose") );
fTtauIDagainstMuonMedium->push_back( lepton.tauID("againstMuonMedium") );
fTtauIDagainstMuonTight->push_back( lepton.tauID("againstMuonTight") );
fTtauIDagainstMuonLoose2->push_back( lepton.tauID("againstMuonLoose2") );
fTtauIDagainstMuonMedium2->push_back( lepton.tauID("againstMuonMedium2") );
fTtauIDagainstMuonTight2->push_back( lepton.tauID("againstMuonTight2") );
fTtauIDagainstMuonLoose3->push_back( lepton.tauID("againstMuonLoose3") );
fTtauIDagainstMuonTight3->push_back( lepton.tauID("againstMuonTight3") );
fTtauIDagainstMuonMVAraw->push_back( lepton.tauID("againstMuonMVAraw") );
fTtauIDagainstMuonLooseMVA->push_back( lepton.tauID("againstMuonLooseMVA") );
fTtauIDagainstMuonMediumMVA->push_back( lepton.tauID("againstMuonMediumMVA") );
fTtauIDagainstMuonTightMVA->push_back( lepton.tauID("againstMuonTightMVA") );
//END

  return;
}

template <>
void LeptonFillerPat<pat::Electron>::getSpecific(const pat::Electron& lepton){
  // speficic for PFElectrons
  fTID95->push_back( lepton.electronID("simpleEleId95cIso") );
  fTID90->push_back( lepton.electronID("simpleEleId90cIso") );
  fTID85->push_back( lepton.electronID("simpleEleId85cIso") );
  fTID80->push_back( lepton.electronID("simpleEleId80cIso") );
  return;
}

template <>
void LeptonFillerPat<pat::Muon>::getSpecific(const pat::Muon& lepton){
  // speficic for PFMuon
  fTMuNMatches->push_back( lepton.numberOfMatches() );
  fTPtErr     ->push_back( lepton.globalTrack()->ptError() );
  return;
}

#endif
