// -*- C++ -*-
//
// Package:    HiggsAnalyzer
// Class:      HiggsAnalyzer
// 
/**\class HiggsAnalyzer HiggsAnalyzer.cc UserCode/HiggsAnalyzer/plugins/HiggsAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Peter Hansen
//         Created:  Thu, 11 Sep 2014 19:53:37 GMT
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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TH1F.h"
#include "TH2F.h"

//
// class declaration
//

class HiggsAnalyzer : public edm::EDAnalyzer {
   public:
      explicit HiggsAnalyzer(const edm::ParameterSet&);
      ~HiggsAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      std::string m_genSource;
      std::string m_jetSource;

      TH1F * m_gen_photon_pt;
      TH1F * m_gen_photon_eta;
      TH1F * m_gen_photon_type;
      TH1F * m_gen_diphoton_mass;
      std::map<float,TH1F*> m_diphoton_mass;

      float m_gen_photon_pt_cut;

      enum photonType
      {
          NO_PHOTON,
          BARREL_PHOTON,
          ENDCAP_PHOTON
      };

      photonType get_photonType(float eta);
      template<class R>
      bool is_matched(const reco::GenParticle &gen, const R &rec, float dr2_limit);
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
HiggsAnalyzer::HiggsAnalyzer(const edm::ParameterSet& iConfig)

{
    //now do what ever initialization is needed
    m_genSource        = iConfig.getUntrackedParameter<std::string>("genSource");
    m_jetSource        = iConfig.getUntrackedParameter<std::string>("jetSource");

    edm::Service<TFileService> fs;
    m_gen_photon_pt      =  fs->make<TH1F>("gen_photon_pt",      ";p_{T};Photons",            20,  0,    125);
    m_gen_photon_eta     =  fs->make<TH1F>("gen_photon_eta",     ";#eta;Photons",             20,  -4,   4);
    m_gen_photon_type     =  fs->make<TH1F>("gen_photon_type",    ";type;Photons",             3,   0,    3);
    m_gen_diphoton_mass  =  fs->make<TH1F>("gen_diphoton_mass",  ";M_{#gamma#gamma};Events",  20,  123,  127);
        
    m_diphoton_mass[0.05f] = fs->make<TH1F>("diphoton_mass05",  ";M_{#gamma#gamma};Events",  20,  100,  150);
    m_diphoton_mass[0.1f] = fs->make<TH1F>("diphoton_mass1",  ";M_{#gamma#gamma};Events",  20,  100,  150);
    m_diphoton_mass[0.2f] = fs->make<TH1F>("diphoton_mass2",  ";M_{#gamma#gamma};Events",  20,  100,  150);
    
    m_gen_photon_pt_cut = 10.0;

    m_gen_photon_type->GetXaxis()->SetBinLabel(1,"None");
    m_gen_photon_type->GetXaxis()->SetBinLabel(2,"Barrel");
    m_gen_photon_type->GetXaxis()->SetBinLabel(3,"Endcap");

}


HiggsAnalyzer::~HiggsAnalyzer()
{
 
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)


}


//
// member functions
//

// ------------ method called for each event  ------------
void
HiggsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   edm::Handle<edm::View<reco::GenParticle> > genParticles;
   iEvent.getByLabel(edm::InputTag(m_genSource), genParticles);

   const reco::GenParticle * barrel_photon = NULL;
   const reco::GenParticle * endcap_photon = NULL;
   photonType type0 = NO_PHOTON;
   photonType type1 = NO_PHOTON;

   for( edm::View<reco::GenParticle>::const_iterator genParticle =  genParticles->begin();
           genParticle != genParticles->end();
           ++genParticle)
   {
       if(genParticle->pdgId() == 25 && genParticle->numberOfDaughters() == 2) //Grab the Higgs that decays to 2 particles
       {
           const reco::GenParticle * photon0 = dynamic_cast<const reco::GenParticle *>(genParticle->daughter(0));
           const reco::GenParticle * photon1 = dynamic_cast<const reco::GenParticle *>(genParticle->daughter(1));

           if( photon0->pdgId() != 22 || photon0->pdgId() != 22) // Not Higgs to gg
               continue;

           type0 = get_photonType(photon0->eta());
           type1 = get_photonType(photon1->eta());

           if(type0 == ENDCAP_PHOTON && type1 == BARREL_PHOTON)
           {
               endcap_photon = photon0;
               barrel_photon = photon1;
           }
           else if(type0 == BARREL_PHOTON && type1 == ENDCAP_PHOTON)
           {
               barrel_photon = photon0;
               endcap_photon = photon1;
           }
           break; //Found photons, leave the loop
       }
   }

   if(!barrel_photon || !endcap_photon)
       return;
   m_gen_photon_pt->Fill(barrel_photon->pt());
   m_gen_photon_pt->Fill(endcap_photon->pt());

   m_gen_photon_eta->Fill(barrel_photon->eta());
   m_gen_photon_eta->Fill(endcap_photon->eta());


   m_gen_photon_type->Fill(type0);
   m_gen_photon_type->Fill(type1);

   m_gen_diphoton_mass->Fill( (barrel_photon->p4() + endcap_photon->p4()).mass() );

   edm::Handle<edm::View<reco::PFCluster> > hgcee_clusters;
   iEvent.getByLabel(edm::InputTag("particleFlowClusterHGCEE"), hgcee_clusters);

   std::map<float,math::PtEtaPhiMLorentzVector> clusterp4;
   for(std::map<float,TH1F*>::iterator it = m_diphoton_mass.begin();
           it != m_diphoton_mass.end(); ++it)
   {
       clusterp4[it->first] = math::PtEtaPhiMLorentzVector(0.0f, 0.0f, 0.0f, 0.0f);
   }

   for( edm::View<reco::PFCluster>::const_iterator cluster =  hgcee_clusters->begin();
           cluster != hgcee_clusters->end();
           ++cluster)
   {
       for(std::map<float,TH1F*>::iterator it = m_diphoton_mass.begin();
               it != m_diphoton_mass.end(); ++it)
       {
           if(is_matched(*endcap_photon,*cluster,it->first))
               clusterp4[it->first] += math::PtEtaPhiMLorentzVector(cluster->pt(), cluster->eta(), cluster->phi(), 0 );
       }
   }
   if(clusterp4[0].pt())
       std::cout  << "We found it" << std::endl;
   else
   {
       std::cout  << "Failure!" ;
       std::cout << endcap_photon->eta() << std::endl;
   }

   edm::Handle<edm::View<reco::Photon> > reco_photons;
   iEvent.getByLabel(edm::InputTag("gedPhotons"), reco_photons);

   const reco::Photon * reco_barrel_photon = NULL;
   for( edm::View<reco::Photon>::const_iterator reco_photon =  reco_photons->begin();
           reco_photon != reco_photons->end();
           ++reco_photon)
   {
       if(is_matched(*barrel_photon,*reco_photon, 0.1))
           reco_barrel_photon = &(*reco_photon);
   }
   if(reco_barrel_photon)
       std::cout  << "We found it" << std::endl;
   else
   {
       std::cout  << "Failure!" ;
       std::cout << barrel_photon->eta() << std::endl;
   }

   for(std::map<float,TH1F*>::iterator it = m_diphoton_mass.begin();
           it != m_diphoton_mass.end(); ++it)
       if(clusterp4[it->first].pt() && reco_barrel_photon)
           it->second->Fill( (reco_barrel_photon->p4() + clusterp4[it->first]).mass() );

   //edm::Handle<edm::View<reco::PFCluster> > ecal_clusters;
   //iEvent.getByLabel(edm::InputTag("particleFlowClusterECAL"), ecal_clusters);

}

HiggsAnalyzer::photonType HiggsAnalyzer::get_photonType(float eta)
{
    if(fabs(eta) < 1.5) return BARREL_PHOTON;
    else if(fabs(eta) < 3.0) return ENDCAP_PHOTON;
    else return NO_PHOTON;
}

template<class R>
bool HiggsAnalyzer::is_matched(const reco::GenParticle &gen, const R &rec, float dr2_limit)
{
    bool ret = false;
    float dr2 = deltaR2(gen, rec);

    if(dr2 < dr2_limit)
    {
        std::cout << "Found cluster with dr2: " << dr2 << std::endl;
        ret = true;
    }
    return ret;
}

// ------------ method called once each job just before starting event loop  ------------
void 
HiggsAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HiggsAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
HiggsAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
HiggsAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
HiggsAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
HiggsAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HiggsAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HiggsAnalyzer);
