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
#include <algorithm> 

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

      TH1F * m_gen_photon_pt1;
      TH1F * m_gen_photon_pt2;
      TH1F * m_gen_photon_eta;
      TH1F * m_gen_photon_type;
      TH1F * m_gen_photon_dr;
      TH1F * m_gen_photon_dphi;
      TH1F * m_gen_photon_deta;

      TH1F * m_gen_filt_diphoton_mass;

      TH1F * m_gen_filt_photon_pt1;
      TH1F * m_gen_filt_photon_pt2;
      TH1F * m_gen_filt_photon_eta;
      TH1F * m_gen_filt_photon_dr;
      TH1F * m_gen_filt_photon_dphi;
      TH1F * m_gen_filt_photon_deta;

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
    m_gen_photon_pt1           =  fs->make<TH1F>("gen_photon_pt1",           ";p_{T1};Photons",            20,          0,    125);
    m_gen_photon_pt2           =  fs->make<TH1F>("gen_photon_pt2",           ";p_{T2};Photons",            20,          0,    125);
    m_gen_photon_eta          =  fs->make<TH1F>("gen_photon_eta",          ";#eta;Photons",             20,          -4,   4);
    m_gen_photon_type         =  fs->make<TH1F>("gen_photon_type",         ";Event type;Events",       4,           0,    4);
    m_gen_photon_dr           =  fs->make<TH1F>("gen_photon_dr",           ";#Delta R;Photons",         20,          0,    10);
    m_gen_photon_dphi           =  fs->make<TH1F>("gen_photon_dphi",           ";#Delta #phi;Photons",         20,          0,    4);
    m_gen_photon_deta           =  fs->make<TH1F>("gen_photon_deta",           ";#Delta #eta;Photons",         20,          0,    10);

    m_gen_filt_photon_pt1      =  fs->make<TH1F>("gen_filt_photon_pt1",      ";p_{T1};Photons",            20,          0,    125);
    m_gen_filt_photon_pt2      =  fs->make<TH1F>("gen_filt_photon_pt2",      ";p_{T2};Photons",            20,          0,    125);
    m_gen_filt_photon_eta     =  fs->make<TH1F>("gen_filt_photon_eta",     ";#eta;Photons",             20,          -4,   4);
    m_gen_filt_photon_dr      =  fs->make<TH1F>("gen_filt_photon_dr",      ";#Delta R;Photons",         20,          0,    10);
    m_gen_filt_diphoton_mass  =  fs->make<TH1F>("gen_filt_diphoton_mass",  ";M_{#gamma#gamma};Events",  20,          123,  133);
    m_gen_filt_photon_dphi    =  fs->make<TH1F>("gen_filt_photon_dphi",    ";#Delta #phi;Photons",         20,          0,    4);
    m_gen_filt_photon_deta    =  fs->make<TH1F>("gen_filt_photon_deta",    ";#Delta #eta;Photons",         20,          0,    10);

    m_diphoton_mass[0.05f]    =  fs->make<TH1F>("diphoton_mass05",         ";M_{#gamma#gamma};Events",  20,          100,  150);
    m_diphoton_mass[0.01f]    =  fs->make<TH1F>("diphoton_mass01",         ";M_{#gamma#gamma};Events",  20,          100,  150);
    m_diphoton_mass[0.02f]    =  fs->make<TH1F>("diphoton_mass02",         ";M_{#gamma#gamma};Events",  20,          100,  150);
    
    m_gen_photon_pt_cut = 10.0;

    m_gen_photon_type->GetXaxis()->SetBinLabel(1,"None");
    m_gen_photon_type->GetXaxis()->SetBinLabel(2,"B+B");
    m_gen_photon_type->GetXaxis()->SetBinLabel(3,"E+E");
    m_gen_photon_type->GetXaxis()->SetBinLabel(4,"E+B");
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

   std::cout << "Starting analyzer\n";
   std::cout << "Checking " << genParticles->size() << " genParticles\n";
   for( edm::View<reco::GenParticle>::const_iterator genParticle =  genParticles->begin();
           genParticle != genParticles->end();
           ++genParticle)
   {
       //if(genParticle->pdgId() == 25)
       //{
       //    std::cout << "Found a Higgs with " << genParticle->numberOfDaughters() << " daughters\n";

       //    for(size_t i=0; i < genParticle->numberOfDaughters(); ++i)
       //    {
       //        const reco::GenParticle * daughter = dynamic_cast<const reco::GenParticle *>(genParticle->daughter(i));
       //        std::cout << "Daughter " << i << " has ID " << daughter->pdgId() << std::endl;
       //        std::cout << "Daughter " << i << " has status " << daughter->status() << std::endl;
       //        std::cout << "Daughter has " << daughter->numberOfDaughters() << " daughters\n";

       //        for(size_t j=0; j < daughter->numberOfDaughters(); ++j)
       //        {
       //            const reco::GenParticle * sub_daughter = dynamic_cast<const reco::GenParticle *>(daughter->daughter(j));
       //            std::cout << "sub Daughter " << j << " has ID " << sub_daughter->pdgId() << std::endl;
       //            std::cout << "sub Daughter " << j << " has status " << sub_daughter->status() << std::endl;
       //        }
       //    }
       //}

       if(genParticle->pdgId() == 25 && genParticle->numberOfDaughters() == 2) //Grab the Higgs that decays to 2 particles
       {
           std::cout << "Found Higgs with 2(3?) daughters\n";
           const reco::GenParticle * photon0 = dynamic_cast<const reco::GenParticle *>(genParticle->daughter(0));
           const reco::GenParticle * photon1 = dynamic_cast<const reco::GenParticle *>(genParticle->daughter(1));

           if( photon0->pdgId() != 22 || photon0->pdgId() != 22) // Not Higgs to gg
           {
               std::cout << "daughters are not photons\n";
               continue;
           }

           type0 = get_photonType(photon0->eta());
           type1 = get_photonType(photon1->eta());

           m_gen_photon_pt1->Fill(std::max(photon0->pt(),photon1->pt()));
           m_gen_photon_pt2->Fill(std::min(photon0->pt(),photon1->pt()));

           m_gen_photon_eta->Fill(photon0->eta());
           m_gen_photon_eta->Fill(photon1->eta());

           if(type0 == BARREL_PHOTON)
               if(type1 == BARREL_PHOTON)
                   m_gen_photon_type->Fill(1);
               else if(type1 == ENDCAP_PHOTON)
                   m_gen_photon_type->Fill(3);
               else
                   m_gen_photon_type->Fill(0);
           else if(type0 == ENDCAP_PHOTON)
               if(type1 == BARREL_PHOTON)
                   m_gen_photon_type->Fill(3);
               else if(type1 == ENDCAP_PHOTON)
                   m_gen_photon_type->Fill(2);
               else
                   m_gen_photon_type->Fill(0);
           else
               m_gen_photon_type->Fill(0);

           m_gen_photon_dr->Fill(deltaR(*photon0,*photon1));
           m_gen_photon_dphi->Fill(fabs(photon0->phi() - photon1->phi()));
           m_gen_photon_deta->Fill(fabs(photon0->eta() - photon1->eta()));

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
           std::cout << "Found photons\n";
           break; //Found photons, leave the loop
       }
   }

   if(!barrel_photon || !endcap_photon)
   {
       std::cout << "Did not find Barrel or Endcap photon\n";
       return;
   }
   m_gen_filt_photon_pt1->Fill(std::max(barrel_photon->pt(),endcap_photon->pt()));
   m_gen_filt_photon_pt2->Fill(std::min(barrel_photon->pt(),endcap_photon->pt()));

   m_gen_filt_photon_eta->Fill(barrel_photon->eta());
   m_gen_filt_photon_eta->Fill(endcap_photon->eta());

   m_gen_filt_photon_dr->Fill(deltaR(*barrel_photon,*endcap_photon));
   m_gen_filt_photon_dphi->Fill(fabs(barrel_photon->phi() - endcap_photon->phi()));
   m_gen_filt_photon_deta->Fill(fabs(barrel_photon->eta() - endcap_photon->eta()));

   m_gen_filt_diphoton_mass->Fill( (barrel_photon->p4() + endcap_photon->p4()).mass() );

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
    else if(fabs(eta) < 3.0 && fabs(eta) > 1.6) return ENDCAP_PHOTON;
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
