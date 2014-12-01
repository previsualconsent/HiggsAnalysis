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

      TH1F * m_gen_photon_type;

      //Histograms by Event type
      //Gen
      std::map<int, TH1F *> m_gen_photon_pt1;
      std::map<int, TH1F *> m_gen_photon_pt2;
      std::map<int, TH1F *> m_gen_photon_eta;
      std::map<int, TH1F *> m_gen_photon_dr;
      std::map<int, TH1F *> m_gen_photon_dphi;
      std::map<int, TH1F *> m_gen_photon_deta;

      //reco
      std::map<int, TH1F *> m_diphoton_mass;

      float m_gen_photon_pt_cut;

      enum photonType
      {
          NO_PHOTON,
          BARREL_PHOTON,
          ENDCAP_PHOTON
      };

      //Histograms by photon type
      //(gen-reco)/gen
      std::map<photonType, TH1F *> m_photon_pt_diff;
      std::map<photonType, TH1F *> m_photon_en_diff;

      //gen
      std::map<photonType, TH1F *> m_gen_photon_pt;
      std::map<photonType, TH1F *> m_gen_photon_en;

      //get reconstructed p4 matched to gen
      math::PtEtaPhiMLorentzVector get_p4(const edm::Event& iEvent, const reco::GenParticle & gen, photonType type);


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

    m_gen_photon_type         =  fs->make<TH1F>("gen_photon_type",         ";Event type;Events",       4,           0,    4);

    std::string names[5] = { "None", "B+B", "E+E", "E+B", "All"};
    for ( int i = 1; i < 5; ++i)
    {
      m_gen_photon_pt1 [i]  =  fs->make<TH1F>(("gen_photon_pt1"  + names[i]).c_str(), ("gen_photon_pt1"  + names[i] + ";p_{T1};Photons").c_str()         , 20, 0,   125);
      m_gen_photon_pt2 [i]  =  fs->make<TH1F>(("gen_photon_pt2"  + names[i]).c_str(), ("gen_photon_pt2"  + names[i] + ";p_{T2};Photons").c_str()         , 20, 0,   125);
      m_gen_photon_eta [i]  =  fs->make<TH1F>(("gen_photon_eta"  + names[i]).c_str(), ("gen_photon_eta"  + names[i] + ";#eta;Photons").c_str()           , 20, -4,  4);
      m_gen_photon_dr  [i]  =  fs->make<TH1F>(("gen_photon_dr"   + names[i]).c_str(), ("gen_photon_dr"   + names[i] + ";#Delta R;Photons").c_str()       , 20, 0,   10);
      m_gen_photon_dphi[i]  =  fs->make<TH1F>(("gen_photon_dphi" + names[i]).c_str(), ("gen_photon_dphi" + names[i] + ";#Delta #phi;Photons").c_str()    , 20, 0,   4);
      m_gen_photon_deta[i]  =  fs->make<TH1F>(("gen_photon_deta" + names[i]).c_str(), ("gen_photon_deta" + names[i] + ";#Delta #eta;Photons").c_str()    , 20, 0,   10);
      m_diphoton_mass  [i]  =  fs->make<TH1F>(("diphoton_mass"   + names[i]).c_str(), ("diphoton_mass"   + names[i] + ";M_{#gamma#gamma};Events").c_str(), 20, 100, 140);
    }
    m_photon_pt_diff[BARREL_PHOTON]  =  fs->make<TH1F>("pt_diff_barrel", "pt_diff_barrel;#Delta p_T;Photons" , 31, -1, 1);
    m_photon_pt_diff[ENDCAP_PHOTON]  =  fs->make<TH1F>("pt_diff_endcap", "pt_diff_endcap;#Delta p_T;Photons" , 31, -1, 1);

    m_gen_photon_pt[BARREL_PHOTON]  =  fs->make<TH1F>("gen_pt_barrel", "gen_pt_barrel;p_T;Photons" , 40, 0, 130);
    m_gen_photon_pt[ENDCAP_PHOTON]  =  fs->make<TH1F>("gen_pt_endcap", "gen_pt_endcap;p_T;Photons" , 40, 0, 130);

    m_gen_photon_en[BARREL_PHOTON]  =  fs->make<TH1F>("gen_en_barrel", "gen_en_barrel;E;Photons" , 100, 0, 300);
    m_gen_photon_en[ENDCAP_PHOTON]  =  fs->make<TH1F>("gen_en_endcap", "gen_en_endcap;E;Photons" , 100, 0, 300);

    m_photon_en_diff[BARREL_PHOTON]  =  fs->make<TH1F>("en_diff_barrel", "en_diff_barrel;#Delta E;Photons" , 31, -1, 1);
    m_photon_en_diff[ENDCAP_PHOTON]  =  fs->make<TH1F>("en_diff_endcap", "en_diff_endcap;#Delta E;Photons" , 31, -1, 1);

    m_gen_photon_pt_cut = 10.0;

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

   const reco::GenParticle * photon0 = NULL;
   const reco::GenParticle * photon1 = NULL;
   photonType type0 = NO_PHOTON;
   photonType type1 = NO_PHOTON;
   int eventType;
    

   //std::cout << "Starting analyzer\n";
   //std::cout << "Checking " << genParticles->size() << " genParticles\n";
   for( edm::View<reco::GenParticle>::const_iterator genParticle =  genParticles->begin();
           genParticle != genParticles->end();
           ++genParticle)
   {
       if(genParticle->pdgId() == 25 && iEvent.id().event() == 351)
       {
           std::cout << "Found a Higgs with " << genParticle->numberOfDaughters() << " daughters\n";

           for(size_t i=0; i < genParticle->numberOfDaughters(); ++i)
           {
               const reco::GenParticle * daughter = dynamic_cast<const reco::GenParticle *>(genParticle->daughter(i));
               std::cout << "Daughter " << i << " has ID " << daughter->pdgId() << std::endl;
               std::cout << "Daughter " << i << " has status " << daughter->status() << std::endl;
               std::cout << "Daughter has " << daughter->numberOfDaughters() << " daughters\n";

               for(size_t j=0; j < daughter->numberOfDaughters(); ++j)
               {
                   const reco::GenParticle * sub_daughter = dynamic_cast<const reco::GenParticle *>(daughter->daughter(j));
                   std::cout << "sub Daughter " << j << " has ID " << sub_daughter->pdgId() << std::endl;
                   std::cout << "sub Daughter " << j << " has status " << sub_daughter->status() << std::endl;
               }
           }
       }

       if(genParticle->pdgId() == 25 && genParticle->numberOfDaughters() == 2) //Grab the Higgs that decays to 2 particles
       {
           //std::cout << "Found Higgs with 2(3?) daughters\n";
           photon0 = dynamic_cast<const reco::GenParticle *>(genParticle->daughter(0));
           photon1 = dynamic_cast<const reco::GenParticle *>(genParticle->daughter(1));

           if( photon0->pdgId() != 22 || photon0->pdgId() != 22 || photon0->status() != 1 || photon1->status() != 1) // Not Higgs to gg
           {
               //std::cout << "daughters are not stable photons\n";
               continue;
           }

           type0 = get_photonType(photon0->eta());
           type1 = get_photonType(photon1->eta());

           m_gen_photon_pt1[4]->Fill(std::max(photon0->pt(),photon1->pt()));
           m_gen_photon_pt2[4]->Fill(std::min(photon0->pt(),photon1->pt()));

           m_gen_photon_eta[4]->Fill(photon0->eta());
           m_gen_photon_eta[4]->Fill(photon1->eta());


           if(type0 == BARREL_PHOTON)
               if(type1 == BARREL_PHOTON)
                   eventType = 1;
               else if(type1 == ENDCAP_PHOTON)
                   eventType = 3;
               else
                   eventType = 0;
           else if(type0 == ENDCAP_PHOTON)
               if(type1 == BARREL_PHOTON)
                   eventType = 3;
               else if(type1 == ENDCAP_PHOTON)
                   eventType = 2;
               else
                   eventType = 0;
           else
               eventType = 0;

           m_gen_photon_type->Fill(eventType);

           m_gen_photon_dr[4]->Fill(deltaR(*photon0,*photon1));
           m_gen_photon_dphi[4]->Fill(fabs(photon0->phi() - photon1->phi()));
           m_gen_photon_deta[4]->Fill(fabs(photon0->eta() - photon1->eta()));

           //std::cout << "Found photons\n";
           break; //Found photons, leave the loop
       }
   }

   if(type0 == NO_PHOTON || type1 == NO_PHOTON)
   {
       //std::cout << "Did not find photons in endcap/barrel\n";
       return;
   }

   m_gen_photon_pt1[eventType]->Fill(std::max(photon0->pt(),photon1->pt()));
   m_gen_photon_pt2[eventType]->Fill(std::min(photon0->pt(),photon1->pt()));

   m_gen_photon_eta[eventType]->Fill(photon0->eta());
   m_gen_photon_eta[eventType]->Fill(photon1->eta());

   m_gen_photon_dr[eventType]->Fill(deltaR(*photon0,*photon1));
   m_gen_photon_dphi[eventType]->Fill(fabs(photon0->phi() - photon1->phi()));
   m_gen_photon_deta[eventType]->Fill(fabs(photon0->eta() - photon1->eta()));


   m_gen_photon_pt[type0]->Fill(photon0->pt());
   m_gen_photon_pt[type1]->Fill(photon1->pt());

   m_gen_photon_en[type0]->Fill(photon0->energy());
   m_gen_photon_en[type1]->Fill(photon1->energy());

   math::PtEtaPhiMLorentzVector p0 = get_p4(iEvent, *photon0, type0);
   math::PtEtaPhiMLorentzVector p1= get_p4(iEvent, *photon1, type1); 


   m_diphoton_mass[eventType]->Fill( (p0 + p1).mass() );
   m_diphoton_mass[4]->Fill( (p0 + p1).mass() );

   m_photon_pt_diff[type0]->Fill( (photon0->pt() - p0.pt())/photon0->pt() );
   m_photon_pt_diff[type1]->Fill( (photon1->pt() - p1.pt())/photon1->pt() );

   m_photon_en_diff[type0]->Fill( (photon0->energy() - p0.energy())/photon0->energy() );
   m_photon_en_diff[type1]->Fill( (photon1->energy() - p1.energy())/photon1->energy() );

}

math::PtEtaPhiMLorentzVector HiggsAnalyzer::get_p4(const edm::Event& iEvent, const reco::GenParticle & gen, photonType type)
{
   math::PtEtaPhiMLorentzVector p4(0.0f, 0.0f, 0.0f, 0.0f);
   if( type == ENDCAP_PHOTON)
   {
      edm::Handle<edm::View<reco::PFCluster> > hgcee_clusters;
      iEvent.getByLabel(edm::InputTag("particleFlowClusterHGCEE"), hgcee_clusters);

      for( edm::View<reco::PFCluster>::const_iterator cluster =  hgcee_clusters->begin();
            cluster != hgcee_clusters->end();
            ++cluster)
      {
         if(is_matched(gen,*cluster,0.05))
            p4 += math::PtEtaPhiMLorentzVector(cluster->pt(), cluster->eta(), cluster->phi(), 0 );
      }
      if( (gen.pt() - p4.pt() )/gen.pt() > .2)
      {
         std::cout << "[" <<  iEvent.id().event() << "]Pt difference high at eta= " << gen.eta() << ": " << gen.pt() << ' ' << p4.pt() <<  std::endl;
      }
      return p4;
   }
   if( type == BARREL_PHOTON)
   {

      edm::Handle<edm::View<reco::Photon> > reco_photons;
      iEvent.getByLabel(edm::InputTag("gedPhotons"), reco_photons);

      for( edm::View<reco::Photon>::const_iterator reco_photon =  reco_photons->begin();
            reco_photon != reco_photons->end();
            ++reco_photon)
      {
         if(is_matched(gen,*reco_photon, 0.1))
         {
            if(p4.pt() < reco_photon->pt())
               p4 = reco_photon->p4();
         }
      }
   }
   return p4;
}

HiggsAnalyzer::photonType HiggsAnalyzer::get_photonType(float eta)
{
    if(fabs(eta) < 1.5) return BARREL_PHOTON;
    else if(fabs(eta) < 2.8 && fabs(eta) > 1.7) return ENDCAP_PHOTON;
    else return NO_PHOTON;
}

template<class R>
bool HiggsAnalyzer::is_matched(const reco::GenParticle &gen, const R &rec, float dr2_limit)
{
    bool ret = false;
    float dr2 = deltaR2(gen, rec);

    if(dr2 < dr2_limit)
    {
        //std::cout << "Found cluster with dr2: " << dr2 << std::endl;
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
