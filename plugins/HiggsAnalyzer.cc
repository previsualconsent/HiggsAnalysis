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
#include "DataFormats/VertexReco/interface/VertexFwd.h"

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
      int m_n_unmatched;
      float m_mass;
      std::string m_genSource;
      std::string m_jetSource;

      float m_linCorr;

      TH1F * m_gen_event_type;

      //Histograms by Event type and cut
      //Gen
      std::map<int, std::map<int, TH1F *> > m_gen_photon_pt1;
      std::map<int, std::map<int, TH1F *> > m_gen_photon_pt2;
      std::map<int, std::map<int, TH1F *> > m_gen_photon_eta;
      std::map<int, std::map<int, TH1F *> > m_gen_photon_dr;
      std::map<int, std::map<int, TH1F *> > m_gen_photon_dphi;
      std::map<int, std::map<int, TH1F *> > m_gen_photon_deta;

      //reco
      std::map<int, std::map<int, TH1F *> > m_diphoton_mass;

      float m_gen_photon_pt_cut;

      enum photonType
      {
          NO_PHOTON,
          BARREL_PHOTON,
          ENDCAP_PHOTON
      };

      struct hgg_event_t
      {
         const edm::Event * iEvent;
         int eventType;
         math::XYZPoint vertex;
         const reco::GenParticle * photon[2];
         photonType type[2];
         math::PtEtaPhiMLorentzVector p[2];
         float r9[2];
      } m_event_info;

      //Histograms by photon type
      //(gen-reco)/gen
      std::map<int, std::map<photonType, TH1F *> > m_photon_pt_diff;
      std::map<int, std::map<photonType, TH1F *> > m_photon_en_diff;

      std::map<int, TH1F *> m_photon_endcap_r9;

      //gen
      std::map<int, std::map<photonType, TH1F *> > m_gen_photon_pt;
      std::map<int, std::map<photonType, TH1F *> > m_gen_photon_en;
      std::map<int, std::map<photonType, TH2F *> > m_gen_vs_reco_en;

      //n vertices

      TH1F * m_n_vertices;
      TH1F * m_n_vertices_good;

      //get reconstructed p4 matched to gen
      float get_mass();
      bool r9_cut(float cut_value);
      bool pt_cut(float cut_value);
      math::PtEtaPhiMLorentzVector get_p4(const edm::Event& iEvent, const reco::GenParticle & gen, photonType type);
      const reco::SuperCluster *  get_endcap_cluster(const edm::Event& iEvent, const reco::GenParticle & gen);
      const reco::Photon * get_barrel_photon(const edm::Event& iEvent, const reco::GenParticle & gen);


      void fillCutPlots(int cuti);
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
    m_genSource  = iConfig.getUntrackedParameter<std::string>("genSource");
    m_jetSource  = iConfig.getUntrackedParameter<std::string>("jetSource");
    m_linCorr  = iConfig.getUntrackedParameter<double>("linCorr");

    m_n_unmatched = 0;

    edm::Service<TFileService> fs;

    m_gen_event_type         =  fs->make<TH1F>("gen_event_type",         ";Event type;Events",       5,           0,    5);
    m_n_vertices = fs->make<TH1F>("n_vert", ";vertices;Events",180,0,180);
    m_n_vertices_good = fs->make<TH1F>("n_vert_good", ";vertices;Events",180,0,180);

    std::string names[6] = { "All", "!E|B", "B+B", "E+E", "E+B", "E-E"};
    std::string cuts[4] = { "None","reco", "pt", "r9"};

    for ( int i = 1; i < 6; ++i)
       m_gen_event_type->GetXaxis()->SetBinLabel(i,names[i].c_str());

    for ( int i = 0; i < 6; ++i)
       for (int j = 0; j < 4; j++)
       {
          m_gen_photon_pt1 [i][j]  =  fs->make<TH1F>(("gen_photon_pt1"  + names[i] + "_" + cuts[j]).c_str(), ("gen_photon_pt1"  + names[i] + "_" + cuts[j] + ";p_{T1};Photons").c_str()         , 20, 0,   125);
          m_gen_photon_pt2 [i][j]  =  fs->make<TH1F>(("gen_photon_pt2"  + names[i] + "_" + cuts[j]).c_str(), ("gen_photon_pt2"  + names[i] + "_" + cuts[j] + ";p_{T2};Photons").c_str()         , 20, 0,   125);
          m_gen_photon_eta [i][j]  =  fs->make<TH1F>(("gen_photon_eta"  + names[i] + "_" + cuts[j]).c_str(), ("gen_photon_eta"  + names[i] + "_" + cuts[j] + ";#eta;Photons").c_str()           , 20, -4,  4);
          m_gen_photon_dr  [i][j]  =  fs->make<TH1F>(("gen_photon_dr"   + names[i] + "_" + cuts[j]).c_str(), ("gen_photon_dr"   + names[i] + "_" + cuts[j] + ";#Delta R;Photons").c_str()       , 20, 0,   10);
          m_gen_photon_dphi[i][j]  =  fs->make<TH1F>(("gen_photon_dphi" + names[i] + "_" + cuts[j]).c_str(), ("gen_photon_dphi" + names[i] + "_" + cuts[j] + ";#Delta #phi;Photons").c_str()    , 20, 0,   4);
          m_gen_photon_deta[i][j]  =  fs->make<TH1F>(("gen_photon_deta" + names[i] + "_" + cuts[j]).c_str(), ("gen_photon_deta" + names[i] + "_" + cuts[j] + ";#Delta #eta;Photons").c_str()    , 20, 0,   10);

       }

    for ( int i = 0; i < 6; ++i)
      for (int j = 0; j < 4; j++)
         m_diphoton_mass  [i][j]  =  fs->make<TH1F>(("diphoton_mass"   + names[i] + "_" + cuts[j] ).c_str(), ("diphoton_mass_"   + names[i] + "_" + cuts[j] + ";M_{#gamma#gamma};Events/GeV").c_str(), 50, 100, 150);



    for( int j = 0; j < 4; j++)
    {
       m_photon_pt_diff[j][BARREL_PHOTON]  =  fs->make<TH1F>(("pt_diff_barrel_" + cuts[j]).c_str(), "pt_diff_barrel;#Delta p_T;Photons" , 31, -1, 1);
       m_photon_pt_diff[j][ENDCAP_PHOTON]  =  fs->make<TH1F>(("pt_diff_endcap_" + cuts[j]).c_str(), "pt_diff_endcap;#Delta p_T;Photons" , 31, -1, 1);

       m_gen_photon_pt[j][BARREL_PHOTON]  =  fs->make<TH1F>(("gen_pt_barrel_" + cuts[j]).c_str(), "gen_pt_barrel;p_T;Photons" , 40, 0, 130);
       m_gen_photon_pt[j][ENDCAP_PHOTON]  =  fs->make<TH1F>(("gen_pt_endcap_" + cuts[j]).c_str(), "gen_pt_endcap;p_T;Photons" , 40, 0, 130);

       m_gen_photon_en[j][BARREL_PHOTON]  =  fs->make<TH1F>(("gen_en_barrel_" + cuts[j]).c_str(), "gen_en_barrel;E;Photons" , 100, 0, 300);
       m_gen_photon_en[j][ENDCAP_PHOTON]  =  fs->make<TH1F>(("gen_en_endcap_" + cuts[j]).c_str(), "gen_en_endcap;E;Photons" , 100, 0, 300);

       m_gen_vs_reco_en[j][BARREL_PHOTON]  =  fs->make<TH2F>(("gen_vs_roco_en_barrel_" + cuts[j]).c_str(), "gen_vs_reco_en_barrel;Gen En;Reco En" , 100, 0, 500, 100, 0, 500);
       m_gen_vs_reco_en[j][ENDCAP_PHOTON]  =  fs->make<TH2F>(("gen_vs_roco_en_endcap_" + cuts[j]).c_str(), "gen_vs_reco_en_endcap;Gen En;Reco En" , 100, 0, 500, 100, 0, 500);

       m_photon_en_diff[j][BARREL_PHOTON]  =  fs->make<TH1F>(("en_diff_barrel_" + cuts[j]).c_str(), "en_diff_barrel;#Delta E;Photons" , 31, -1, 1);
       m_photon_en_diff[j][ENDCAP_PHOTON]  =  fs->make<TH1F>(("en_diff_endcap_" + cuts[j]).c_str(), "en_diff_endcap;#Delta E;Photons" , 31, -1, 1);

       m_photon_endcap_r9[j] = fs->make<TH1F>(("photon_endcap_r9_" + cuts[j]).c_str(),"Endcap Photon 'R9';R_9;Photons", 100, 0, 1);
    }

    m_gen_photon_pt_cut = 20.0;

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
    
   edm::Handle<reco::VertexCollection> reco_vertices;
   iEvent.getByLabel("offlinePrimaryVertices", reco_vertices);
   int num = 0;
   int num_good = 0;
   for( reco::VertexCollection::const_iterator vertex =  reco_vertices->begin();
         vertex!= reco_vertices->end();
         ++vertex)
   {
      num++;
      if ( // Criteria copied from twiki
            !(vertex->isFake())
            && (vertex->ndof() > 4)
            && (fabs(vertex->z()) <= 24.0)
            && (vertex->position().Rho() <= 2.0)
         ) {
         num_good++;
      }
   }
   m_n_vertices->Fill(num);
   m_n_vertices_good->Fill(num_good);

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
           m_event_info.vertex = genParticle->vertex();

           if( photon0->pdgId() != 22 || photon0->pdgId() != 22 || photon0->status() != 1 || photon1->status() != 1) // Not Higgs to gg
           {
               //std::cout << "daughters are not stable photons\n";
               continue;
           }



           //std::cout << "Found photons\n";
           break; //Found photons, leave the loop
       }
   }

   //Fill gen plots for all events

   type0 = get_photonType(photon0->eta());
   type1 = get_photonType(photon1->eta());

    //std::string names[5] = { "All", "!E|B", "B+B", "E+E", "E+B"};
   if(type0 == BARREL_PHOTON)
      if(type1 == BARREL_PHOTON)
         eventType = 2;  //B+B
      else if(type1 == ENDCAP_PHOTON)
         eventType = 4; // B+E
      else
         eventType = 1; //!B|E
   else if(type0 == ENDCAP_PHOTON)
      if(type1 == BARREL_PHOTON)
         eventType = 4; // B+E
      else if(type1 == ENDCAP_PHOTON)
         if( photon0->eta() * photon1->eta() > 0)
            eventType = 3; // E+E
         else
            eventType = 5; // E-E
      else
         eventType = 1; //!B|E
   else
      eventType = 1; //!B|E

   m_gen_event_type->Fill(eventType-1);

   m_event_info.iEvent=&iEvent;
   m_event_info.eventType=eventType;
   m_event_info.photon[0]=photon0;
   m_event_info.photon[1]=photon1;
   m_event_info.type[0]=type0;
   m_event_info.type[1]=type1;
   m_event_info.r9[0]=-1.0f;
   m_event_info.r9[1]=-1.0f;
   m_mass = 0.0f;


   int cuti = 0;
   fillCutPlots(cuti++);

   if(type0 == NO_PHOTON || type1 == NO_PHOTON)
   {
       //std::cout << "Did not find photons in endcap/barrel\n";
       return;
   }

   m_mass = get_mass();

   if(!m_mass) return;

   fillCutPlots(cuti++);

   if( pt_cut(20.0f)) return;

   fillCutPlots(cuti++);

   if( r9_cut(.85)) return;

   fillCutPlots(cuti++);



}
void HiggsAnalyzer::fillCutPlots(int cuti)
{
   m_gen_photon_pt1[0][cuti]->Fill(std::max(m_event_info.photon[0]->pt(),m_event_info.photon[1]->pt()));
   m_gen_photon_pt2[0][cuti]->Fill(std::min(m_event_info.photon[0]->pt(),m_event_info.photon[1]->pt()));

   m_gen_photon_eta[0][cuti]->Fill(m_event_info.photon[0]->eta());
   m_gen_photon_eta[0][cuti]->Fill(m_event_info.photon[1]->eta());

   m_gen_photon_dr[0][cuti]->Fill(deltaR(*m_event_info.photon[0],*m_event_info.photon[1]));
   m_gen_photon_dphi[0][cuti]->Fill(fabs(m_event_info.photon[0]->phi() - m_event_info.photon[1]->phi()));
   m_gen_photon_deta[0][cuti]->Fill(fabs(m_event_info.photon[0]->eta() - m_event_info.photon[1]->eta()));

   int eventType = m_event_info.eventType;
   m_gen_photon_pt1[eventType][cuti]->Fill(std::max(m_event_info.photon[0]->pt(),m_event_info.photon[1]->pt()));
   m_gen_photon_pt2[eventType][cuti]->Fill(std::min(m_event_info.photon[0]->pt(),m_event_info.photon[1]->pt()));

   m_gen_photon_eta[eventType][cuti]->Fill(m_event_info.photon[0]->eta());
   m_gen_photon_eta[eventType][cuti]->Fill(m_event_info.photon[1]->eta());

   m_gen_photon_dr[eventType][cuti]->Fill(deltaR(*m_event_info.photon[0],*m_event_info.photon[1]));
   m_gen_photon_dphi[eventType][cuti]->Fill(fabs(m_event_info.photon[0]->phi() - m_event_info.photon[1]->phi()));
   m_gen_photon_deta[eventType][cuti]->Fill(fabs(m_event_info.photon[0]->eta() - m_event_info.photon[1]->eta()));

   for(int i=0; i<2; i++)
   {
      if(m_event_info.type[i] == NO_PHOTON) continue;
      m_gen_vs_reco_en[cuti][m_event_info.type[i]]->Fill(
            m_event_info.photon[i]->energy(), m_event_info.p[i].energy());

      m_photon_pt_diff[cuti][m_event_info.type[i]]->Fill(
            (m_event_info.photon[i]->pt() - m_event_info.p[i].pt())/m_event_info.photon[i]->pt() );

      m_photon_en_diff[cuti][m_event_info.type[i]]->Fill(
            (m_event_info.photon[i]->energy() - m_event_info.p[i].energy())/m_event_info.photon[i]->energy() );

      m_gen_photon_pt[cuti][m_event_info.type[i]]->Fill(m_event_info.photon[i]->pt());
      m_gen_photon_en[cuti][m_event_info.type[i]]->Fill(m_event_info.photon[i]->energy());

   }


   if(m_mass)
   {
      m_diphoton_mass[0][cuti]->Fill(m_mass);
      m_diphoton_mass[eventType][cuti]->Fill(m_mass);

   for(int i=0; i<2; i++)
      if(m_event_info.type[i] == ENDCAP_PHOTON)
         m_photon_endcap_r9[cuti]->Fill( m_event_info.r9[i] );

   }
}

float HiggsAnalyzer::get_mass()
{
   for(int i=0; i<2; i++)
   {
      if(m_event_info.type[i] == ENDCAP_PHOTON)
      {
         const reco::SuperCluster * pToCluster = get_endcap_cluster(*m_event_info.iEvent, *m_event_info.photon[i]);
         
         if(!pToCluster) 
         {
            std::cout << "no Cluster found" << std::endl;
            m_n_unmatched++;
            return 0.0f;
         }

         math::XYZPoint cluster_pos(pToCluster->position());
         cluster_pos.SetZ( cluster_pos.Z() - m_event_info.vertex.Z());
         float en = pToCluster->seed()->energy();
         if (m_linCorr) en /= m_linCorr;

         float pt = en/cosh(cluster_pos.eta());
         m_event_info.p[i] = math::PtEtaPhiMLorentzVector(pt, cluster_pos.eta(), cluster_pos.phi(), 0.0f);
         m_event_info.r9[i] =  pToCluster->seed()->energy() / pToCluster->energy();

      }
      else if(m_event_info.type[i] == BARREL_PHOTON)
      {
         const reco::Photon * pToPhoton = get_barrel_photon(*m_event_info.iEvent, *m_event_info.photon[i]);
         if(!pToPhoton) 
         {
            std::cout << "no Photon found ("<< m_event_info.iEvent->id().event() << ")" << std::endl;

            math::PtEtaPhiMLorentzVector p4(m_event_info.photon[i]->p4());
            std::cout << "(" << p4.pt() << ',' << p4.eta() << ',' << p4.phi() << ')' << std::endl;
            m_n_unmatched++;
            return 0.0f;
         }
         math::XYZPoint photon_pos(pToPhoton->caloPosition());
         photon_pos.SetZ(photon_pos.Z() - m_event_info.vertex.Z());

         m_event_info.p[i] = math::PtEtaPhiMLorentzVector( pToPhoton->energy() / cosh( photon_pos.eta() ), photon_pos.eta(), photon_pos.phi(), 0.0f);
         m_event_info.r9[i] =  1.0f;
      }

   }

   return (m_event_info.p[0] + m_event_info.p[1]).mass();
}

bool HiggsAnalyzer::r9_cut(float cut_value)
{
   //true if one r9 is less than cut value
   bool ret = m_event_info.r9[0] < cut_value || m_event_info.r9[1] < cut_value;

   if(ret && m_event_info.eventType == 2)
      std::cout << "r9cut debug:" 
         << m_event_info.r9[0] 
         << m_event_info.r9[1] 
         << m_event_info.type[0] 
         << m_event_info.type[1] 
         << std::endl;

   return m_event_info.r9[0] < cut_value || m_event_info.r9[1] < cut_value ;
}

bool HiggsAnalyzer::pt_cut(float cut_value)
{
   // true if one pt is less than cut_value
   return m_event_info.p[0].pt() < cut_value ||  m_event_info.p[1].pt() < cut_value;
}

const reco::SuperCluster *  HiggsAnalyzer::get_endcap_cluster(const edm::Event& iEvent, const reco::GenParticle & gen)
{
   //edm::Handle<edm::View<reco::PFCluster> > hgcee_clusters;
   //iEvent.getByLabel(edm::InputTag("particleFlowClusterHGCEE"), hgcee_clusters);

   //for( edm::View<reco::PFCluster>::const_iterator cluster =  hgcee_clusters->begin();
   //      cluster != hgcee_clusters->end();
   //      ++cluster)
   //{
   //   if(is_matched(gen,*cluster,0.05))
   //      p4 += math::PtEtaPhiMLorentzVector(cluster->pt(), cluster->eta(), cluster->phi(), 0 );
   //}
   //if( (gen.pt() - p4.pt() )/gen.pt() > .2)
   //{
   //   std::cout << "[" <<  iEvent.id().event() << "]Pt difference high at eta= " << gen.eta() << ": " << gen.pt() << ' ' << p4.pt() <<  std::endl;
   //}

   edm::Handle<reco::SuperClusterCollection> emPFClusters;
   iEvent.getByLabel(edm::InputTag("particleFlowSuperClusterHGCEE"),emPFClusters);

   const reco::SuperCluster * pToCluster = 0;

   for(reco::SuperClusterCollection::const_iterator c_it = emPFClusters->begin();
         c_it!=emPFClusters->end();
         ++c_it)
   {
      if(!is_matched(gen, *c_it, 0.1)) continue;
      if(pToCluster == 0) pToCluster = &(*c_it);
      else if (pToCluster->seed()->energy() < c_it->seed()->energy()) pToCluster = &(*c_it);
   }
   return pToCluster;

}

const reco::Photon * HiggsAnalyzer::get_barrel_photon(const edm::Event& iEvent, const reco::GenParticle & gen)
{
   math::PtEtaPhiMLorentzVector p4(0.0f, 0.0f, 0.0f, 0.0f);
   edm::Handle<edm::View<reco::Photon> > reco_photons;
   iEvent.getByLabel(edm::InputTag("photons"), reco_photons);

   const reco::Photon * pToPhoton = 0;

   for( edm::View<reco::Photon>::const_iterator reco_photon =  reco_photons->begin();
         reco_photon != reco_photons->end();
         ++reco_photon)
   {
      if(is_matched(gen,*reco_photon, 0.1))
      {
         if(pToPhoton == 0) pToPhoton = &(*reco_photon);
         else if(pToPhoton->energy() < reco_photon->energy()) pToPhoton = &(*reco_photon);
      }
   }
   return pToPhoton;
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

   std::cout << "Have " << m_n_unmatched << " unmatched particles\n";
   
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
