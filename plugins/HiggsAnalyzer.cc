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
#include <sstream>

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
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

//#include "DetectorDescription/OfflineDBLoader/interface/GeometryInfoDump.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "UserCode/HGCanalysis/interface/HGCAnalysisTools.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

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
      int m_n_converted;
      float m_mass;
      std::string m_genSource;
      std::string m_jetSource;


      //Geometry info
      unsigned m_n_layers_ee;
      float m_mip_ee;

      std::vector<float> m_layerZ;
      edm::ESHandle<HGCalGeometry> geom;

      float m_linCorr;

      TH1F * m_gen_event_type;
      TH2F * m_interaction_point;
      TH2F * m_interaction_point_cut;

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

      struct part_tree_info_t
      {
         double etamax;
         double phimax;
         double truth_tanangle_x;
         double truth_tanangle_y;
         double truth_E;
         int type;
         std::vector<std::vector<double> > Exy;
      };

      struct tree_info_t
      {
         double truth_vtx_x;
         double truth_vtx_y;
         double truth_vtx_z;
         part_tree_info_t photon[2];

      } m_tree_info;

      TFile * m_tree_file;
      TTree * m_tree;
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


      void fillTree();
      void getMaximumCellFromGeom(const double & phimax,const double & etamax,std::vector<double> & xmax,std::vector<double> & ymax);
      void getMaximumCell(const HGCRecHitCollection &rechitvec,const double & phimax,const double & etamax,std::vector<double> & xmax,std::vector<double> & ymax);
      unsigned getLayer(float z);
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
    m_n_converted = 0;

    m_mip_ee = 55.1;

    edm::Service<TFileService> fs;

    m_gen_event_type         =  fs->make<TH1F>("gen_event_type",         ";Event type;Events",       5,           0,    5);
    m_n_vertices = fs->make<TH1F>("n_vert", ";vertices;Events",180,0,180);
    m_n_vertices_good = fs->make<TH1F>("n_vert_good", ";vertices;Events",180,0,180);

    std::string names[6] = { "All", "!E|B", "B+B", "E+E", "E+B", "E-E"};
    int n_cuts = 5;
    std::string cuts[5] = { "None","converted","reco", "pt", "r9"};

    for ( int i = 1; i < 6; ++i)
       m_gen_event_type->GetXaxis()->SetBinLabel(i,names[i].c_str());

    for ( int i = 0; i < 6; ++i)
       for (int j = 0; j < n_cuts; j++)
       {
          m_gen_photon_pt1 [i][j]  =  fs->make<TH1F>(("gen_photon_pt1"  + names[i] + "_" + cuts[j]).c_str(), ("gen_photon_pt1"  + names[i] + "_" + cuts[j] + ";p_{T1};Photons").c_str()         , 20, 0,   125);
          m_gen_photon_pt2 [i][j]  =  fs->make<TH1F>(("gen_photon_pt2"  + names[i] + "_" + cuts[j]).c_str(), ("gen_photon_pt2"  + names[i] + "_" + cuts[j] + ";p_{T2};Photons").c_str()         , 20, 0,   125);
          m_gen_photon_eta [i][j]  =  fs->make<TH1F>(("gen_photon_eta"  + names[i] + "_" + cuts[j]).c_str(), ("gen_photon_eta"  + names[i] + "_" + cuts[j] + ";#eta;Photons").c_str()           , 20, -4,  4);
          m_gen_photon_dr  [i][j]  =  fs->make<TH1F>(("gen_photon_dr"   + names[i] + "_" + cuts[j]).c_str(), ("gen_photon_dr"   + names[i] + "_" + cuts[j] + ";#Delta R;Photons").c_str()       , 20, 0,   10);
          m_gen_photon_dphi[i][j]  =  fs->make<TH1F>(("gen_photon_dphi" + names[i] + "_" + cuts[j]).c_str(), ("gen_photon_dphi" + names[i] + "_" + cuts[j] + ";#Delta #phi;Photons").c_str()    , 20, 0,   4);
          m_gen_photon_deta[i][j]  =  fs->make<TH1F>(("gen_photon_deta" + names[i] + "_" + cuts[j]).c_str(), ("gen_photon_deta" + names[i] + "_" + cuts[j] + ";#Delta #eta;Photons").c_str()    , 20, 0,   10);

       }

    for ( int i = 0; i < 6; ++i)
      for (int j = 0; j < n_cuts; j++)
         m_diphoton_mass  [i][j]  =  fs->make<TH1F>(("diphoton_mass"   + names[i] + "_" + cuts[j] ).c_str(), ("diphoton_mass_"   + names[i] + "_" + cuts[j] + ";M_{#gamma#gamma};Events/GeV").c_str(), 50, 100, 150);



    for( int j = 0; j < n_cuts; j++)
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

    m_interaction_point = fs->make<TH2F>("int_point","int_point;z;r",400,-400,400,100,0,200);
    m_interaction_point_cut = fs->make<TH2F>("int_point_cut","int_point_cut;z;r",400,-400,400,100,0,200);


    
    m_tree_file = new TFile("HiggsTree.root","recreate");
    m_tree = new TTree("EcellsSR4", "Tree to save energies in each cell of SR4");
    m_tree->SetDirectory(m_tree_file);


    std::string index[2] = {"0","1"}; 
    for( int i = 0; i < 2; i++)
    {
       m_tree->Branch(("etamax_" + index[i]).c_str(), &m_tree_info.photon[i].etamax);
       m_tree->Branch(("phimax_" + index[i]).c_str(), &m_tree_info.photon[i].phimax);
       m_tree->Branch(("truth_tanangle_x_" + index[i]).c_str(), &m_tree_info.photon[i].truth_tanangle_x);
       m_tree->Branch(("truth_tanangle_y_" + index[i]).c_str(), &m_tree_info.photon[i].truth_tanangle_y);
       m_tree->Branch(("truth_E_" + index[i]).c_str(), &m_tree_info.photon[i].truth_E);
       m_tree->Branch(("type_" + index[i]).c_str(), &m_tree_info.photon[i].type);

    }

    m_tree->Branch("event_type",  &m_event_info.eventType);
    m_tree->Branch("truth_vtx_x", &m_tree_info.truth_vtx_x);
    m_tree->Branch("truth_vtx_y", &m_tree_info.truth_vtx_y);
    m_tree->Branch("truth_vtx_z", &m_tree_info.truth_vtx_z);
}


HiggsAnalyzer::~HiggsAnalyzer()
{
 
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

   m_tree_file->cd();
   m_tree->Write();
   m_tree_file->Close();

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HiggsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    //Geometry

    iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",geom);
    const HGCalTopology &topo=geom->topology();
    const HGCalDDDConstants &dddConst=topo.dddConstants();
    m_n_layers_ee = dddConst.layers(true);
    if(m_layerZ.size() == 0)
    {
       for(size_t ilay=0; ilay<m_n_layers_ee; ilay++)
       {
          uint32_t mySubDet(ForwardSubdetector::HGCEE);
          uint32_t recoDetId = (uint32_t)HGCEEDetId(ForwardSubdetector(mySubDet),0,ilay+1,1,0,1);
          //require to be on the same side of the generated particle
          const GlobalPoint pos( std::move( geom->getPosition(recoDetId) ) );
          m_layerZ.push_back(fabs(pos.z()));
       }

       // Setup Size of Exy and TTree
       std::vector<double> init(25,0);
       for(unsigned i = 0; i < 2; ++i)
       {
          m_tree_info.photon[i].Exy.resize(m_n_layers_ee,init);
          for (unsigned iL(0); iL < m_n_layers_ee; ++iL)
          {
             std::ostringstream label;
             for (unsigned iy(0);iy<5;++iy)
             {
                for (unsigned ix(0);ix<5;++ix)
                {
                   unsigned idx = 5*iy+ix;
                   label.str("");
                   label << "E" << i << "_" << iL << "_" << idx;
                   m_tree->Branch(label.str().c_str(), &m_tree_info.photon[i].Exy[iL][idx]);
                }
             }
          }
       }
    }

    //Reinit Exy
    if(m_n_layers_ee)
       for(unsigned i = 0; i < 2; ++i)
          for (unsigned iL(0); iL < m_n_layers_ee; ++iL)
             for (unsigned idx(0);idx<25;++idx)
                m_tree_info.photon[i].Exy[iL][idx] = 0.0;

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
       photon0 = NULL;
       photon1 = NULL;
       //if(genParticle->pdgId() == 25 && iEvent.id().event() == 47)
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
           //std::cout << "Found Higgs with 2(3?) daughters\n";
           photon0 = dynamic_cast<const reco::GenParticle *>(genParticle->daughter(0));
           photon1 = dynamic_cast<const reco::GenParticle *>(genParticle->daughter(1));
           m_event_info.vertex = genParticle->vertex();

           if( photon0->pdgId() != 22 || photon0->pdgId() != 22 || photon0->status() != 1 || photon1->status() != 1) // Not Higgs to gg
           {
               std::cout << "daughters are not stable photons\n";
               photon0 = NULL;
               photon1 = NULL;
               continue;
           }



           //std::cout << "Found photons\n";
           break; //Found photons, leave the loop
       }
   }
   //Quit if we didn't find photons
   if(!photon0 || !photon1) return;

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

   //Lets get rid of converted photons first

   //Geant4 collections
   edm::Handle<std::vector<SimTrack> > SimTk;
   iEvent.getByLabel("g4SimHits",SimTk);
   edm::Handle<std::vector<SimVertex> > SimVtx;
   iEvent.getByLabel("g4SimHits",SimVtx); 
   edm::Handle<std::vector<int> > barcodes;
   iEvent.getByLabel("genParticles",barcodes); 


   //get trackIds for photons
   int trackId0 = 0;
   int trackId1 = 0;

   //First find genPartIndex
   int genPartIndex0 = 0;
   int genPartIndex1 = 0;
   int i =0;
   for(const reco::GenParticle &genParticle : *genParticles)
   //for( edm::View<reco::GenParticle>::const_iterator genParticle =  genParticles->begin();
           //genParticle != genParticles->end();
           //++genParticle)
   {
      i++;
      if( &genParticle == photon0) genPartIndex0 = i;
      if( &genParticle == photon1) genPartIndex1 = i;
      if (genPartIndex0 && genPartIndex1) break;
   }
   //loop over tracks to find track id for correct genPartIndex
   for (const SimTrack &simTrk : *SimTk) 
   {
      if(simTrk.genpartIndex() == genPartIndex0) trackId0 = simTrk.trackId();
      if(simTrk.genpartIndex() == genPartIndex1) trackId1 = simTrk.trackId();
      if (trackId0 && trackId1) break;
   }

   G4InteractionPositionInfo intInfo=getInteractionPosition(SimTk.product(),SimVtx.product(),trackId0);
   math::XYZVectorD hitPos0=intInfo.pos;
   m_interaction_point->Fill(hitPos0.z(),hitPos0.rho());

   intInfo=getInteractionPosition(SimTk.product(),SimVtx.product(),trackId1);
   math::XYZVectorD hitPos1=intInfo.pos;
   m_interaction_point->Fill(hitPos1.z(),hitPos1.rho());

   int converted = 0;
   //Check if this is a converted photon
   if(fabs(hitPos0.z()) < 295.0f && hitPos0.rho() < 112.0f)
   {
      std::cout << "Event " << iEvent.id().event() << " has converted photon(pt=" << photon0->pt() << " GeV)\n";
      converted++;
   }
   if(fabs(hitPos1.z()) < 295.0f && hitPos1.rho() < 112.0f)
   {
      std::cout << "Event " << iEvent.id().event() << " has converted photon(pt=" << photon1->pt() << " GeV)\n";
      converted++;
   }
   if(converted)
   {
      m_n_converted += converted;
      return;
   }
   m_interaction_point_cut->Fill(hitPos0.z(),hitPos0.rho());
   m_interaction_point_cut->Fill(hitPos1.z(),hitPos1.rho());

   fillCutPlots(cuti++);
   m_mass = get_mass();
   if(!m_mass) return;

   fillTree();

   fillCutPlots(cuti++);

   if( pt_cut(20.0f)) return;

   fillCutPlots(cuti++);

   if( r9_cut(.85)) return;

   fillCutPlots(cuti++);



}

void HiggsAnalyzer::getMaximumCellFromGeom(const double & phimax,const double
      & etamax,std::vector<double> & xmax,std::vector<double> & ymax){

   for (unsigned iL(0); iL<m_n_layers_ee;++iL){
      double theta = 2*atan(exp(-1.*etamax));
      double rho = m_layerZ[iL]/cos(theta);
      xmax[iL] = rho*sin(theta)*cos(phimax);
      ymax[iL] = rho*sin(theta)*sin(phimax);

      if (xmax[iL]>0) xmax[iL]=static_cast<int>((xmax[iL]+4.999999)/10.)*10;
      else xmax[iL]=static_cast<int>((xmax[iL]-4.999999)/10.)*10;
      if (ymax[iL]>0) ymax[iL]=static_cast<int>((ymax[iL]+4.999999)/10.)*10;
      else ymax[iL]=static_cast<int>((ymax[iL]-4.999999)/10.)*10;

   }//loop on layers

}

unsigned HiggsAnalyzer::getLayer(float z)
{

   //match nearest HGC layer in Z
   float dzMin(99999999.);
   int layer = 0;
   for(size_t ilay=0; ilay<m_layerZ.size(); ilay++)
   {
      float dz=fabs(fabs(m_layerZ[ilay])-fabs(z));
      if(dz>dzMin) continue;
      dzMin=dz;
      layer = ilay;
   }
   return layer;
}

void HiggsAnalyzer::getMaximumCell(const HGCRecHitCollection &rechitvec,const
 double & phimax,const double & etamax,std::vector<double> &
 xmax,std::vector<double> & ymax){

   std::vector<double> xmaxgeom;
   xmaxgeom.resize(m_n_layers_ee,0);
   std::vector<double> ymaxgeom;
   ymaxgeom.resize(m_n_layers_ee,0);
   getMaximumCellFromGeom(phimax,etamax,xmaxgeom,ymaxgeom);

   //choose cell with maximum energy from 3*3 array around geom pos of max.
   std::vector<double> Emax;
   Emax.resize(m_n_layers_ee,0);

   for (unsigned iH(0); iH<(rechitvec).size(); ++iH){//loop on rechits
     const HGCRecHit & lHit = (rechitvec)[iH];
     const GlobalPoint pos = geom->getPosition(lHit.id());
     if(etamax*pos.eta() < 0 ) continue;

     unsigned layer = getLayer(pos.z());
     double posx = pos.x();
     double posy = pos.y();
     double energy = lHit.energy();

 //15 assumes unit is mm ! If cm in CMSSW -> change for 1.5 ...

     if (fabs(posx-xmaxgeom[layer]) <= 1.5 &&
         fabs(posy-ymaxgeom[layer]) <= 1.5){

       if (energy>Emax[layer]){
         Emax[layer] = energy;
         xmax[layer] = posx;
         ymax[layer] = posy;
       }
     }

   }//loop on rechits

 }



void HiggsAnalyzer::fillTree()
{

   //std::cout << "FillTree\n";
   for( int i = 0; i < 2; ++i)
   {
      m_tree_info.photon[i].etamax = m_event_info.p[i].eta();
      m_tree_info.photon[i].phimax = m_event_info.p[i].phi();
      m_tree_info.photon[i].truth_tanangle_x =  m_event_info.photon[i]->px() / m_event_info.photon[i]->pz();
      m_tree_info.photon[i].truth_tanangle_y =  m_event_info.photon[i]->py() / m_event_info.photon[i]->pz();
      m_tree_info.photon[i].truth_E = m_event_info.photon[i]->energy();
      m_tree_info.photon[i].type = (int) m_event_info.type[i];

      //std::vector<std::vector<double> > Exy;

      if(m_event_info.type[i] == ENDCAP_PHOTON)
      {

         //std::cout << "Grab RecHits\n";
         edm::Handle<HGCRecHitCollection> rechitvec;
         m_event_info.iEvent->getByLabel(edm::InputTag("HGCalRecHit","HGCEERecHits"),rechitvec);


         std::vector<double> xmax;
         xmax.resize(m_n_layers_ee,0);
         std::vector<double> ymax;
         ymax.resize(m_n_layers_ee,0);

         //std::cout << "getMaximumCell\n";
         getMaximumCell(*rechitvec,m_tree_info.photon[i].phimax,m_tree_info.photon[i].etamax,xmax,ymax);

         double steplarge = 25+0.1;//+0.1 to accomodate double precision

         //std::cout << "Loop Over RecHits\n";
         for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on rechits
            const HGCRecHit & lHit = (*rechitvec)[iH];
            const GlobalPoint pos = geom->getPosition(lHit.id());

            if (pos.eta() * m_tree_info.photon[i].etamax < 0) continue;

            double energy = lHit.energy()/m_mip_ee;//in MIP already...
            unsigned layer = getLayer(pos.z());
            double posx = pos.x();
            double posy = pos.y();

            if (fabs(posx-xmax[layer]) < steplarge &&
                  fabs(posy-ymax[layer]) < steplarge){
               //double leta = lHit.eta();
               double lCor =  0;
               //puDensity_.getDensity(leta,layer,geomConv_.cellSizeInCm(layer,leta),nPU);
               energy = std::max(0.,energy - lCor);

               int ix = (posx-xmax[layer])/10.;
               int iy = (posy-ymax[layer])/10.;
               unsigned idx = 0;
               if ((ix > 2 || ix < -2) || (iy>2 || iy<-2)) {
                  std::cout << " error, check ix=" << ix << " iy=" << iy  << " "
                            << "posx,y-max=" << posx-xmax[layer] << " " 
                            << posy-ymax[layer] << " steplarge " <<
                            steplarge << std::endl;
                  continue;
               }
               else
                  idx = 5*(iy+2)+(ix+2);
               //std::cout << "Store energy at layer=" << layer << " idx=" << idx << std::endl;

               m_tree_info.photon[i].Exy.at(layer).at(idx) = energy;
            }

         }//loop on rechits
      }
   }

   m_tree_info.truth_vtx_x = m_event_info.vertex.x();
   m_tree_info.truth_vtx_y = m_event_info.vertex.y();
   m_tree_info.truth_vtx_z = m_event_info.vertex.z();

   m_tree->Fill();

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
   std::cout << "Have " << m_n_converted << " converted particles\n";
   
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
