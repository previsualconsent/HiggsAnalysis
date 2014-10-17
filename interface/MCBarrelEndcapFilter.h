#ifndef MCParticlePairFilter_h
#define MCParticlePairFilter_h
// -*- C++ -*-
//
// Package:    MCBarrelEndcapFilter
// Class:      MCBarrelEndcapFilter
// 
/* 

 Description: filter events based on the Pythia particle information

 Implementation: inherits from generic EDFilter
     
*/
//
// Original Author:  Fabian Stoeckli
//         Created:  Mon Sept 11 10:57:54 CET 2006
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


//
// class decleration
//

class MCBarrelEndcapFilter : public edm::EDFilter {
   public:
      explicit MCBarrelEndcapFilter(const edm::ParameterSet&);
      ~MCBarrelEndcapFilter();


      virtual bool filter(edm::Event&, const edm::EventSetup&);
   private:
      // ----------memeber function----------------------

      // ----------member data ---------------------------
      
       std::string label_;
       int barrelID;
       int endcapID;
       double ptMinBarrel;
       double ptMinEndcap;
       double etaMinBarrel;  
       double etaMaxBarrel;
       double etaMinEndcap;  
       double etaMaxEndcap;
       int statusBarrel;
       int statusEndcap;
       double minInvMass;
       double maxInvMass;
       
};
#endif

