#include "UserCode/HiggsAnalysis/interface/MCBarrelEndcapFilter.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <iostream>

using namespace edm;
using namespace std;


MCBarrelEndcapFilter::MCBarrelEndcapFilter(const edm::ParameterSet& iConfig) :
    label_(iConfig.getUntrackedParameter("moduleLabel",std::string("generator"))),
       barrelID(iConfig.getUntrackedParameter("barrelID", 22)),
       endcapID(iConfig.getUntrackedParameter("endcapID", 22)),
       ptMinBarrel(iConfig.getUntrackedParameter("ptMinBarrel", 10.0)),
       ptMinEndcap(iConfig.getUntrackedParameter("ptMinEndcap", 10.0)),
       etaMinBarrel(iConfig.getUntrackedParameter("etaMinBarrel", 0.0)),
       etaMaxBarrel(iConfig.getUntrackedParameter("etaMaxBarrel", 1.5)),
       etaMinEndcap(iConfig.getUntrackedParameter("etaMinEndcap", 1.6)),
       etaMaxEndcap(iConfig.getUntrackedParameter("etaMaxEndcap", 2.9)),
       statusBarrel(iConfig.getUntrackedParameter("statusBarrel", 1)),
       statusEndcap(iConfig.getUntrackedParameter("statusEndcap", 1)),
       minInvMass(iConfig.getUntrackedParameter("minInvMass", 122.0)),
       maxInvMass(iConfig.getUntrackedParameter("maxInvMass", 128.0)),
       verbose(iConfig.getUntrackedParameter("verbose", false))
{
    //here do whatever other initialization is needed
}


MCBarrelEndcapFilter::~MCBarrelEndcapFilter()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}


// ------------ method called to skim the data  ------------
bool MCBarrelEndcapFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    bool accepted = false;
    Handle<HepMCProduct> evt;
    iEvent.getByLabel(label_, evt);

    vector<HepMC::GenParticle*> barrel_passed;
    vector<HepMC::GenParticle*> endcap_passed;


    const HepMC::GenEvent * myGenEvent = evt->GetEvent();

    if(verbose) cout << "Start Event\n";

    for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin();
            p != myGenEvent->particles_end(); ++p ) {


        if(verbose) cout << abs((*p)->pdg_id()) << ' ' 
            << (*p)->momentum().perp() << ' ' 
            << (*p)->momentum().eta() << ' ' 
            << (*p)->status() << endl;
        // check for barrel conditions
        if ((abs((*p)->pdg_id()) == abs(barrelID) || barrelID == 0) &&
                (*p)->momentum().perp() > ptMinBarrel && (*p)->momentum().eta() > etaMinBarrel 
                && (*p)->momentum().eta() < etaMaxBarrel && ((*p)->status() == statusBarrel || statusBarrel == 0))
        { 
            if(verbose) cout << "Found Barrel\n";
            // passed barrel conditions ...
            // ... now check pair-conditions with endcap passed particles
            //
            unsigned int i=0;
            double invmass =0.;
            while(!accepted && i<endcap_passed.size()) {
                invmass = ( math::PtEtaPhiMLorentzVector((*p)->momentum()) + math::PtEtaPhiMLorentzVector(endcap_passed[i]->momentum())).mass();
                if(verbose) cout << "Invariant Mass is" << invmass << endl;
                if(invmass > minInvMass && invmass < maxInvMass) {
                    accepted = true;
                }
                i++;
            }    
            // if we found a matching pair quit the loop
            if(accepted) break;
            else{
                barrel_passed.push_back(*p);   // else remember the particle to have passed barrel conditions
            }
        }

        // check for endcap conditions

        if ((abs((*p)->pdg_id()) == abs(endcapID) || endcapID == 0) && 
                (*p)->momentum().perp() > ptMinEndcap && (*p)->momentum().eta() > etaMinEndcap
                && (*p)->momentum().eta() < etaMaxEndcap && ((*p)->status() == statusEndcap || statusEndcap == 0)) { 
            if(verbose) cout << "Found Endcap\n";
            // passed endcap conditions ...
            // ... now check pair-conditions with barrel passed particles vector
            unsigned int i=0;
            double invmass =0.;
            while(!accepted && i<barrel_passed.size()) {
                if((*p) != barrel_passed[i]) {
                    if(verbose) cout << "Checking to match Barrel\n";
                    invmass = ( math::PtEtaPhiMLorentzVector((*p)->momentum()) + math::PtEtaPhiMLorentzVector(barrel_passed[i]->momentum())).mass();
                    if(verbose) cout << "Invariant Mass is" << invmass << endl;
                    if(invmass > minInvMass && invmass < maxInvMass) {
                        accepted = true;
                    }
                    i++;
                }
            }    
            // if we found a matching pair quit the loop
            if(accepted) break;
            else {
                endcap_passed.push_back(*p);   // else remember the particle to have passed endcap conditions
            }
        }
    }

    if(verbose)
    {
        cout << "This event ";
        if (accepted){ cout <<  "passed "; } else {cout << "failed ";}
        cout << "the filter\n";
    }

    return accepted;

}

//define this as a plug-in
DEFINE_FWK_MODULE(MCBarrelEndcapFilter);
