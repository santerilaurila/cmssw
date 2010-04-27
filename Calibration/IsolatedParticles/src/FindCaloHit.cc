#include "Calibration/IsolatedParticles/interface/FindCaloHit.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include <iostream>

namespace spr {

  std::vector<EcalRecHitCollection::const_iterator> find(edm::Handle<EcalRecHitCollection>& hits, DetId thisDet, bool debug) {

    std::vector<EcalRecHitCollection::const_iterator> hit;
    hit.push_back( hits->find(thisDet) );
    return hit;
  }

  std::vector<HBHERecHitCollection::const_iterator> find(edm::Handle<HBHERecHitCollection>& hits, DetId thisDet, bool debug) {
    std::vector<HBHERecHitCollection::const_iterator> hit;
    hit.push_back( hits->find(thisDet) );
    return hit;
  }

  std::vector<edm::PCaloHitContainer::const_iterator> find(edm::Handle<edm::PCaloHitContainer>& hits, DetId thisDet, bool debug) {

    std::vector<edm::PCaloHitContainer::const_iterator> hit;

    edm::PCaloHitContainer::const_iterator ihit;
    for (ihit=hits->begin(); ihit!=hits->end(); ihit++) {
      DetId detId(ihit->id());
      if (detId == thisDet) {
        hit.push_back(ihit);
      }
    }
  
    return hit;
  }
}
