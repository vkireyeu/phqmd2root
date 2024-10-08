#include "PHQMDevent.h"

namespace PHQMD{

/** Standard constructor
*/
  Particle::Particle(){
    fID      = 0;
    fPDG     = 0;
    fCharge  = 0;
    fChannel = 0;
    fBarDens = 0;
    fEnDens  = 0;
    fIsInMst = 0;
  }

  Particle::Particle(const Particle& in){
    *this = in;
  }

  Particle::~Particle(){
  }



  Baryon::Baryon(){
    fPHSD_POS = 0;
    fPDG      = 0;
    fPHSD_ID  = 0;
    fFID      = 0;
    fFSIZE    = 0;
    fPREG     = 0;
    fPCH      = 0;
    fEBND     = 0;
  }

  Baryon::Baryon(const Baryon& in){
    *this = in;
  }

  Baryon::~Baryon(){
  }



  Fragment::Fragment(){
    fFID      = 0;
    fPDG      = 0;
  }

  Fragment::Fragment(const Fragment& in){
    *this = in;
  }

  Fragment::~Fragment(){
  }

}
