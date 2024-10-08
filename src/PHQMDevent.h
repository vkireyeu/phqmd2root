#ifndef PHQMDEVENT
#define PHQMDEVENT

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TBranch.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TObject.h>
#include <TNamed.h>
#include <TMath.h>
#include <Math/Vector4D.h>

#include <array>
#include <vector>
#include <map>

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <math.h>

namespace PHQMD{

/**  \class PHQMD::Particle
  PHQMD::Particle is a class to store the PHSD/PHQMD particle information ('phsd.dat' file analogue)
  \author Viktar Kireyeu
*/
  class Particle{
    private:
      int                      fID; ///< Position in the PHSD vector
      int                     fPDG; ///< Monte Carlo particle code, [see the scheme](https://pdg.lbl.gov/2022/reviews/rpp2022-rev-monte-carlo-numbering.pdf)
      int                  fCharge; ///< Charge
      int                 fChannel; ///< Production channel
      float               fBarDens; ///< Baryonic density at the freeze-out time at the particle position
      float                fEnDens; ///< Energy density at the freeze-out time at the particle position
      bool                fIsInMst; ///< Flag to store the information about the particle's entrance into the MST clusterization
      ROOT::Math::XYZTVector    fR; ///< Freeze-out coordinates [4-vector](https://root.cern.ch/doc/master/group__GenVector.html)
      ROOT::Math::PxPyPzEVector fP; ///< Momentum [4-vector](https://root.cern.ch/doc/master/group__GenVector.html)

    public:
      Particle(); ///< Default constructor
      virtual ~Particle(); ///< Default destructor
      Particle(const Particle& in); ///< Copy constructor

      inline int                      GetId () const {return      fID;}
      inline int                     GetPdg () const {return     fPDG;}
      inline int                  GetCharge () const {return  fCharge;}
      inline int                 GetChannel () const {return fChannel;}
      inline float               GetBarDens () const {return fBarDens;}
      inline float                GetEnDens () const {return  fEnDens;}
      inline bool                   IsInMst () const {return fIsInMst;}
      inline ROOT::Math::XYZTVector    GetR () const {return       fR;}
      inline ROOT::Math::PxPyPzEVector GetP () const {return       fP;}

      inline void SetId      (int                       val) {fID      = val;}
      inline void SetPdg     (int                       val) {fPDG     = val;}
      inline void SetCharge  (int                       val) {fCharge  = val;}
      inline void SetChannel (int                       val) {fChannel = val;}
      inline void SetBarDens (float                     val) {fBarDens = val;}
      inline void SetEnDens  (float                     val) {fEnDens  = val;}
      inline void SetInMst   (bool                      val) {fIsInMst = val;}
      inline void SetR       (ROOT::Math::XYZTVector    val) {fR       = val;}
      inline void SetP       (ROOT::Math::PxPyPzEVector val) {fP       = val;}
  };



/**  \class PHQMD::Baryon
  PHQMD::Baryon is a class to store the MST baryons information ('fort.791' of 'fort.891' file analogue)
  \author Viktar Kireyeu
*/
  class Baryon {
    private:
      int                fPHSD_POS; ///< Position of the corrsponding particle in the Event::particles vector
      int                     fPDG; ///< Monte Carlo particle code, [see the scheme](https://pdg.lbl.gov/2022/reviews/rpp2022-rev-monte-carlo-numbering.pdf)
      int                 fPHSD_ID; ///< Position in the PHSD vector
      int                     fFID; ///< ID of the cluster, to which this baryon is assigned
      int                   fFSIZE; ///< Size of the cluster, to which this baryon is assigned
      int                    fPREG; ///< Production region of the baryon
      int                     fPCH; ///< Production channel
      float                  fEBND; ///< Binding energy per nuceon of the cluster, to which this baryon is assigned
      ROOT::Math::XYZTVector    fR; ///< Coordinates [4-vector](https://root.cern.ch/doc/master/group__GenVector.html)
      ROOT::Math::PxPyPzMVector fP; ///< Momentum [4-vector](https://root.cern.ch/doc/master/group__GenVector.html)

    public:
      Baryon();
      virtual ~Baryon();
      Baryon(const Baryon& in);

      inline int                 GetPhsdPos () const {return fPHSD_POS;}
      inline int                     GetPdg () const {return      fPDG;}
      inline int                  GetPhsdId () const {return  fPHSD_ID;}
      inline int                     GetFid () const {return      fFID;}
      inline int                    GetSize () const {return    fFSIZE;}
      inline int              GetProdRegion () const {return     fPREG;}
      inline int             GetProdChannel () const {return      fPCH;}
      inline float                 GetEbind () const {return     fEBND;}
      inline bool                   IsBound () const {return fEBND < -1E-5 ? 1 : 0;}
      inline ROOT::Math::XYZTVector    GetR () const {return        fR;}
      inline ROOT::Math::PxPyPzMVector GetP () const {return        fP;}

      inline void SetPhsdPos     (int                    val) {fPHSD_POS = val;}
      inline void SetPdg         (int                    val) {fPDG      = val;}
      inline void SetPhsdId      (int                    val) {fPHSD_ID  = val;}
      inline void SetFid         (int                    val) {fFID      = val;}
      inline void SetFsize       (int                    val) {fFSIZE    = val;}
      inline void SetProdRegion  (int                    val) {fPREG     = val;}
      inline void SetProdChannel (int                    val) {fPCH      = val;}
      inline void SetEbind       (float                  val) {fEBND     = val;}
      inline void SetR        (ROOT::Math::XYZTVector    val) {fR        = val;}
      inline void SetP        (ROOT::Math::PxPyPzMVector val) {fP        = val;}
  };



/**  \class PHQMD::Fragment
  PHQMD::Fragment is a class to store the MST clusters information ('fort.790' file analogue)
  \author Viktar Kireyeu
*/
  class Fragment {
    private:
      int                     fFID; ///< ID of the cluster
      int                     fPDG; ///< Monte Carlo particle code, [see the scheme](https://pdg.lbl.gov/2022/reviews/rpp2022-rev-monte-carlo-numbering.pdf)
      ROOT::Math::XYZTVector    fR; ///< Coordinates [4-vector](https://root.cern.ch/doc/master/group__GenVector.html)
      ROOT::Math::PxPyPzMVector fP; ///< Momentum [4-vector](https://root.cern.ch/doc/master/group__GenVector.html)

    public:
      Fragment();
      virtual ~Fragment();
      Fragment(const Fragment& in);

      inline int                     GetPdg () const {return      fPDG;}
      inline int                     GetFid () const {return      fFID;}
      inline ROOT::Math::XYZTVector    GetR () const {return        fR;}
      inline ROOT::Math::PxPyPzMVector GetP () const {return        fP;}
      inline int                       GetA () const {return  (fPDG %   10000) /   10;}
      inline int                       GetZ () const {return  (fPDG /   10000) % 1000;}
      inline int                       GetL () const {return  (fPDG /10000000) %   10;}

      inline void SetPdg         (int                    val) {fPDG      = val;}
      inline void SetFid         (int                    val) {fFID      = val;}
      inline void SetR        (ROOT::Math::XYZTVector    val) {fR        = val;}
      inline void SetP        (ROOT::Math::PxPyPzMVector val) {fP        = val;}
  };



/**  \class PHQMD::Event
  PHQMD::Event is a class to store the PHQMD::Particle, PHQMD::Baryon and PHQMD::Fragment vectors
  together with the event-wise information.
  \author Viktar Kireyeu
*/
  class Event {
    private:
      float                                     fB; ///< Impact parameter
      int                              fNparticles; ///< Number of particles in the event
      int                           fNparticipants; ///< Number of participants in the event
      std::array<float,4>                     fPsi; ///< Participant plane angle of the n-th harmonic: 2,3,4,5
      std::array<float,4>                     fEcc; ///< Eccentricity of the n-th harmonic: 2,3,4,5
      std::vector<Particle>              particles; ///< Vector to store PHQMD::Particles
      std::vector<std::vector<Baryon>>   mstbsteps; ///< Vector to store vectors (for the each MST time step) of PHQMD::Baryon
      std::vector<std::vector<Fragment>> mstfsteps; ///< Vector to store vectors (for the each MST time step) of PHQMD::Fragment

    public:
      inline double                  GetB       () const {return             fB;}
      inline int            GetNparticles       () const {return    fNparticles;}
      inline int         GetNparticipants       () const {return fNparticipants;}
      inline std::array<float,4>   GetPsi       () const {return           fPsi;}
      inline std::array<float,4>   GetEcc       () const {return           fEcc;}
      inline std::vector<Particle> *GetParticles ()       {return      &particles;}
      inline std::vector<std::vector<Baryon>>   GetMstBSteps () const {return mstbsteps;}
      inline std::vector<std::vector<Fragment>> GetMstFSteps () const {return mstfsteps;}

      inline void SetB	           (double               val) {fB                  = val;}
      inline void SetNparticles    (int                  val) {fNparticles         = val;}
      inline void SetNparticipants (int                  val) {fNparticipants      = val;}
      inline void SetPsi           (std::array<float,4>  val) {fPsi                = val;}
      inline void SetEcc           (std::array<float,4>  val) {fEcc                = val;}
      inline void AddParticle      (Particle            part) {particles.push_back(part);}
      inline void AddMstBStep      (std::vector<Baryon>   step) { mstbsteps.push_back(step);}
      inline void AddMstFStep      (std::vector<Fragment> step) { mstfsteps.push_back(step);}
  };



/**  \class PHQMD::Run
  PHQMD::Run is a class to store the PHQMD run settings ('inputPHSD' file analogue)
  \author Viktar Kireyeu
*/
  class Run {
    private:
      int    fNevents; /// Number of events in the run
      int       fAtar; /// Target mass
      int       fZtar; /// Protons in target
      int      fAproj; /// Projectile mass
      int      fZproj; /// Protons in projectile
      float     fTkin; /// Kinetic energy (GeV)
      float     fBmin; /// BMIN:   minimal impact parameter in fm ! no effect for p+A
      float     fBmax; /// BMAX:   maximal impact parameter in fm ! no effect for p+A
      float   fDeltaB; /// DeltaB: impact parameter step in fm (used only if IBweight_MC=0)
      int      fIGLUE; /// IGLUE: =1 with partonic QGP phase (PHSD mode); =0 - HSD mode
      int        fIBW; /// IBweight_MC: =0 constant step in B =DBIMP; =1 choose B by Monte-Carlo ISUBS times in [Bmin,Bmax]
      int        fNUM; /// NUM:    optimized number of parallel ensambles ("events")
      int        fSUB; /// ISUBS:  number of subsequent runs
      int       fSEED; /// Random seed
      float   fFinalT; /// FINALT: final time of calculation in fm/c
      int       fILOW; /// Output verbosity level
      int       fIdil; /// Idilept: =0 no dileptons; =1 electron pair; =2  muon pair
      int        fICQ; /// ICQ: =0 free rho's, =1 dropping mass, =2 broadening, =3 drop.+broad.
      int      fIHARD; /// IHARD: =1 with charm and bottom; =0 - without
      int      fIDQPM; /// IDQPM: ???
      int      fIUSER; ///  =1 for general users : use default /optimized settings; = 0 for PHSD team
      int    fINUCLEI; ///  INUCLEI  =1 reactions with deuterons (starting from 4d33166a)
      int     fIPHQMD; /// IPHQMD=1: propagation with QMD dynamics; =0 with HSD/PHSD dynamics
      int      fISACA; /// ISACA: enable or disable SACA output (note: SACA or MST is controled by iflagsaca)
      float    fTSACA; /// TSACA: starting time for SACA
      float   fDTSACA; /// DTSACA, time step for SACA calculations
      int     fNTSACA; /// NTSACA, Number of SACA timesteps ( =>  tmax=tsaca+dtsaca*ntsaca must be < FINALT !)
      int    fFLGSACA; /// iflagsaca: =1 SACA(=FRIGA) analysis(+MST), =0 MST, no SACA --!!!! FRIGA
      int        fYuk; /// tageyuk=1  ! if Yukawa potential requested in SACA       
      int        fAsy; /// tageasy=1  ! if asymmetry energy requested in SACA 
      int       fPair; /// tagepair=1 ! tructure effects (pairing,...) requested in SACA 
      int       fCoul; /// tagecoul=0 ! activates the coulomb energy for fragment selection
      float     fAsy0; /// vasy0=23.3 ! asymmetry potential energy at normal density (MeV) in SACA 
      float    fEPair; /// eta_pairing=0.0 ! pairing potential exponant (0.->only forbids unbound isotopes, 1., 0.65, 0.35 or 0.25) in SACA 
      int        fEOS; /// iqmdeos  ! EoS for QMD option IPHQMD=1 ; =0: hard EOS without M.D.I; =1 soft EoS
      int    fResSACA; /// IFLAG_Res_SACA ! =1 include ALL resonances with their decay to SACA; =2 - only nucleons; =3 nucleons and hyperons
      int    fWigDens; /// IfragWigDen  ! =0: no; =1 yes, light clusters formation according to the Wigner density


    public:

      inline int GetNevents () const {return fNevents;}
      inline Int_t     GetAtar() const {return  fAtar;}
      inline Int_t     GetZtar() const {return  fZtar;}
      inline Int_t    GetAproj() const {return fAproj;}
      inline Int_t    GetZproj() const {return fZproj;}

      inline float   GetTkin() const {return    fTkin;}
      inline float   GetBmin() const {return    fBmin;}
      inline float   GetBmax() const {return    fBmax;}
      inline float GetDeltaB() const {return  fDeltaB;}
      inline int    GetIGlue() const {return   fIGLUE;}
      inline int      GetIBW() const {return     fIBW;}
      inline int      GetNum() const {return     fNUM;}
      inline int      GetSub() const {return     fSUB;}
      inline float GetFinalT() const {return  fFinalT;}
      inline int     GetIdil() const {return    fIdil;}
      inline int      GetICQ() const {return     fICQ;}
      inline int    GetIHARD() const {return   fIHARD;}
      inline int    GetIDQPM() const {return   fIDQPM;}
      inline int  GetINUCLEI() const {return fINUCLEI;}
      inline int   GetIPHQMD() const {return  fIPHQMD;}
      inline int    GetISACA() const {return   fISACA;}
      inline float  GetTSACA() const {return   fTSACA;}
      inline float GetDTSACA() const {return  fDTSACA;}
      inline int   GetNTSACA() const {return  fNTSACA;}
      inline int  GetFLGSACA() const {return fFLGSACA;}
      inline int      GetYuk() const {return     fYuk;}
      inline int      GetAsy() const {return     fAsy;}
      inline int     GetPair() const {return    fPair;}
      inline int     GetCoul() const {return    fCoul;}
      inline float   GetAsy0() const {return    fAsy0;}
      inline float  GetEPair() const {return   fEPair;}
      inline int      GetEOS() const {return     fEOS;}
      inline int  GetResSACA() const {return fResSACA;}
      inline int  GetWigDens() const {return fWigDens;}
      inline int     GetSeed() const {return    fSEED;}
      inline int     GetIlow() const {return    fILOW;}
      inline int    GetIuser() const {return   fIUSER;}

      inline float   GetPlab() const {return TMath::Sqrt(fTkin*fTkin + 2.*0.938*fTkin);}
      inline float    GetSRT() const {return TMath::Sqrt(4.*(0.938*0.938) + 2.*0.938*fTkin);}

      inline void SetNevents (int   val) {fNevents = val;}
      inline void SetAtar    (int   val) {fAtar    = val;}
      inline void SetZtar    (int   val) {fZtar    = val;}
      inline void SetAproj   (int   val) {fAproj   = val;}
      inline void SetZproj   (int   val) {fZproj   = val;}
      inline void SetTkin    (float val) {fTkin    = val;}
      inline void SetBmin    (float val) {fBmin    = val;}
      inline void SetBmax    (float val) {fBmax    = val;}
      inline void SetDeltaB  (float val) {fDeltaB  = val;}
      inline void SetIGLUE   (int   val) {fIGLUE   = val;}
      inline void SetIBW     (int   val) {fIBW     = val;}
      inline void SetNum     (int   val) {fNUM     = val;}
      inline void SetSub     (int   val) {fSUB     = val;}
      inline void SetFinalT  (int   val) {fFinalT  = val;}
      inline void SetIdil    (int   val) {fIdil    = val;}
      inline void SetICQ     (int   val) {fICQ     = val;}
      inline void SetIHARD   (int   val) {fIHARD   = val;}
      inline void SetIDQPM   (int   val) {fIDQPM   = val;}
      inline void SetINUCLEI (int   val) {fINUCLEI = val;}
      inline void SetIPHQMD  (int   val) {fIPHQMD  = val;}
      inline void SetISACA   (int   val) {fISACA   = val;}
      inline void SetTSACA   (float val) {fTSACA   = val;}
      inline void SetDTSACA  (float val) {fDTSACA  = val;}
      inline void SetNTSACA  (int   val) {fNTSACA  = val;}
      inline void SetFLGSACA (int   val) {fFLGSACA = val;}
      inline void SetYuk     (int   val) {fYuk     = val;}
      inline void SetAsy     (int   val) {fAsy     = val;}
      inline void SetPair    (int   val) {fPair    = val;}
      inline void SetCoul    (int   val) {fCoul    = val;}
      inline void SetAsy0    (float val) {fAsy0    = val;}
      inline void SetEPair   (float val) {fEPair   = val;}
      inline void SetEOS     (int   val) {fEOS     = val;}
      inline void SetResSACA (int   val) {fResSACA = val;}
      inline void SetWigDens (int   val) {fWigDens = val;}
      inline void SetSeed    (int   val) {fSEED    = val;}
      inline void SetIlow    (int   val) {fILOW    = val;}
      inline void SetIuser   (int   val) {fIUSER   = val;}
  };

}
#endif
