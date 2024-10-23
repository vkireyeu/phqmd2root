#include "PHQMDevent.h"

using namespace PHQMD;

int main(int argc, char* argv[]){


  char       fbuffer[512];
  const char *filename;
  const char *f791name;
  const char *inputPHSD;
  gzFile     fgzFile;
  gzFile     f791File;

  int irun = 0;
  if(argc < 5){
    printf("usage: %s <inputPHSD> <phsd.dat.gz> <fort.791> <out.root>\n", argv[0]);
    return 1;
  }
  else{
    inputPHSD = argv[1];
    filename  = argv[2];
    f791name  = argv[3];
  }

  TFile *OUTfile = new TFile(argv[4], "RECREATE");
  OUTfile -> SetCompressionLevel(7);

// ---------------------- Run header (inputPHSD) -------------------------------
  TTree *input    = new TTree("inputPHSD", "Run header");
  std::vector<Run> *brun = new std::vector<Run>;
  input -> Branch("run", &brun);

  int n_lines = 0;
  int PHSD_INPUT_VERSION = 4;

  fgzFile = gzopen(inputPHSD, "rb");
  if(!fgzFile){
    printf("-E- Can not open file: %s\n", inputPHSD);
    exit(-1);
  }

  while (gzgets(fgzFile, fbuffer, sizeof(fbuffer)) != NULL){
    ++n_lines;
  }
  gzrewind(fgzFile);
  Run *run = new Run();

  if     (n_lines == 36) PHSD_INPUT_VERSION = 5;
  else if(n_lines == 39) PHSD_INPUT_VERSION = 6;
  printf("PHSD input file format: %2d\n", PHSD_INPUT_VERSION);
  int   ival;
  float fval;
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetAtar(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetZtar(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetAproj(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetZproj(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%f", &fval); run -> SetTkin(fval);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%f", &fval); run -> SetBmin(fval);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%f", &fval); run -> SetBmax(fval);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%f", &fval); run -> SetDeltaB(fval);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetNum(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetSub(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetSeed(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetIGLUE(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%f", &fval); run -> SetFinalT(fval);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetIlow(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetIdil(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetICQ(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetIHARD(ival);
  if(PHSD_INPUT_VERSION > 4) { // IDQPM -- new version of the PHQMD
    gzgets(fgzFile, fbuffer, 256);
    sscanf(fbuffer, "%d", &ival);
    run -> SetIDQPM(ival);
  }
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetIBW(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetIuser(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetINUCLEI(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetIPHQMD(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetISACA(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%f", &fval); run -> SetTSACA(fval);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%f", &fval); run -> SetDTSACA(fval);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetNTSACA(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetFLGSACA(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetYuk(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetAsy(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetPair(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetCoul(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%f", &fval); run -> SetAsy0(fval);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%f", &fval); run -> SetEPair(fval);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetEOS(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetResSACA(ival);
  gzgets(fgzFile, fbuffer, 256); sscanf(fbuffer, "%d", &ival); run -> SetWigDens(ival);
  gzclose(fgzFile);


// ------------------------ phsd.dat stuff -----------------------------
  int    Nparticles, evSub, evNum;
  int    chrg, chnl;
  double pB;
  double px,  py,  pz, pe;
  double fx,  fy,  fz, ft;
  double fpx, fpy, fpz;
  double bd, ed;
  long   PDG;
  int    PHSD_UNIQ_ID;
  int    Nparticipants;
  std::array<float, 4> psi;
  std::array<float, 4> ecc;

  std::vector<Event> *events = new std::vector<Event>;

  printf("%s: ", filename);
  fgzFile  = gzopen(filename, "rb");
  if(!fgzFile){
    printf("-E- Can not open file: %s\n", filename);
    exit(1);
  }

  int    EventNr = 1;
  while (gzgets(fgzFile, fbuffer, 512) != NULL){
    if (gzeof(fgzFile)){
      gzclose(fgzFile);
      fgzFile=NULL;
      break;
    }
    Event *evt = new Event();
    sscanf(fbuffer, "%d %d %d %lf", &Nparticles, &evSub, &evNum, &pB);
    gzgets(fgzFile, fbuffer, 512);
    sscanf(fbuffer, "%d %f %f %f %f %f %f %f %f", &Nparticipants, &psi[0], &ecc[0], &psi[1], &ecc[1], &psi[2], &ecc[2], &psi[3], &ecc[3]);
    evt -> SetNparticipants(Nparticipants);
    evt -> SetNparticles(Nparticles);
    evt -> SetB(pB);
    evt -> SetNum(evNum);
    evt -> SetSub(evSub);
    evt -> SetPsi(psi);
    evt -> SetEcc(ecc);
    for (int i = 0; i < Nparticles; ++i){
      gzgets(fgzFile, fbuffer, 512);
      sscanf(fbuffer, "%ld %d %lf %lf %lf %lf %d %*d %d %*f %*f %*f %lf %lf %lf %lf %lf %lf %lf %*f %lf %lf",
                       &PDG, &chrg, &px, &py, &pz, &pe, &chnl, &PHSD_UNIQ_ID, &fx, &fy, &fz, &ft, &fpx, &fpy, &fpz, &bd, &ed);
      Particle *part = new Particle();
      part -> SetId(PHSD_UNIQ_ID);
      part -> SetPdg(PDG);
      part -> SetCharge(chrg);
      part -> SetChannel(chnl);
      part -> SetBarDens(bd);
      part -> SetEnDens(ed);
      part -> SetP(ROOT::Math::PxPyPzEVector(px, py, pz, pe));
      part -> SetR(ROOT::Math::XYZTVector(fx, fy, fz, ft));
      evt  -> AddParticle(*part);
    }
    events -> push_back(*evt);
    ++EventNr;
  }
  gzclose(fgzFile);
// At this step all events from the phsd.dat file are stored
// in the 'std::vector<Event> *events'


  const int NumOutEvents   = events -> size();
  printf(" %d events\n", NumOutEvents);
  run   -> SetNevents(NumOutEvents);
  brun  -> push_back(*run);
  input -> Fill();
  brun  -> clear();

// ------------------------ fort.791(891) stuff -----------------------------
  std::map<int, std::vector<int>> steps_map;
  steps_map.insert(std::pair<int, std::vector<int> >(NumOutEvents, std::vector<int>()));

  std::map<int, std::map<int, std::map<int, std::vector<int>>>> tree_clusters_map;
  std::map<int, std::map<int, std::map<int, std::vector<int>>>> tree_clusters_map_mstbar;
//  std::map<int, std::map<int, std::vector<int>>> timestep_clusters_map;
//  std::map<int, std::vector<int>> event_clusters_map;

  int   NUM    = run ->    GetNum();
  int   SUB    = run ->    GetSub();
  float TSACA  = run ->  GetTSACA();
  float DTSACA = run -> GetDTSACA();
  int   NTSACA = run -> GetNTSACA();

  std::vector<std::vector<Baryon>> *timesteps = new std::vector<std::vector<Baryon>>;
  std::vector<Baryon> *baryons = new std::vector<Baryon>;
  Baryon *baryon;

  printf("%s: ", f791name);
  fflush(stdout);
  f791File = gzopen(f791name, "rb");
  if(!f791File){
    printf("-E- Can not open file: %s\n", f791name);
    exit(1);
  }

  int Nbaryons = 0;
  int curNUM = 0;
  int curSUB = 0;
  int curT  = 0;

  double x, y, z, t;
  double m;
  int Bid;
  int FragId, FragSize;
  int ProdRegion, ProdChannel;
  double ProdTime, EbindN;
  int curEntry = 0;
  while (gzgets(f791File, fbuffer, 512) != NULL){
    if (gzeof(f791File)){
      gzclose(f791File);
      f791File=NULL;
      break;
    }
    sscanf(fbuffer, "%d %d %*f %*f %d", &curNUM, &curSUB, &curT);
    int curEvt = curNUM + NUM*(curSUB-1) -1;
    Event *event = &events -> at(curEvt);
    gzgets(f791File, fbuffer, 512);
    gzgets(f791File, fbuffer, 512);
    sscanf(fbuffer, "%d", &Nbaryons);
    for (int i = 0; i < Nbaryons; ++i){
      baryon = new Baryon();
      gzgets(f791File, fbuffer, 512);
      //                  bid px  py  pz  x   y   z   m  fid fs unq       reg ch  t  ebnd
      sscanf(fbuffer, "%*d %d %lf %lf %lf %lf %lf %lf %lf %d %d %d %*d %*d %d %d %lf %lf",
                           &Bid,
                           &px, &py, &pz,
                           &x, &y,  &z,
                           &m, &FragId, &FragSize, &PHSD_UNIQ_ID,
                           &ProdRegion, &ProdChannel, &t, &EbindN);
      baryon -> SetPhsdId(PHSD_UNIQ_ID);
      baryon -> SetP(ROOT::Math::PxPyPzMVector(px, py, pz, m));
      baryon -> SetR(ROOT::Math::XYZTVector(x, y, z, t));
      baryon -> SetFid(FragId);
      baryon -> SetFsize(FragSize);
      baryon -> SetProdRegion(ProdRegion);
      baryon -> SetProdChannel(ProdChannel);
      baryon -> SetEbind(EbindN);
      std::vector<Particle> *particles = event -> GetParticles();
      for (int iterPart = 0; iterPart < event -> GetNparticles(); ++iterPart){
        if (particles->at(iterPart).GetPdg() > 999   && // Mesons should not go to the baryons vector
            particles->at(iterPart).GetPdg() < 9999  && //   as well as deuterons (kinetic d pdg = 100121)
            particles->at(iterPart).GetId() == PHSD_UNIQ_ID){
          baryon -> SetPdg(particles->at(iterPart).GetPdg());
          baryon -> SetPhsdPos(iterPart);
          if(! particles->at(iterPart).IsInMst()) particles->at(iterPart).SetInMst(1);
          break;
        }
        else{
          baryon -> SetPdg(0);
          baryon -> SetPhsdPos(-1);
        }
      }
      if(baryon -> GetPdg() == 0){ // PHSD particle does not exist in the final phsd.dat file!
        baryon -> SetPhsdPos(-1);
        if      (Bid ==  1) baryon -> SetPdg(2212);
        else if (Bid ==  0) baryon -> SetPdg(2212);
        else if (Bid == 17) baryon -> SetPdg(3122);
        else if (Bid == 18) baryon -> SetPdg(3212);
      }
      if(FragSize > 1 && EbindN < -1E-5){
        tree_clusters_map_mstbar[curEvt][curT][FragId].push_back(i);
        if(baryon -> GetPhsdPos() >= 0){
          tree_clusters_map[curEvt][curT][FragId].push_back(baryon -> GetPhsdPos());
        }
      }
      baryons  -> push_back(*baryon);
    }
    steps_map[curEvt].push_back(curEntry);
    timesteps -> push_back(*baryons);
    baryons -> clear();
    ++curEntry;
  }
  gzclose(f791File);
  printf(" DONE\n");



// -----------------------------------------------------------------------
// ------------------------ final tree stuff -----------------------------
// -----------------------------------------------------------------------
  printf("Final tree construction: ");
  fflush(stdout);
  TTree *tree_final = new TTree("Events", "PHQMD events tree");
  std::vector<Event> *events_final = new std::vector<Event>;
  tree_final -> Branch("event", &events_final);

  for (int iterEvt = 0; iterEvt < NumOutEvents; ++iterEvt){
    Event *event = &events -> at(iterEvt);
    //~ printf("Event: %4d\n", iterEvt);
    for (int iterStep = 0; iterStep < steps_map[iterEvt].size(); ++iterStep){
      //~ if (iterEvt == 322) printf(" ----- Step: %4d, stored clusters: %d\n", iterStep, event -> GetClusters() -> size());
      std::vector<Fragment> mstfstep;
      std::vector<Fragment> stable_clusters_array;
      event -> AddMstBStep(timesteps -> at(steps_map[iterEvt][iterStep]));

      std::map<int, std::vector<int>> event_clusters_map_mstbar = tree_clusters_map_mstbar[iterEvt][iterStep+1];
      for (auto const& [FID, val] : event_clusters_map_mstbar){
        Fragment cl;
        int  cZ   = 0;
        int  cN   = 0;
        int  cLS  = 0;
        double  M = 0;
        double  X = 0;
        double  Y = 0;
        double  Z = 0;
        double  T = 0;
        ROOT::Math::PxPyPzMVector cl_p;
        for (int cl_bar = 0; cl_bar < val.size(); ++cl_bar){
          Baryon    cbar = event -> GetMstBSteps().at(iterStep).at(val[cl_bar]);
          int bpdg = cbar.GetPdg();
          cl_p += cbar.GetP();
          if     (bpdg == 2212) ++cZ;
          else if(bpdg == 2112) ++cN;
          else if(bpdg == 3122 || bpdg == 3212) ++cLS;
          X += cbar.GetR().X() * cbar.GetP().M();
          Y += cbar.GetR().Y() * cbar.GetP().M();
          Z += cbar.GetR().Z() * cbar.GetP().M();
          M += cbar.GetP().M();
          T = cbar.GetR().T();
        }
        cN += cLS;
        int cA = cN + cZ;
        cl.SetPdg(10*10E7 + cLS*10E6 + cZ*10E3 + cA*10);
        cl.SetFid(FID);
        cl.SetP(cl_p);
        X /= M;
        Y /= M;
        Z /= M;
        cl.SetR(ROOT::Math::XYZTVector(X, Y, Z, T));
        mstfstep.push_back(cl);
      }

      // Stable (final) clusters after new 'stabilization' method
      std::map<int, std::vector<int>> event_clusters_map = tree_clusters_map[iterEvt][iterStep+1];
      for (auto const& [FID, val] : event_clusters_map){
        std::vector<int> phsd_baryons;
        int  cZ   = 0;
        int  cN   = 0;
        int  cLS  = 0;
        double TF = 0;
        bool stbl = 1;
        for (int cl_bar = 0; cl_bar < val.size(); ++cl_bar){
          int phsd_position = val[cl_bar];
          Particle cpart;
          if(phsd_position == -1){
            stbl = 0;
          }
          if(stbl){
            cpart = event -> GetParticles()->at(val[cl_bar]);
            phsd_baryons.push_back(val[cl_bar]);
          }

          int ppdg = cpart.GetPdg();
          if     (ppdg == 2212) ++cZ;
          else if(ppdg == 2112) ++cN;
          else if(ppdg == 3122 || ppdg == 3212) ++cLS;

          double frz_time = -1.0;
          if(stbl){
            frz_time = cpart.GetR().T();
            if(frz_time > TSACA + DTSACA * iterStep){
              stbl = 0;
            }
            if(frz_time > TF){
              TF = frz_time;
            }
          }
        } // loop over baryons in the current cluster

        if(!stbl) continue;

        cN += cLS;
        int cA = cN + cZ;

        std::vector<Fragment> *stored_clusters = event -> GetClusters();
        bool already_stored = 0;
        std::vector<int> current_structure = val;
        sort(current_structure.begin(), current_structure.end());
        double match = 0.0;

        std::vector<int> marked_for_removal;
          for (int iterClust = 0; iterClust < event -> GetClusters() -> size(); iterClust++){
            Fragment stored_cluster = event -> GetClusters() ->at(iterClust);
            std::vector<int> stored_structure = stored_cluster.GetStructure();
            sort(stored_structure.begin(), stored_structure.end());

            std::vector<int> intersection;
            std::set_intersection(
              stored_structure.begin(), stored_structure.end(),
              current_structure.begin(), current_structure.end(),
              std::back_inserter(intersection));
            if(intersection.size() == 0) continue;

            double match = (double)intersection.size() / (double)cA;
            if (match > 0 && current_structure.size() > stored_structure.size()){
              marked_for_removal.push_back(iterClust);
            }
            else if (match > 0 && current_structure.size() <= stored_structure.size()){
              already_stored = 1;
              continue;
            }
          } // loop over stored clusters

        if(already_stored) continue;

        if(marked_for_removal.size() > 0){
          std::sort(marked_for_removal.begin(), marked_for_removal.end(), std::greater{});
          for(int marked = 0; marked < marked_for_removal.size(); marked++){
            event -> RemoveCluster(marked_for_removal[marked]);
          }
        }
        marked_for_removal.clear();

        Fragment *stable_cluster = new Fragment();
        ROOT::Math::PxPyPzMVector cl_pf;
        stable_cluster -> SetFid(FID);
        stable_cluster -> SetStructure(current_structure);
        double XF = 0;
        double YF = 0;
        double ZF = 0;

        for (int cl_bar = 0; cl_bar < current_structure.size(); ++cl_bar){
          Particle cpart = event -> GetParticles()->at(current_structure[cl_bar]);
          double bx = cpart.GetP().Px() / cpart.GetP().E();
          double by = cpart.GetP().Py() / cpart.GetP().E();
          double bz = cpart.GetP().Pz() / cpart.GetP().E();
          double dt = TF - cpart.GetR().T();
          double x  = cpart.GetR().X() + bx * dt;
          double y  = cpart.GetR().Y() + by * dt;
          double z  = cpart.GetR().Z() + bz * dt;
          XF += x * cpart.GetP().M();
          YF += y * cpart.GetP().M();
          ZF += z * cpart.GetP().M();
          cl_pf +=  cpart.GetP();
        }
        XF /= cl_pf.M();
        YF /= cl_pf.M();
        ZF /= cl_pf.M();
        stable_cluster -> SetPdg(10*10E7 + cLS*10E6 + cZ*10E3 + cA*10);
        stable_cluster -> SetFreezeoutR(ROOT::Math::XYZTVector(XF, YF, ZF, TF));
        stable_cluster -> SetFreezeoutP(cl_pf);
        event -> AddCluster(*stable_cluster);
        phsd_baryons.clear();
      }// loop over current clusters
      event -> AddMstFStep(mstfstep);
      mstfstep.clear();
    } // loop over steps
    std::vector<Fragment> *stored_clusters = event -> GetClusters();
    events_final -> push_back(*event);
    tree_final -> Fill();
    events_final -> clear();
  }
  printf(" DONE\n");

  OUTfile -> Write();
  OUTfile -> Close();
  return 0;
}
