#include "PHQMDevent.h"

using namespace PHQMD;

struct pInfo{
  int    pdg;
  int charge;
  double  px;
  double  py;
  double  pz;
  double   e;
  double   x;
  double   y;
  double   z;
  double   t;
  int   uniq;
};

int main(int argc, char* argv[]){
  const char*  input_file_name;
  const char* output_file_name;
  if(argc < 3){
    printf("usage: %s <input.root> <output.dat>\n", argv[0]);
    return 1;
  }
  else{
    input_file_name  = argv[1];
    output_file_name = argv[2];
  }

// ---------------------------------------------------------------------
// Input file check
  std::unique_ptr<TFile> inFile(TFile::Open(input_file_name, "READ"));
  if (!inFile || inFile -> IsZombie()) {
    printf("-MST2DAT ERR: Input file: %s does not exist\n", input_file_name);
    return 3;
  }
  else printf("-MST2DAT OK: Input file %s\n", input_file_name);


// ----------------------------- ASCII ---------------------------------
  FILE *outFile;                                            // ASCII file for 'phsd.dat'-like output
  if ((outFile = fopen(output_file_name, "w"))==NULL) {
    printf("-MST2DAT ERR: Cannot open output file %s.\n", output_file_name);
    exit(1);
  }

// ----------------------------- Header ---------------------------------
  TTree *run = (TTree*) inFile -> Get("inputPHSD");
  std::vector<Run> *headers = new std::vector<Run>;
  run -> SetBranchAddress("run", &headers);
  run -> GetEntry(0);
  Run *header = &headers -> at(0);

// ------------------------ Main events loop ----------------------------
  TTree *tree = (TTree*) inFile -> Get("Events");
  std::vector<Event> *event = new std::vector<Event>;
  tree -> SetBranchAddress("event", &event);

  int NEVENTS = tree -> GetEntries();
  printf("Tree: %d events\n", NEVENTS);

  for (int iterEvt = 0; iterEvt < NEVENTS; ++iterEvt){
    tree -> GetEntry(iterEvt);
    Event *ev = &event -> at(0);
    //~ printf("%4d %8.3f  %4d  %4d\n", iterEvt, ev -> GetB(), ev -> GetNparticles(), ev -> GetNparticipants());

    int evNparticipants = ev -> GetNparticipants();
    std::array<float,4> evPsi = ev -> GetPsi();
    std::array<float,4> evEcc = ev -> GetEcc();

    //~ std::vector<Fragment> *clusters = ev -> GetClusters();
    std::vector<pInfo> vClusters;
    std::vector<int> bound_baryons;
    std::vector<int> particles_in_mst;

    for(int cls = 0; cls < ev -> GetClusters() -> size(); ++cls){
      Fragment cluster = ev -> GetClusters() -> at(cls);
      vClusters.push_back({cluster.GetPdg(), cluster.GetZ(),
                          cluster.GetFreezeoutP().Px(), cluster.GetFreezeoutP().Py(), cluster.GetFreezeoutP().Pz(),
                          cluster.GetFreezeoutP().E(),
                          cluster.GetFreezeoutR().X(), cluster.GetFreezeoutR().Y(), cluster.GetFreezeoutR().Z(),
                          cluster.GetFreezeoutR().T(),
                          cls});
      std::vector<int> cluster_structure = cluster.GetStructure();
      bound_baryons.insert(bound_baryons.end(), cluster_structure.begin(), cluster_structure.end());
    }


    std::vector<std::vector<Baryon>> MstBSteps = ev -> GetMstBSteps();
    int last_step = MstBSteps.size() - 1;
    std::vector<Baryon> baryons = MstBSteps[last_step];
    std::vector<pInfo> vBaryons;
    for (int bar = 0; bar < baryons.size(); ++bar){
      Baryon baryon = baryons[bar];
      int phsd_position = baryon.GetPhsdPos();
      if(std::find(bound_baryons.begin(), bound_baryons.end(), phsd_position) == bound_baryons.end()) {
        Particle particle = ev -> GetParticles() -> at(phsd_position);
        vBaryons.push_back({particle.GetPdg(), particle.GetCharge(),
                            particle.GetP().Px(), particle.GetP().Py(), particle.GetP().Pz(),
                            particle.GetP().E(),
                            particle.GetR().X(), particle.GetR().Y(), particle.GetR().Z(),
                            particle.GetR().T(),
                            particle.GetId()});
        particles_in_mst.push_back(phsd_position);
      }
    }

    std::vector<Particle> *particles = ev -> GetParticles();
    std::vector<pInfo> vParticles;
    for(int part = 0; part < particles -> size(); ++part){
      Particle particle = particles -> at(part);
      //~ if(particle.IsInMst()) continue;
      //~ if(particle.IsInMst()) continue;
      if(std::find(bound_baryons.begin(), bound_baryons.end(), part) != bound_baryons.end()) continue;
      if(std::find(particles_in_mst.begin(), particles_in_mst.end(), part) != particles_in_mst.end()) continue;
      vParticles.push_back({particle.GetPdg(), particle.GetCharge(),
                            particle.GetP().Px(), particle.GetP().Py(), particle.GetP().Pz(),
                            particle.GetP().E(),
                            particle.GetR().X(), particle.GetR().Y(), particle.GetR().Z(),
                            particle.GetR().T(),
                            particle.GetId()});
    }
    bound_baryons.clear();
    particles_in_mst.clear();

    int NOutPart = vClusters.size() + vBaryons.size() + vParticles.size();

// ----- Constructing the full event
    fprintf(outFile, "%6d %10d %6d  %10E  %f  1\n",
            NOutPart, ev -> GetSub(), ev -> GetNum(),
            ev -> GetB(), header -> GetTkin());  // Event header filling - 2 lines

    fprintf(outFile, "%6d %E %E %E %E %E %E %E %E\n", evNparticipants,
                                                      evPsi[0], evEcc[0],
                                                      evPsi[1], evEcc[1],
                                                      evPsi[2], evEcc[2],
                                                      evPsi[3], evEcc[3]);

    for (auto& it :  vClusters) fprintf(outFile, "%12d %6d % 3.8E % 3.8E % 3.8E % 3.8E % 3.8E % 3.8E % 3.8E % 3.8E %12d \n", it.pdg, it.charge, it.px, it.py, it.pz, it.e, it.x, it.y, it.z, it.t, it.uniq);
    for (auto& it :   vBaryons) fprintf(outFile, "%12d %6d % 3.8E % 3.8E % 3.8E % 3.8E % 3.8E % 3.8E % 3.8E % 3.8E %12d \n", it.pdg, it.charge, it.px, it.py, it.pz, it.e, it.x, it.y, it.z, it.t, it.uniq);
    for (auto& it : vParticles) fprintf(outFile, "%12d %6d % 3.8E % 3.8E % 3.8E % 3.8E % 3.8E % 3.8E % 3.8E % 3.8E %12d \n", it.pdg, it.charge, it.px, it.py, it.pz, it.e, it.x, it.y, it.z, it.t, it.uniq);
  }


  inFile -> Close();
  fclose(outFile);
  return 0;
}
