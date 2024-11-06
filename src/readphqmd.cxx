#include "PHQMDevent.h"

using namespace PHQMD;

int main(int argc, char* argv[]){
  if(argc < 3){
    printf("usage: %s <in.root> <evt>\n", argv[0]);
    return 1;
  }

  TFile *INfile = new TFile(argv[1], "READ");
  int nevt = atoi(argv[2]);

  TTree *tree         = (TTree*) INfile -> Get("Events");
  std::vector<Event> *event = new std::vector<Event>;
  tree -> SetBranchAddress("event", &event);

  int NEVENTS = tree -> GetEntries();
  printf("Tree: %d events\n", NEVENTS);


// ----------------------------- Header ---------------------------------
  TTree *run = (TTree*) INfile -> Get("inputPHSD");
  std::vector<Run> *headers = new std::vector<Run>;
  run -> SetBranchAddress("run", &headers);
  run -> GetEntry(0);
  Run *header = &headers -> at(0);

  float TSACA  = header ->  GetTSACA();
  float DTSACA = header -> GetDTSACA();

  tree -> GetEntry(nevt);
  Event *ev = &event -> at(0);
  printf("%4d %8.3f  %4d  %4d\n", nevt, ev -> GetB(), ev -> GetNparticles(), ev -> GetNparticipants());
  std::vector<Particle> *particles = ev -> GetParticles();
  int phsd_size = particles -> size();

  std::vector<std::vector<Fragment>> MstFSteps = ev -> GetMstFSteps();
  std::vector<std::vector<Baryon>>   MstBSteps = ev -> GetMstBSteps();
  for (int iterStep = 0; iterStep < MstFSteps.size(); ++iterStep){
    printf("---------------------------------------------------------\n");
    printf("Step %3d, time = %8.2f fm/c\n", iterStep, TSACA + iterStep*DTSACA);
    printf("---------------------------------------------------------\n");
    std::vector<Fragment> clusters = MstFSteps[iterStep];
    std::vector<Baryon>   baryons  = MstBSteps[iterStep];
    
    for (int iterCls = 0; iterCls < clusters.size(); ++iterCls){
      std::vector<int> structure = clusters[iterCls].GetStructure();
      printf("%3d       %d\n", iterCls, clusters[iterCls].GetPdg());
      for(int iterBar = 0; iterBar < structure.size(); ++iterBar){
        int phsd_pos = baryons[structure[iterBar]].GetPhsdPos();
        if(phsd_pos == -1){
          printf("%d [%6d], Does not exist in the final PHSD event!\n", 
                  baryons[structure[iterBar]].GetPdg(),
                  baryons[structure[iterBar]].GetPhsdId());
        }
        else{
          printf("%d [%6d], t_freezeout = %8.2f,     px = %9.5f\n", 
                  particles -> at(phsd_pos).GetPdg(),
                  particles -> at(phsd_pos).GetId(),
                  particles -> at(phsd_pos).GetR().T(),
                  particles -> at(phsd_pos).GetP().Px());
        }
      }
    }
    printf("\n");
  }


  printf("\n");
  printf("\n");
  printf("\n");
  printf("\n");
  printf("---------------------------------------------------------\n");
  printf("------------------ V.K. Stabilization -------------------\n");
  printf("---------------------------------------------------------\n");
  for(int cls = 0; cls < ev -> GetClusters() -> size(); ++cls){
    Fragment cluster = ev -> GetClusters() -> at(cls);
    printf("%3d       %d\n", cls, cluster.GetPdg());
    std::vector<int> structure = cluster.GetStructure();
    for(int iterBar = 0; iterBar < structure.size(); ++iterBar){
      printf("%d [%6d], t_freezeout = %8.2f,     px = %9.5f\n", 
              particles -> at(structure[iterBar]).GetPdg(),
              particles -> at(structure[iterBar]).GetId(),
              particles -> at(structure[iterBar]).GetR().T(),
              particles -> at(structure[iterBar]).GetP().Px());
    }
  }



  INfile -> Close();
  return 0;
}
