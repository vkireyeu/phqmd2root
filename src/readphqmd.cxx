#include "PHQMDevent.h"

using namespace PHQMD;

int main(int argc, char* argv[]){
  if(argc < 2){
    printf("usage: %s <in.root>\n", argv[0]);
    return 1;
  }

  TFile *INfile = new TFile(argv[1], "READ");

  TTree *tree         = (TTree*) INfile -> Get("Events");
  std::vector<Event> *event = new std::vector<Event>;
  tree -> SetBranchAddress("event", &event);

  int NEVENTS = tree -> GetEntries();
  printf("Tree: %d events\n", NEVENTS);

  for (int iterEvt = 0; iterEvt < NEVENTS; ++iterEvt){
    tree -> GetEntry(iterEvt);
    Event *ev = &event -> at(0);
    printf("%4d %8.3f  %4d  %4d\n", iterEvt, ev -> GetB(), ev -> GetNparticles(), ev -> GetNparticipants());
    std::vector<Particle> *particles = ev -> GetParticles();
    
    //~ std::vector<std::vector<Baryon>>   MstBSteps = ev -> GetMstBSteps();
    //~ for (int iterStep = 0; iterStep < MstBSteps.size(); ++iterStep){
      //~ printf(" (B) Step %d\n", iterStep);
      //~ std::vector<Baryon> baryons = MstBSteps[iterStep];
      //~ for (int iterBar = 0; iterBar < baryons.size(); ++iterBar){
        //~ if(baryons[iterBar].GetPhsdPos() > -1){
          //~ printf("%d(%d)=%d ", baryons[iterBar].GetPhsdId(), particles -> at(baryons[iterBar].GetPhsdPos()).GetId(), baryons[iterBar].IsBound() );
        //~ }
        //~ else{
          //~ printf("%d(%d)=%d ", baryons[iterBar].GetPhsdId(), baryons[iterBar].GetPhsdPos(), baryons[iterBar].IsBound() );
        //~ }
      //~ }
      //~ printf("\n");
    //~ }
    
    std::vector<std::vector<Fragment>> MstFSteps = ev -> GetMstFSteps();
    for (int iterStep = 0; iterStep < MstFSteps.size(); ++iterStep){
      printf(" (F) Step %d\n", iterStep);
      std::vector<Fragment> clusters = MstFSteps[iterStep];
      for (int iterCls = 0; iterCls < clusters.size(); ++iterCls){
        printf("%d[%d:%d]=%f ", clusters[iterCls].GetPdg(), clusters[iterCls].GetFid(), clusters[iterCls].IsStable(), clusters[iterCls].GetFreezeoutR().T());
      }
      printf("\n");
    }
  }


  INfile -> Close();
  return 0;
}
