#ifdef __CINT__

#pragma link C++ class PHQMD::Particle +;
#pragma link C++ class std::vector<PHQMD::Particle> +;

#pragma link C++ class PHQMD::Baryon +;
#pragma link C++ class std::vector<PHQMD::Baryon> +;
#pragma link C++ class std::vector<std::vector<PHQMD::Baryon>> +;

#pragma link C++ class PHQMD::Fragment +;
#pragma link C++ class std::vector<PHQMD::Fragment> +;
#pragma link C++ class std::vector<std::vector<PHQMD::Fragment>> +;

#pragma link C++ class PHQMD::Event +;
#pragma link C++ class std::vector<PHQMD::Event> +;

#pragma link C++ class PHQMD::Run +;
#pragma link C++ class std::vector<PHQMD::Run> +;

#endif
