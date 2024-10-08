[TOC]

# PHQMDEvent
A new ROOT TTree based format to store the PHQMD events

<p align="center">
   <img src="https://github.com/vkireyeu/phqmd2root/blob/main/mst_evolution_1fm.gif">
   <img src="https://github.com/vkireyeu/phqmd2root/blob/main/mst_evolution_6fm.gif">
   <img src="https://github.com/vkireyeu/phqmd2root/blob/main/mst_evolution_10fm.gif">
</p>

- [Tree structure](#tree)
- [Implemended classes](#classes)
- [Compilation](#compile)
- [Usage](#usage)
- [Documentation](#documentation)

## Tree structure {#tree}
```
Events
|--event[fNevents]
   |-- fB
   |-- fNparticles
   |-- fNparticipants
   |-- fPsi[4]
   |-- fEcc[4]
   |-- particles[fNparticles]
       |-- fID
       |-- fPDG
       |-- fCharge
       |-- fChannel
       |-- fBarDens
       |-- fEnDens
       |-- fIsInMst
       |-- fR[4] (X, Y, Z, T)
       |-- fP[4] (Px, Py, Pz, E)
   |-- mstbsteps[fNTSACA+1]
       |-- fPHSD_POS
       |-- fPDG
       |-- fPHSD_ID
       |-- fFID
       |-- fFSIZE
       |-- fPREG
       |-- fPCH
       |-- fEBND
       |-- fR[4] (X, Y, Z, T)
       |-- fP[4] (Px, Py, Pz, M)
   |-- mstfsteps[fNTSACA+1]
       |-- fFID
       |-- fPDG
       |-- fR[4] (X, Y, Z, T)
       |-- fP[4] (Px, Py, Pz, M)
inputPHSD
|--run
   |-- fNevents
   |-- fAtar
   |-- fZtar
   |-- fAproj
   |-- fZproj
   |-- fTkin
   |-- fBmin
   |-- fBmax
   |-- fDeltaB
   |-- fIGLUE
   |-- fIBW
   |-- fNUM
   |-- fSUB
   |-- fSEED
   |-- fFinalT
   |-- fILOW
   |-- fIdil
   |-- fICQ
   |-- fIHARD
   |-- fIDQPM
   |-- fIUSER
   |-- fINUCLEI
   |-- fIPHQMD
   |-- fISACA
   |-- fTSACA
   |-- fDTSACA
   |-- fNTSACA
   |-- fFLGSACA
   |-- fYuk
   |-- fAsy
   |-- fPair
   |-- fCoul
   |-- fAsy0
   |-- fEPair
   |-- fEOS
   |-- fResSACA
   |-- fWigDens
```


## Implemended classes {#classes}
- \ref PHQMD::Particle()
- \ref PHQMD::Baryon()
- \ref PHQMD::Fragment()
- \ref PHQMD::Event()
- \ref PHQMD::Run()

## Compilation {#compile}
In the source directory:  
```
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local
make
make install
```
`-DCMAKE_INSTALL_PREFIX` is an optional argument to the CMake.


## Usage {#usage}
The `phqmd2root` converter is shipped within the source code.  
To convert the PHQMD ASCII files to the single ROOT file use this program:  
``` phqmd2root <inputPHSD> <phsd.dat.gz> <fort.791> <out.root> ```  
Where `<out.root>` is the name of your new ROOT file.

The `readphqmd` program is an example of how to access and read the stored information:  
``` readphqmd <out.root> ```  


## Documentation {#documentation}
To generate the PHQMDEvent classes documentation one can use the Doxygen genetator:
``` doxygen doxygen_config ```  
The generated documentation will be placed in the newly created folder: ```html/index.html```.
