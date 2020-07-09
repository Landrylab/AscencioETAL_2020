% % Process excel files from the fluorometry competition experiments
%%Calculate fluorescence ratios and selectioncoefficients
% Fill Input arguments

Oneset=1; %1, one set of plates, 0, two set of plates. USER 
addpath('..\..\..\SCRIPTOMA\')
%-Folder with the excel files -
expdirectory='MoBYall\'; %USER

%-Specify number of plates {mutant,reference}-
pltCell={[1:26]}; %--USER--
mut_strain=['rfp']; %--USER--
ref_strain=['cfp']; %--USER--
mutRefCell={{mut_strain,ref_strain}}; %first mut, then ref.

%-Datos para sustituir nombres de variables en las estructuras-
newfields{1}='od';
newfields{2}='rfp'; %--USER-- 
newfields{3}='cfp'; %--USER--
%newfields{4}='yfp'; %--USER-- optional, uncomment if you were use 3 fluorophores
newfields{4}='t'; %--USER--, only check the number
newfields{5}='id'; %--USER--, only check the number


%-Specify type of plate (96, microplate or 384, omnytray)- 
pltSize=96; %USER

% 2. Get raw data. 

%Parsea\ all data from excel files
ParsingDatosA;
bgdataA=bgdata;
save rawA bgdataA
clear bgdata
if Oneset==0
    mutRefCell={{mut_strain,ref_strain}}; %primero mut, luego ref.
    ParsingDatosB;
    bgdataB=bgdata;
    save rawB bgdataB
    clear bgdata
    bgdataAll=combineData(bgdataA,bgdataB);
    save rawAll bgdataAll       
else 
    bgdataAll=bgdataA;

save rawAll bgdataAll
  end
%%
%load rawAll
%close all
plsA=[1:11];%Plates SET1 Nominal
plsB=[12:13];%Plates SET2 NaCl
wrpl=[12,91,63,27,55,80];%well with ref\ref in each well
ww=1:6; wref=[ww ww+12 ww+24 ww+36 ww+48 ww+60 ww+72 ww+84];%plate ref\ref
alli=1:96;
xi=setdiff(alli,wrpl);
%bkgFile1='bkgFile_nominal.xls';bkgFile2='bkgFile_nominal.xls';
bkgFile1='bkgFile_nacl.xls';bkgFile2='bkgFile_nacl.xls';
%% General Parameters. 
% OD y and time used for calculate R/C ratios
useOd = 0.2
useT=10;
% Signal threshold R\C
detectTol=4;
%Dilution factor applied each competition day.
%dilFactor = [32 32 32 32 32 32 32];
dilFactor = [16 16 16 16 16 16]
longev=0;%1 longevity, 0 Competition
pltSize=96;
%% Assingn names to data
nameFile1='CarrCollectionMobyAll.xls'; %maps experiment plates to collection plates
nameFile2='MOBY_all_genenames.xlsx'; %maps collection plates to gene names.
bgdataAll=AssignNames(bgdataAll,nameFile1,nameFile2); % OPCIONAL 1 BKG PLT

%% Calculate BKG and normalize raw data. 
[polRefBg1,polMutBg1,polRefBg2,polMutBg2]=backgroundCurvesAJG2_2exp(bgdataAll,bkgFile1,bkgFile2);%saca 4 polinomios 2 para GABA y 2 Para Glut
bgdataAll=bkgSubtract_2exp(bgdataAll,polRefBg1,polMutBg1,polRefBg2,polMutBg2,plsA, plsB);
%% Estimate RoC useOd
bgdataAll=logarithmicCocientAJG(bgdataAll,useOd,detectTol,longev);%Competencia

%% Estimate S for competition
bgdataAll=selectionCoefficients(bgdataAll,dilFactor,longev)%Competencias

%% Robust Coeff (s4) y error ajuste
%bgdataAll=selectionCoefficientsAJG3Robustfit(bgdataAll,dilFactor,longev);
%
%% 16 post analysis: RoC and s3,s4 plate by plate (96 wells per figure.)
longevityPostNoSWAP_DIAS(bgdataAll,8,12,[6]);

%%
%bgdataAll=SpecificFixesGAB(mutbkg,refbkg,bgdataAll);

%% Histograms RawS
[rawsNominal,rawsNaCl]=longevityPostNoSWAP2wellsRaw(bgdataAll,plsA, plsB);%Datos KO, no incluye NaNs

[fnominal,bnominal,fncl,bncl]=rawdata_pdf(rawsNominal,rawsNaCl, 10); % The final number is the cut off
%% Data normalization

%bgdataAll=exploreresultsGabaGlut1_dcutoff(bgdataAll,plsA,plsB,wref,wrpl,30,xi);%wref=placa control, wrpl == ref per plate, plc=placas control

save dataAll_nacl bgdataAll
