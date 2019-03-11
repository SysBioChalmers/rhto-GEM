%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [model,GAMpol] = changeOtherComp(model,expData,i)
%
%   Modified from SLIMEr, under MIT License:
%   https://github.com/SysBioChalmers/SLIMEr/blob/master/models/scaleAbundancesInModel.m
%
% Benjamin J. Sanchez. Last update: 2018-09-04
% Eduard Kerkhoven. Last modified: 2019-02-19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,GAMpol] = changeOtherComp(model,expData,i)

if nargin<3
    i=1;
end

%Components of biomass:
%        id             MW [g/mol]  class     name
comps = {'s_0404'	89.09       'P'     % A     Alanine         ala
         's_0542'    121.16      'P'     % C     Cysteine        cys
         's_0432'    133.11      'P'     % D     Aspartic acid   asp
         's_0748'    147.13      'P'     % E     Glutamic acid   glu
         's_1314'    165.19      'P'     % F     Phenylalanine   phe
         's_0757'    75.07       'P'     % G     Glycine         gly
         's_0832'    155.15      'P'     % H     Histidine       his
         's_0847'    131.17      'P'     % I     Isoleucine      ile
         's_1099'    146.19      'P'     % K     Lysine          lys
         's_1077'    131.17      'P'     % L     Leucine         leu
         's_1148'    149.21      'P'     % M     Methionine      met
         's_0430'    132.12      'P'     % N     Asparagine      asn
         's_1379'    115.13      'P'     % P     Proline         pro
         's_0747'    146.14      'P'     % Q     Glutamine       gln
         's_0428'    174.2       'P'     % R     Arginine        arg
         's_1428'    105.09      'P'     % S     Serine          ser
         's_1491'    119.12      'P'     % T     Threonine       thr
         's_1561'    117.15      'P'     % V     Valine          val
         's_1527'    204.23      'P'     % W     Tryptophan      trp
         's_1533'    181.19      'P'     % Y     Tyrosine        tyr
         's_0001'	 180.16      'C'     % (1->3)-beta-D-glucan
         's_0004'	 180.16      'C'     % (1->6)-beta-D-glucan
         's_0773'    180.16      'C'     % glycogen
         's_1107'    180.16      'C'     % mannan
         's_1520'    342.296 	 'C'     % trehalose
         's_0423'    347.22      'R'     % AMP
         's_0526'    323.2       'R'     % CMP
         's_0782'    363.22      'R'     % GMP
         's_1545'    324.18      'R'     % UMP
         's_0584'    331.22      'D'     % dAMP
         's_0589'    307.2       'D'     % dCMP
         's_0615'    345.21      'D'     % dGMP
         's_0649'    322.21      'D'     % dTMP
         's_3714'    852.83      'N'     % heme a
         's_1405'    376.36      'N'     % riboflavin
         's_1467'    96.06       'N'};   % sulphate

%Change given abundances in model:
[~,P,~,~,~,~]   = sumBioMass(model);
protPos         = strcmp(model.rxnNames,'protein pseudoreaction');
fP              = expData.protein(i)/P;       	%ratio to scale
isAA            = contains(model.metNames,'tRNA');	%protein components
model.S(isAA,protPos) = full(model.S(isAA,protPos))*fP;

%Compute new biomass and lipid fraction:
[X,~,C,~,~,~] = sumBioMass(model);
delta         = X - 1;                      %difference to balance

%Balance out mass with all sugars:
mets     = comps(strcmp(comps(:,3),'C'),1);
massPre  = C;
massPost = massPre - delta;
fC       = massPost/massPre;
carbPos  = find(strcmp(model.rxnNames,'carbohydrate pseudoreaction'));
for i = 1:length(mets)
    modelPos = strcmp(model.mets,mets{i});
    model.S(modelPos,carbPos) = model.S(modelPos,carbPos)*fC;
end

%Estimate maintenance belonging to polymerization:
[X,P,C,R,D,~] = sumBioMass(model);
GAMpol        = P*37.7 + C*12.8 + R*26.0 + D*26.0;    %Forster 2003 (sup table 8)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%