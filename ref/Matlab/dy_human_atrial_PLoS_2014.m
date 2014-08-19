function [dy] = dy_human_atrial(t, y, modelvariant)
%**************************************************************************
% Human atrial myocyte model
%
% Citation:
%
% JT Koivumäki, G Seemann, MM Maleckar, P Tavi (2014)
% In Silico Screening of the Key Cellular Remodeling Targets in Chronic
% Atrial Fibrillation.
% PLoS Comput Biol 10(5): e1003620
% http://dx.doi.org/10.1371/journal.pcbi.1003620
%
% Modified from:
%
% JT Koivumäki, T Korhonen, P Tavi (2011)
% Impact of Sarcoplasmic Reticulum Calcium Release on Calcium Dynamics and
% Action Potential Morphology in Human Atrial Myocytes: A Computational
% Study.
% PLoS Comput Biol 7(1): e1001067.
% http://dx.doi.org/10.1371/journal.pcbi.1001067
%
% Modifications made to the 2011-model:
%
% * reformulated L-type calcium current
% * SERCA pump with explicit PLB and SLN effects
% * corrected conductances of Ito and IKur (8.25 and 2.25 nS) to be
%   coherent with Maleckar et al. (2009), as slightly deviating parameter
%   values were given in two different sections of the that paper.
%
%**************************************************************************

%**********************************************
% Globals
%**********************************************

global stim INa ICaL It Isus IKr IKs IK1 INab ICab ICaP INaK INaCa If Jrelss Jrel1 Jrel2 Jrel3 J_bulkSERCA1 J_bulkSERCA2 J_bulkSERCA3 J_bulkSERCAss JSRCaleak1 JSRCaleak2 JSRCaleak3 JSRCaleakss

%**********************************************
% Definition of differential variables
%**********************************************

i_V = 1;
i_lastend = i_V;

% INa
i_start = i_lastend + 1;
i_INam = i_start; i_INah1 = i_start + 1; i_INah2 = i_start + 2;
i_lastend = i_INah2;

% ICaL 
i_start = i_lastend + 1;
i_ICaLd = i_start; i_ICaLf1 = i_start + 1; i_ICaLf2 = i_start + 2; i_ICaLfca = i_start + 3;
i_lastend = i_ICaLfca;

% It
i_start = i_lastend + 1;
i_Itr = i_start; i_Its = i_start + 1;
i_lastend = i_Its;

% Isus (Ikur)
i_start = i_lastend + 1;
i_Isusr = i_start; i_Isuss = i_start + 1;
i_lastend = i_Isuss;

% IKs
i_start = i_lastend + 1;
i_IKsn = i_start;
i_lastend = i_IKsn;

% IKr
i_start = i_lastend + 1;
i_IKrpa = i_start;
i_lastend = i_IKrpa;

% If
i_start = i_lastend + 1;
i_Ify = i_start;
i_lastend = i_Ify;

% RyR
i_start = i_lastend + 1;
i_RyRoss = i_start; i_RyRcss = i_start + 1; i_RyRass = i_start + 2;
i_RyRo1 = i_start + 3; i_RyRc1 = i_start + 4; i_RyRa1 =  i_start + 5;
i_RyRo2 = i_start + 6; i_RyRc2 = i_start + 7; i_RyRa2 =  i_start + 8;
i_RyRo3 = i_start + 9; i_RyRc3 = i_start + 10; i_RyRa3 =  i_start + 11;
i_lastend = i_RyRa3;

% SERCA
i_start = i_lastend + 1;
i_SERCACa1 = i_start; % first compartment from center
i_SERCACa2 = i_start + 1; % 2nd compartment from center
i_SERCACa3 = i_start + 2; % 3rd compartment from center
i_SERCACass = i_start + 3; % ss-SERCA
i_lastend = i_SERCACass;

% Nai ja Ki
i_start = i_lastend + 1;
i_Nass = i_start; i_Nai = i_start + 1;
i_Ki = i_start + 2;
i_lastend = i_Ki;

% Cai and CaSR 
i_start = i_lastend + 1;
i_Cass = i_start;
i_Cacenter = i_start + 1;
%i_lastend = i_Cacenter;

%**********************************************
% Stimulus
%**********************************************

% This uses previous datarow in stim-matrix for the value of stimuluscurrent at time t. In other words: interpolating to previous time-point.   

st = stim( find(stim(:,1) <= t & stim(:,1) >= 0) , 2 ); % Find all (time, current)-data which t is smaller or equal to t
if length(st) == 0
    Istim = 0;
else
    Istim = st(end); % use the last value (highest t) % for normal pacing
    %Istim = -(st(end) - y(i_V)) / 0.001; % for voltage clamp simulation, according to rate = 0.001
end

%**********************************************
% Select model variant
%**********************************************

switch modelvariant

    case 'nSR'
        cAF_lcell = 1; cAF_gCaL = 1; cAF_gt = 1; cAF_gsus = 1; cAF_gK1 = 1; cAF_kNaCa = 1; cAF_cpumps = 1; cAF_PLB = 1; cAF_SLN = 1; cAF_phos = 1; cAF_RyR = 1;
        
    case 'cAF_all'
        cAF_lcell = 1.10; cAF_gCaL = 0.41; cAF_gt = 0.38; cAF_gsus = 0.62; cAF_gK1 = 1.62; cAF_kNaCa = 1.50; cAF_cpumps = 0.84; cAF_PLB = 1.18; cAF_SLN = 0.60; cAF_phos = 2; cAF_RyR = 2;
        
    case 'cAF_dilation'
        cAF_lcell = 1.10; cAF_gCaL = 1; cAF_gt = 1; cAF_gsus = 1; cAF_gK1 = 1; cAF_kNaCa = 1; cAF_cpumps = 1; cAF_PLB = 1; cAF_SLN = 1; cAF_phos = 1; cAF_RyR = 1;
    case 'cAF_ICaL'
        cAF_lcell = 1; cAF_gCaL = 0.41; cAF_gt = 1; cAF_gsus = 1; cAF_gK1 = 1; cAF_kNaCa = 1; cAF_cpumps = 1; cAF_PLB = 1; cAF_SLN = 1; cAF_phos = 1; cAF_RyR = 1;
    case 'cAF_IK1'
        cAF_lcell = 1; cAF_gCaL = 1; cAF_gt = 1; cAF_gsus = 1; cAF_gK1 = 1.62; cAF_kNaCa = 1; cAF_cpumps = 1; cAF_PLB = 1; cAF_SLN = 1; cAF_phos = 1; cAF_RyR = 1;
    case 'cAF_Isus'
        cAF_lcell = 1; cAF_gCaL = 1; cAF_gt = 1; cAF_gsus = 0.62; cAF_gK1 = 1; cAF_kNaCa = 1; cAF_cpumps = 1; cAF_PLB = 1; cAF_SLN = 1; cAF_phos = 1; cAF_RyR = 1;
    case 'cAF_Ito'
        cAF_lcell = 1; cAF_gCaL = 1; cAF_gt = 0.38; cAF_gsus = 1; cAF_gK1 = 1; cAF_kNaCa = 1; cAF_cpumps = 1; cAF_PLB = 1; cAF_SLN = 1; cAF_phos = 1; cAF_RyR = 1;
    case 'cAF_NCX'
        cAF_lcell = 1; cAF_gCaL = 1; cAF_gt = 1; cAF_gsus = 1; cAF_gK1 = 1; cAF_kNaCa = 1.50; cAF_cpumps = 1; cAF_PLB = 1; cAF_SLN = 1; cAF_phos = 1; cAF_RyR = 1;
    case 'cAF_RyR'
        cAF_lcell = 1; cAF_gCaL = 1; cAF_gt = 1; cAF_gsus = 1; cAF_gK1 = 1; cAF_kNaCa = 1; cAF_cpumps = 1; cAF_PLB = 1; cAF_SLN = 1; cAF_phos = 1; cAF_RyR = 2;
    case 'cAF_SERCA'
        cAF_lcell = 1; cAF_gCaL = 1; cAF_gt = 1; cAF_gsus = 1; cAF_gK1 = 1; cAF_kNaCa = 1; cAF_cpumps = 0.84; cAF_PLB = 1.18; cAF_SLN = 0.60; cAF_phos = 2; cAF_RyR = 1;

    case 'cAF_no_dilation'
        cAF_lcell = 1; cAF_gCaL = 0.41; cAF_gt = 0.38; cAF_gsus = 0.62; cAF_gK1 = 1.62; cAF_kNaCa = 1.50; cAF_cpumps = 0.84; cAF_PLB = 1.18; cAF_SLN = 0.60; cAF_phos = 2; cAF_RyR = 2;
    case 'cAF_no_ICaL'
        cAF_lcell = 1.10; cAF_gCaL = 1; cAF_gt = 0.38; cAF_gsus = 0.62; cAF_gK1 = 1.62; cAF_kNaCa = 1.50; cAF_cpumps = 0.84; cAF_PLB = 1.18; cAF_SLN = 0.60; cAF_phos = 2; cAF_RyR = 2;
    case 'cAF_no_IK1'
        cAF_lcell = 1.10; cAF_gCaL = 0.41; cAF_gt = 0.38; cAF_gsus = 0.62; cAF_gK1 = 1; cAF_kNaCa = 1.50; cAF_cpumps = 0.84; cAF_PLB = 1.18; cAF_SLN = 0.60; cAF_phos = 2; cAF_RyR = 2;
    case 'cAF_no_Isus'
        cAF_lcell = 1.10; cAF_gCaL = 0.41; cAF_gt = 0.38; cAF_gsus = 1; cAF_gK1 = 1.62; cAF_kNaCa = 1.50; cAF_cpumps = 0.84; cAF_PLB = 1.18; cAF_SLN = 0.60; cAF_phos = 2; cAF_RyR = 2;
    case 'cAF_no_Ito'
        cAF_lcell = 1.10; cAF_gCaL = 0.41; cAF_gt = 1; cAF_gsus = 0.62; cAF_gK1 = 1.62; cAF_kNaCa = 1.50; cAF_cpumps = 0.84; cAF_PLB = 1.18; cAF_SLN = 0.60; cAF_phos = 2; cAF_RyR = 2;
    case 'cAF_no_NCX'
        cAF_lcell = 1.10; cAF_gCaL = 0.41; cAF_gt = 0.38; cAF_gsus = 0.62; cAF_gK1 = 1.62; cAF_kNaCa = 1; cAF_cpumps = 0.84; cAF_PLB = 1.18; cAF_SLN = 0.60; cAF_phos = 2; cAF_RyR = 2;
    case 'cAF_no_RyR'
        cAF_lcell = 1.10; cAF_gCaL = 0.41; cAF_gt = 0.38; cAF_gsus = 0.62; cAF_gK1 = 1.62; cAF_kNaCa = 1.50; cAF_cpumps = 0.84; cAF_PLB = 1.18; cAF_SLN = 0.60; cAF_phos = 2; cAF_RyR = 1;
    case 'cAF_no_SERCA'
        cAF_lcell = 1.10; cAF_gCaL = 0.41; cAF_gt = 0.38; cAF_gsus = 0.62; cAF_gK1 = 1.62; cAF_kNaCa = 1.50; cAF_cpumps = 1; cAF_PLB = 1; cAF_SLN = 1; cAF_phos = 1; cAF_RyR = 2;
        
    otherwise
        disp('Unknown model variant.')
end

%**********************************************
% Numerical parameters
%**********************************************

% Physical & environmental constants
F = 96487;
R = 8314;
T = 306.15; % 33C

Nao = 130;
Nass = y(i_Nass);

Cao = 1.8;

Ko = 5.4;
Ki = y(i_Ki);

Cm = 0.05; % nF

% Cell dilation in cAF
Ddcell = (cAF_lcell - 1)*(20/10) + 1;
Dvcell = cAF_lcell*Ddcell^2;

% Geometry
Vss = 4.99232e-5 * Dvcell; % nL
rjunct = 6.5 * Ddcell; % mum
lcell = 122.051 * cAF_lcell; % mum

% Ca diffusion grid
dx = 1.625 * Ddcell; % mum
rstart = 0 + 0.5*dx;
rend = rjunct - 0.5*dx;
j = round(rstart/dx):1:round(rjunct/dx); % Spatial index of Cai diffusion
j = j';

Aj_nj = pi*rjunct*2*lcell*0.5; % Area between junct and nonjunct
xj_nj = 0.02/2 * Ddcell + dx/2; % diffusion distance from center to center of junct to first njunct
xj_nj_Nai = 0.02/2 * Ddcell + 2*dx; % diffusion distance from center of junct to center of njunct (between 2nd and 3rd njunct)

% Diffusion compartment volumes
Vnonjunct = zeros(length(j),1);
Vnonjunct = (pi.*(j.*dx).^2.*lcell-pi.*((j-1).*dx).^2.*lcell).*1e-6.*0.5; %nL

Vcytosol = sum(Vnonjunct) + Vss;

VSR = 0.05.*Vnonjunct./2*0.9 / Dvcell;

Vnonjunct_Nai = sum(Vnonjunct);

% Non-junct Cai data & CaSR data
Cai = zeros(length(j),1);
Cai = y(i_Cacenter:length(j)+i_Cacenter-1);
CaSR = y(length(j)+i_Cacenter:length(j)*2+i_Cacenter-1);

% Cytosol Ca Buffers
Begta = 0;
%Begta = 10; % 10 mM EGTA for voltage clamp measurements
Bcmdn = 24e-3;
BCa = Bcmdn + Begta;
SLlow = 165;
SLhigh = 13;

KdBegta = 0.12e-3;
KdBcmdn = 2.38e-3;
KdBCa = (Bcmdn*KdBcmdn + Begta*KdBegta) / (Bcmdn + Begta);
KdSLlow = 1.1;
KdSLhigh = 13e-3;

% SR Ca buffers
CSQN =  6.7;
KdCSQN = 0.8;

% Sarcolemmal Na burrering
BNa = 0.49 * 2.31;
KdBNa = 10;


% Ion channel conductances & permeabilities & other parameters
PNa = 0.00182;

ECa_app = 60; kCa = 0.6e-3; gCaL = 15 * cAF_gCaL;

gt = 8.25 * cAF_gt;
gsus = 2.25 * cAF_gsus;
gKs = 1;
gKr = 0.50 * sqrt(Ko/5.4);
gK1 = 3.5 * cAF_gK1;

gNab = 0.060599;
gCab = 0.0952;

INaKmax = 70.8253; kNaKK = 1; kNaKNa = 11;

ICaPmax = 2.0; kCaP = 0.0005;

kNaCa = 0.0084 * cAF_kNaCa; gam = 0.45; dNaCa = 0.0003;

gIf = 1;

% Ca and Na diffusion
DCa = 780;
DCaSR = 44;
DCaBm = 25;
DNa = 0.12;

% SERCA
base_phos = 0.1 * cAF_phos; % Baseline phosphorylation
PLB_SERCA_ratio = cAF_PLB;
SLN_SERCA_ratio = cAF_SLN;
Kmf_PLBKO = 0.15e-3; Kmf_PLB = 0.12e-3; Kmf_SLN = 0.07e-3;
Kmr_PLBKO = 2.5; Kmr_PLB = 0.88; Kmr_SLN = 0.5;
SERCAKmf = Kmf_PLBKO + Kmf_PLB * PLB_SERCA_ratio * (1 - base_phos) + Kmf_SLN * SLN_SERCA_ratio * (1 - base_phos);
SERCAKmr = Kmr_PLBKO - Kmr_PLB * PLB_SERCA_ratio * (1 - base_phos) - Kmr_SLN * SLN_SERCA_ratio * (1 - base_phos);
k4 = 13; % pump rate
k3 = k4 / SERCAKmr^2;
k1 = 1000^2 * k4; 
k2 = k1 * SERCAKmf^2;
cpumps = 30e-3 / Dvcell * cAF_cpumps; % pump concentration in cytosol volume, in cAF SR does not dilate

% RyR
RyRtauadapt = 1;
RyRtauactss = 5e-3; RyRtauinactss = 15e-3;
RyRtauact = 18.75e-3; RyRtauinact = 87.5e-3;

% SR Ca leak
kSRleak = 6e-3; 


%**********************************************
% Analytical equations
%**********************************************

% Reversal potentials
ENa = R*T/F * log ( Nao / Nass );
EK = R*T/F * log ( Ko / Ki );
ECa = R*T/F/2 * log ( Cao / y(i_Cass) );

% INa
INa = PNa .* y(i_INam).^3 .* ( 0.9.*y(i_INah1) + 0.1.*y(i_INah2) ) * Nao * y(i_V) * F^2/(R*T) * ( exp( (y(i_V)-ENa)*F/R/T ) - 1) / ( exp( y(i_V)*F/R/T ) - 1);
INaminf = 1/(1+exp((y(i_V)+27.12)/-8.21));
INahinf = 1/(1+exp((y(i_V)+63.6)/5.3));
INamtau = 0.000042*exp( -((y(i_V)+25.57)/28.8).^2 ) + 0.000024;
INah1tau = 0.03/(1+exp((y(i_V)+35.1)/3.2)) + 0.0003;
INah2tau = 0.12/(1+exp((y(i_V)+35.1)/3.2)) + 0.003;

% ICaL
ICaLfcainf = 1 / ( 1 + (y(i_Cass)/kCa)^2);
ICaL = gCaL * y(i_ICaLd) * y(i_ICaLf2) * y(i_ICaLfca) * (y(i_V) - ECa_app);
ICaLdinf = 1/(1+exp((y(i_V)+9)/-5.8));
ICaLfinf = 0.04 + 0.96 / (1 + exp((y(i_V) + 25.5)/8.4)) + 1 / (1 + exp(-(y(i_V) - 60)/8.0)); % fit Li et al. data
ICaLdtau = 0.00065 * exp(-((y(i_V) + 35)/30).^2) + 0.0005; % Nygren four-fold faster, according to Cavalie et al.
ICaLf1tau = 0.98698.*exp( -((y(i_V)+30.16047)./7.09396).^2 ) + 0.04275./(1+exp((y(i_V)-51.61555)/-80.61331)) + 0.03576./(1+exp((y(i_V)+29.57272)/13.21758)) - 0.00821; % not used in the current equation
ICaLf2tau = 1.34*exp( -((y(i_V)+40)/14.2).^2 ) + 0.04; % average from Li et al. and Christ et al. (55 + 51.8 + 20.6 + 33.9)/4
ICaLfcatau = 2e-3;


% It
It = gt * y(i_Itr) * y(i_Its) * (y(i_V) - EK);
Itrinf = 1/(1+exp((y(i_V)-1)/-11));
Itsinf = 1/(1+exp((y(i_V)+40.5)/11.5));
Itrtau = 0.0035*exp( -((y(i_V)+0)/30).^2 ) + 0.0015;
Itstau = 0.025635*exp( -((y(i_V)+52.45)/15.8827).^2 ) + 0.01414; % Maleckar et al.

% Isus
Isus = gsus * y(i_Isusr) * y(i_Isuss) * (y(i_V) - EK);
Isusrinf = 1/(1 + exp((y(i_V) + 6)/-8.6)); % Maleckar et al.
Isussinf = 1/(1 + exp((y(i_V) + 7.5)/10)); % Maleckar et al.
Isusrtau = 0.009/(1 + exp((y(i_V) + 5)/12)) + 0.0005; % Maleckar et al.
Isusstau = 0.59/(1 + exp((y(i_V) + 60)/10)) + 3.05; % Maleckar et al.

% IKr
IKrpi = 1/(1+exp((y(i_V)+55)/24));
IKr = gKr * y(i_IKrpa) * IKrpi * (y(i_V) - EK);
IKrpainf = 1/(1+exp((y(i_V)+15)/-6));
IKrpatau = 0.21718*exp( -((y(i_V)+20.1376)/22.1996).^2 ) + 0.03118;

% IKs
IKs = gKs * y(i_IKsn) * (y(i_V) - EK);
IKsninf = 1/(1+exp((y(i_V)-19.9)/-12.7));
IKsntau = 0.4*exp( -((y(i_V)-20)/20).^2 ) + 0.7;

% IK1
IK1 = gK1 * Ko.^0.4457 * (y(i_V) - EK) / (1 + exp(1.5*(y(i_V)-EK+3.6)*F/R/T));

% Background leaks
INab = gNab * (y(i_V) - ENa);
ICab = gCab * (y(i_V) - ECa);

% INaK
INaK = INaKmax * Ko./(Ko + kNaKK) * Nass.^1.5/(Nass.^1.5 + kNaKNa.^1.5) * (y(i_V) + 150) / (y(i_V) + 200);

% INaCa
fCaNCX = 1;
INaCa = kNaCa * ( (exp( gam.*y(i_V)*F/R/T ) .* Nass.^3 .* Cao - exp( (gam-1).*y(i_V)*F/R/T ) .* Nao^3 .* y(i_Cass)*fCaNCX ) / ( 1 + dNaCa*(Nao^3 .* y(i_Cass)*fCaNCX + Nass.^3 .* Cao) ) );

% ICaP
ICaP = ICaPmax * y(i_Cass) / (kCaP + y(i_Cass));

% If, Zorn-Pauly LAW fit
Ifyinf = 1 / (1 + exp((y(i_V)+97.82874)/12.48025));
Ifytau = 1 ./ (0.00332.*exp(-y(i_V)./16.54103)+23.71839.*exp(y(i_V)./16.54103));
IfNa = gIf * y(i_Ify)*((0.2677)*(y(i_V)-ENa));
IfK = gIf * y(i_Ify)*((1-0.2677)*(y(i_V)-EK)); 
If = IfK + IfNa;

% Ca buffers
betass = ( 1 + SLlow*KdSLlow./(y(i_Cass) + KdSLlow).^2 + SLhigh*KdSLhigh./(y(i_Cass) + KdSLhigh).^2 + BCa*KdBCa./(y(i_Cass) + KdBCa).^2  ).^(-1);
betai = ( 1 + BCa.*KdBCa./(Cai + KdBCa).^2  ).^(-1);
gammai = BCa.*KdBCa./(Cai + KdBCa).^2;

betaSR = ( 1 + CSQN.*KdCSQN./(CaSR + KdCSQN).^2 ).^(-1); 

betaNass = ( 1 + BNa*KdBNa./(y(i_Nass) + KdBNa).^2 ).^(-1);

%  Diffusion from junct to non-junct
Jj_nj = DCa * Aj_nj / xj_nj * (y(i_Cass)-Cai(end)).*1e-6;

% SERCA fluxes 
J_SERCASR1 = (-k3*CaSR(1).^2*(cpumps-y(i_SERCACa1))+k4*y(i_SERCACa1))*Vnonjunct(1)*2; % in 1 nl volume
J_bulkSERCA1 = (k1*Cai(1).^2*(cpumps-y(i_SERCACa1))-k2*y(i_SERCACa1))*Vnonjunct(1)*2; % in 1 nl volume

J_SERCASR2 = (-k3*CaSR(2).^2*(cpumps-y(i_SERCACa2))+k4*y(i_SERCACa2))*Vnonjunct(2)*2; % in 1 nl volume
J_bulkSERCA2 = (k1*Cai(2).^2*(cpumps-y(i_SERCACa2))-k2*y(i_SERCACa2))*Vnonjunct(2)*2; % in 1 nl volume

J_SERCASR3 = (-k3*CaSR(3).^2*(cpumps-y(i_SERCACa3))+k4*y(i_SERCACa3))*Vnonjunct(3)*2; % in 1 nl volume
J_bulkSERCA3 = (k1*Cai(3).^2*(cpumps-y(i_SERCACa3))-k2*y(i_SERCACa3))*Vnonjunct(3)*2; % in 1 nl volume

J_SERCASRss = (-k3*CaSR(4).^2*(cpumps-y(i_SERCACass))+k4*y(i_SERCACass))*Vss*2; % in 1 nl volume
J_bulkSERCAss = (k1*y(i_Cass).^2*(cpumps-y(i_SERCACass))-k2*y(i_SERCACass))*Vss*2; % in 1 nl volume

% RyR
nuss = 625*Vss / Dvcell; % in cAF SR does not dilate;
RyRSRCass = (1 - 1./(1 +  exp((CaSR(4)-0.3/cAF_RyR)./0.1)));
RyRainfss = 0.505-0.427./(1 + exp((y(i_Cass).*1000-0.29)./0.082));
RyRoinfss = (1 - 1./(1 +  exp((y(i_Cass).*1000-(y(i_RyRass) + 0.22/cAF_RyR))./0.03)));
RyRcinfss = (1./(1 + exp((y(i_Cass).*1000-(y(i_RyRass)+0.02))./0.01)));
Jrelss = nuss * ( y(i_RyRoss) ) * y(i_RyRcss) * RyRSRCass * ( CaSR(4) -  y(i_Cass) ); 

nu1 = 1*Vnonjunct(1) / Dvcell; % in cAF SR does not dilate
RyRSRCa1 = (1 - 1./(1 +  exp((CaSR(1)-0.3/cAF_RyR)./0.1)));
RyRainf1 = 0.505-0.427./(1 + exp((Cai(1).*1000-0.29)./0.082));
RyRoinf1 = (1 - 1./(1 +  exp(( Cai(1).*1000-(y(i_RyRa1) + 0.22/cAF_RyR))./0.03)));
RyRcinf1 = (1./(1 +  exp(( Cai(1).*1000-(y(i_RyRa1)+0.02))./0.01)));
Jrel1 = nu1 * ( y(i_RyRo1) ) * y(i_RyRc1) * RyRSRCa1 * ( CaSR(1) -  Cai(1) ); 

nu2 = 1*Vnonjunct(2) / Dvcell; % in cAF SR does not dilate
RyRSRCa2 = (1 - 1./(1 +  exp((CaSR(2)-0.3/cAF_RyR)./0.1)));
RyRainf2 =  0.505-0.427./(1 + exp((Cai(2).*1000-0.29)./0.082));
RyRoinf2 = (1 - 1./(1 +  exp(( Cai(2).*1000-(y(i_RyRa2) + 0.22/cAF_RyR))./0.03)));
RyRcinf2 = (1./(1 +  exp(( Cai(2).*1000-(y(i_RyRa2)+0.02))./0.01)));
Jrel2 = nu2 * ( y(i_RyRo2) ) * y(i_RyRc2) * RyRSRCa2 * ( CaSR(2) -  Cai(2) ); 

nu3 = 1*Vnonjunct(3) / Dvcell; % in cAF SR does not dilate
RyRSRCa3 = (1 - 1./(1 +  exp((CaSR(3)-0.3/cAF_RyR)./0.1)));
RyRainf3 =  0.505-0.427./(1 + exp((Cai(3).*1000-0.29)./0.082));
RyRoinf3 = (1 - 1./(1 +  exp(( Cai(3).*1000-(y(i_RyRa3) + 0.22/cAF_RyR))./0.03)));
RyRcinf3 = (1./(1 +  exp(( Cai(3).*1000-(y(i_RyRa3)+0.02))./0.01)));
Jrel3 = nu3 * ( y(i_RyRo3) ) * y(i_RyRc3) * RyRSRCa3 * ( CaSR(3) -  Cai(3) ); 

% SR leak fluxes
JSRCaleak1 = kSRleak * ( CaSR(1) - Cai(1) ) * Vnonjunct(1) / Dvcell; % in cAF SR does not dilate
JSRCaleak2 = kSRleak * ( CaSR(2) - Cai(2) ) * Vnonjunct(2) / Dvcell; % in cAF SR does not dilate
JSRCaleak3 = kSRleak * ( CaSR(3) - Cai(3) ) * Vnonjunct(3) / Dvcell; % in cAF SR does not dilate
JSRCaleakss = kSRleak * ( CaSR(4) - y(i_Cass) ) * Vss / Dvcell; % in cAF SR does not dilate

% Cafluxes in 1 nl volume
JCa = zeros(length(j),1);
JCa(1) = -J_bulkSERCA1 + JSRCaleak1 + Jrel1;
JCa(2) = -J_bulkSERCA2 + JSRCaleak2 + Jrel2;
JCa(3) = -J_bulkSERCA3 + JSRCaleak3 + Jrel3;
JCa(4) = Jj_nj;
JCass = -Jj_nj + JSRCaleakss - J_bulkSERCAss + Jrelss;

JSRCa = zeros(length(j),1);
JSRCa(1) = J_SERCASR1 - JSRCaleak1 - Jrel1;
JSRCa(2) = J_SERCASR2 - JSRCaleak2 - Jrel2;
JSRCa(3) = J_SERCASR3 - JSRCaleak3 - Jrel3;
JSRCa(4) = J_SERCASRss - JSRCaleakss - Jrelss;

% Naflux in 1 nl volume
JNa = DNa * Aj_nj / xj_nj_Nai * (y(i_Nass) - y(i_Nai))* 1e-6;


%**********************************************
% Differential equations
%**********************************************

dy = zeros(length(j)*2+i_Cacenter-1,1);

% V
dy(i_V) = (INa + ICaL + It + Isus + IK1 + IKr + IKs + INab + ICab + INaK + ICaP + INaCa + If + Istim)/(-Cm); % currents are in (pA)

% INa
dy(i_INam) = (INaminf - y(i_INam))/INamtau;
dy(i_INah1) = (INahinf - y(i_INah1))/INah1tau;
dy(i_INah2) = (INahinf - y(i_INah2))/INah2tau;

% ICaL
dy(i_ICaLd) = (ICaLdinf - y(i_ICaLd))/ICaLdtau;
dy(i_ICaLf1) = (ICaLfinf - y(i_ICaLf1))/ICaLf1tau;
dy(i_ICaLf2) = (ICaLfinf - y(i_ICaLf2))/ICaLf2tau;
dy(i_ICaLfca) = (ICaLfcainf - y(i_ICaLfca))/ICaLfcatau;

% It
dy(i_Itr) = (Itrinf - y(i_Itr))/Itrtau;
dy(i_Its) = (Itsinf - y(i_Its))/Itstau;

% Isus
dy(i_Isusr) = (Isusrinf - y(i_Isusr))/Isusrtau;
dy(i_Isuss) = (Isussinf - y(i_Isuss))/Isusstau;

% IKs
dy(i_IKsn) = (IKsninf - y(i_IKsn))/IKsntau;

% IKr
dy(i_IKrpa) = (IKrpainf - y(i_IKrpa))/IKrpatau;

% If
dy(i_Ify) = (Ifyinf - y(i_Ify))/Ifytau;

% SERCACa
dy(i_SERCACa1) = 0.5*(-J_SERCASR1 + J_bulkSERCA1)/Vnonjunct(1); 
dy(i_SERCACa2) = 0.5*(-J_SERCASR2 + J_bulkSERCA2)/Vnonjunct(2); 
dy(i_SERCACa3) = 0.5*(-J_SERCASR3 + J_bulkSERCA3)/Vnonjunct(3); 
dy(i_SERCACass) = 0.5*(-J_SERCASRss + J_bulkSERCAss)/Vss; 

% RyR
dy(i_RyRoss) = (RyRoinfss-y(i_RyRoss))./RyRtauactss;
dy(i_RyRcss) = (RyRcinfss-y(i_RyRcss))./RyRtauinactss;
dy(i_RyRass) = (RyRainfss-y(i_RyRass))./RyRtauadapt;
dy(i_RyRo1) = (RyRoinf1-y(i_RyRo1))./RyRtauact;
dy(i_RyRc1) = (RyRcinf1-y(i_RyRc1))./RyRtauinact;
dy(i_RyRa1) = (RyRainf1-y(i_RyRa1))./RyRtauadapt;
dy(i_RyRo2) = (RyRoinf2-y(i_RyRo2))./RyRtauact;
dy(i_RyRc2) = (RyRcinf2-y(i_RyRc2))./RyRtauinact;
dy(i_RyRa2) = (RyRainf2-y(i_RyRa2))./RyRtauadapt;
dy(i_RyRo3) = (RyRoinf3-y(i_RyRo3))./RyRtauact;
dy(i_RyRc3) = (RyRcinf3-y(i_RyRc3))./RyRtauinact;
dy(i_RyRa3) = (RyRainf3-y(i_RyRa3))./RyRtauadapt;

% Nai & Ki
dy(i_Nass) = betaNass * (-JNa/Vss -(INa + INab + 3*INaK + 3*INaCa + IfNa) / (Vss*F));
dy(i_Nai) = JNa/Vnonjunct_Nai;
dy(i_Ki) = -(It + Isus + IK1 + IKr + IKs - 2*INaK + IfK + Istim) / (Vcytosol*F);

% Ca
dy(i_Cass) = betass * ( JCass/Vss + (-ICaL - ICab - ICaP + 2*INaCa) / (2*Vss*F) );

dCaidt = zeros(length(Cai),1); 
dCaidt(1) = betai(1) .* (DCa + gammai(1).*DCaBm) .* ( (Cai(2)-2.*Cai(1)+Cai(1))./dx.^2 + (Cai(2)-Cai(1))./(2.*j(1).*dx.^2) ) - 2.*betai(1).*gammai(1).*DCaBm./(KdBCa + Cai(1)) .* ((Cai(2)-Cai(1))./(2.*dx)).^2 + JCa(1)./Vnonjunct(1).*betai(1); 
dCaidt(2:end-1) = betai(2:end-1) .* (DCa + gammai(2:end-1).*DCaBm) .* ( (Cai(3:end)-2.*Cai(2:end-1)+Cai(1:end-2))./dx.^2 + (Cai(3:end)-Cai(1:end-2))./(2.*j(2:end-1).*dx.^2) ) - 2.*betai(2:end-1).*gammai(2:end-1).*DCaBm./(KdBCa + Cai(2:end-1)) .* ((Cai(3:end)-Cai(1:end-2))./(2.*dx)).^2 + JCa(2:end-1)./Vnonjunct(2:end-1).*betai(2:end-1); 
dCaidt(end) = betai(end) .* (DCa + gammai(end).*DCaBm) .* ( (Cai(end)-2.*Cai(end)+Cai(end-1))./dx.^2 + (Cai(end)-Cai(end-1))./(2.*j(end).*dx.^2) ) - 2.*betai(end).*gammai(end).*DCaBm./(KdBCa + Cai(end)) .* ((Cai(end)-Cai(end-1))./(2.*dx)).^2 + JCa(end)./Vnonjunct(end).*betai(end); 

dy(i_Cacenter:length(Cai)+i_Cacenter-1) = dCaidt;

dCaSRdt = zeros(length(CaSR),1); 
dCaSRdt(1) = betaSR(1) .* (DCaSR) .* ( (CaSR(2)-2.*CaSR(1)+CaSR(1))./dx.^2 + (CaSR(2)-CaSR(1))./(2.*j(1).*dx.^2) ) + JSRCa(1)./VSR(1).*betaSR(1); 
dCaSRdt(2:end-1) = betaSR(2:end-1) .* (DCaSR) .* ( (CaSR(3:end)-2.*CaSR(2:end-1)+CaSR(1:end-2))./dx.^2 + (CaSR(3:end)-CaSR(1:end-2))./(2.*j(2:end-1).*dx.^2) ) + JSRCa(2:end-1)./VSR(2:end-1).*betaSR(2:end-1); 
dCaSRdt(end) = betaSR(end) .* (DCaSR) .* ( (CaSR(end)-2.*CaSR(end)+CaSR(end-1))./dx.^2 + (CaSR(end)-CaSR(end-1))./(2.*j(end).*dx.^2) ) + JSRCa(end)./VSR(end).*betaSR(end); 

dy(i_Cacenter+length(j):length(j).*2+i_Cacenter-1) = dCaSRdt;
