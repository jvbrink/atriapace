#include <cmath>
#include <cstring>
#include <stdexcept>
// Gotran generated C/C++ code for the "hAM_KSMT_nSR" model

// Init state values
void init_state_values(double* states)
{
  states[0] = -75.2; // V;
  states[1] = 0.002854; // INam;
  states[2] = 0.9001; // INah1;
  states[3] = 0.9003; // INah2;
  states[4] = 1.104e-05; // ICaLd;
  states[5] = 0.9988; // ICaLf1;
  states[6] = 0.9988; // ICaLf2;
  states[7] = 0.9226; // ICaLfca;
  states[8] = 0.0009795; // Itr;
  states[9] = 0.9535; // Its;
  states[10] = 0.0003196; // Isusr;
  states[11] = 0.9726; // Isuss;
  states[12] = 0.004421; // IKsn;
  states[13] = 4.354e-05; // IKrpa;
  states[14] = 0.055; // Ify;
  states[15] = 5.934e-05; // RyRoss;
  states[16] = 0.9999; // RyRcss;
  states[17] = 0.2453; // RyRass;
  states[18] = 0.0001916; // RyRo1;
  states[19] = 0.997; // RyRc1;
  states[20] = 0.2033; // RyRa1;
  states[21] = 0.0001485; // RyRo2;
  states[22] = 0.9985; // RyRc2;
  states[23] = 0.2106; // RyRa2;
  states[24] = 0.000101; // RyRo3;
  states[25] = 0.9993; // RyRc3;
  states[26] = 0.2229; // RyRa3;
  states[27] = 0.005351; // SERCACa1;
  states[28] = 0.005213; // SERCACa2;
  states[29] = 0.00502; // SERCACa3;
  states[30] = 0.004891; // SERCACass;
  states[31] = 8.513; // Nass;
  states[32] = 8.851; // Nai;
  states[33] = 135.2; // Ki;
  states[34] = 0.0001737; // Cass;
  states[35] = 0.0001672; // Cai1;
  states[36] = 0.000167; // Cai2;
  states[37] = 0.000168; // Cai3;
  states[38] = 0.0001719; // Cai4;
  states[39] = 0.5755; // CaSR1;
  states[40] = 0.5664; // CaSR2;
  states[41] = 0.5535; // CaSR3;
  states[42] = 0.5428; // CaSR4;
}

// Default parameter values
void init_parameters_values(double* parameters)
{
  parameters[0] = -2820; // stim_amplitude;
  parameters[1] = 1; // stim_duration;
  parameters[2] = 0; // stim_offset;
  parameters[3] = 1000; // stim_period;
  parameters[4] = 130; // Nao;
  parameters[5] = 1.8; // Cao;
  parameters[6] = 5.4; // Ko;
  parameters[7] = 0.05; // Cm;
  parameters[8] = 0.024; // BCa;
  parameters[9] = 165; // SLlow;
  parameters[10] = 13; // SLhigh;
  parameters[11] = 0.00238; // KdBCa;
  parameters[12] = 1.1; // KdSLlow;
  parameters[13] = 0.013; // KdSLhigh;
  parameters[14] = 6.7; // CSQN;
  parameters[15] = 0.8; // KdCSQN;
  parameters[16] = 1.1319; // BNa;
  parameters[17] = 10; // KdBNa;
  parameters[18] = 0.0018; // PNa;
  parameters[19] = 60; // ECa_app;
  parameters[20] = 2; // kCan;
  parameters[21] = 0.0006; // kCa;
  parameters[22] = 1; // gKs;
  parameters[23] = 3.45; // gK1;
  parameters[24] = 0.060599; // gNab;
  parameters[25] = 0.0952; // gCab;
  parameters[26] = 70.8253; // INaKmax;
  parameters[27] = 1; // kNaKK;
  parameters[28] = 11; // kNaKNa;
  parameters[29] = 2.0; // ICaPmax;
  parameters[30] = 0.0005; // kCaP;
  parameters[31] = 0.45; // gam;
  parameters[32] = 0.0003; // dNaCa;
  parameters[33] = 1; // gIf;
  parameters[34] = 780; // DCa;
  parameters[35] = 44; // DCaSR;
  parameters[36] = 25; // DCaBm;
  parameters[37] = 0.12; // DNa;
  parameters[38] = 0.00015; // Kmf_PLBKO;
  parameters[39] = 0.00012; // Kmf_PLB;
  parameters[40] = 7e-05; // Kmf_SLN;
  parameters[41] = 2.5; // Kmr_PLBKO;
  parameters[42] = 0.88; // Kmr_PLB;
  parameters[43] = 0.5; // Kmr_SLN;
  parameters[44] = 13; // k4;
  parameters[45] = 0.006; // kSRleak;
}

// State index
int state_index(const char name[])
{
  // State names
  char names[][10] = {"V", "INam", "INah1", "INah2", "ICaLd", "ICaLf1",
    "ICaLf2", "ICaLfca", "Itr", "Its", "Isusr", "Isuss", "IKsn", "IKrpa",
    "Ify", "RyRoss", "RyRcss", "RyRass", "RyRo1", "RyRc1", "RyRa1", "RyRo2",
    "RyRc2", "RyRa2", "RyRo3", "RyRc3", "RyRa3", "SERCACa1", "SERCACa2",
    "SERCACa3", "SERCACass", "Nass", "Nai", "Ki", "Cass", "Cai1", "Cai2",
    "Cai3", "Cai4", "CaSR1", "CaSR2", "CaSR3", "CaSR4"};

  int i;
  for (i=0; i<43; i++)
  {
    if (strcmp(names[i], name)==0)
    {
      return i;
    }
  }
  return -1;
}

// Parameter index
int parameter_index(const char name[])
{
  // Parameter names
  char names[][15] = {"stim_amplitude", "stim_duration", "stim_offset",
    "stim_period", "Nao", "Cao", "Ko", "Cm", "BCa", "SLlow", "SLhigh",
    "KdBCa", "KdSLlow", "KdSLhigh", "CSQN", "KdCSQN", "BNa", "KdBNa", "PNa",
    "ECa_app", "kCan", "kCa", "gKs", "gK1", "gNab", "gCab", "INaKmax",
    "kNaKK", "kNaKNa", "ICaPmax", "kCaP", "gam", "dNaCa", "gIf", "DCa",
    "DCaSR", "DCaBm", "DNa", "Kmf_PLBKO", "Kmf_PLB", "Kmf_SLN", "Kmr_PLBKO",
    "Kmr_PLB", "Kmr_SLN", "k4", "kSRleak"};

  int i;
  for (i=0; i<46; i++)
  {
    if (strcmp(names[i], name)==0)
    {
      return i;
    }
  }
  return -1;
}

// Compute the right hand side of the hAM_KSMT_nSR ODE
void rhs(const double* states, const double t, const double* parameters,
  double* values)
{

  // Assign states
  const double V = states[0];
  const double INam = states[1];
  const double INah1 = states[2];
  const double INah2 = states[3];
  const double ICaLd = states[4];
  const double ICaLf1 = states[5];
  const double ICaLf2 = states[6];
  const double ICaLfca = states[7];
  const double Itr = states[8];
  const double Its = states[9];
  const double Isusr = states[10];
  const double Isuss = states[11];
  const double IKsn = states[12];
  const double IKrpa = states[13];
  const double Ify = states[14];
  const double RyRoss = states[15];
  const double RyRcss = states[16];
  const double RyRass = states[17];
  const double RyRo1 = states[18];
  const double RyRc1 = states[19];
  const double RyRa1 = states[20];
  const double RyRo2 = states[21];
  const double RyRc2 = states[22];
  const double RyRa2 = states[23];
  const double RyRo3 = states[24];
  const double RyRc3 = states[25];
  const double RyRa3 = states[26];
  const double SERCACa1 = states[27];
  const double SERCACa2 = states[28];
  const double SERCACa3 = states[29];
  const double SERCACass = states[30];
  const double Nass = states[31];
  const double Nai = states[32];
  const double Ki = states[33];
  const double Cass = states[34];
  const double Cai1 = states[35];
  const double Cai2 = states[36];
  const double Cai3 = states[37];
  const double Cai4 = states[38];
  const double CaSR1 = states[39];
  const double CaSR2 = states[40];
  const double CaSR3 = states[41];
  const double CaSR4 = states[42];

  // Assign parameters
  const double stim_amplitude = parameters[0];
  const double stim_duration = parameters[1];
  const double stim_offset = parameters[2];
  const double stim_period = parameters[3];
  const double Nao = parameters[4];
  const double Cao = parameters[5];
  const double Ko = parameters[6];
  const double Cm = parameters[7];
  const double BCa = parameters[8];
  const double SLlow = parameters[9];
  const double SLhigh = parameters[10];
  const double KdBCa = parameters[11];
  const double KdSLlow = parameters[12];
  const double KdSLhigh = parameters[13];
  const double CSQN = parameters[14];
  const double KdCSQN = parameters[15];
  const double BNa = parameters[16];
  const double KdBNa = parameters[17];
  const double PNa = parameters[18];
  const double ECa_app = parameters[19];
  const double kCan = parameters[20];
  const double kCa = parameters[21];
  const double gKs = parameters[22];
  const double gK1 = parameters[23];
  const double gNab = parameters[24];
  const double gCab = parameters[25];
  const double INaKmax = parameters[26];
  const double kNaKK = parameters[27];
  const double kNaKNa = parameters[28];
  const double ICaPmax = parameters[29];
  const double kCaP = parameters[30];
  const double gam = parameters[31];
  const double dNaCa = parameters[32];
  const double gIf = parameters[33];
  const double DCa = parameters[34];
  const double DCaSR = parameters[35];
  const double DCaBm = parameters[36];
  const double DNa = parameters[37];
  const double Kmf_PLBKO = parameters[38];
  const double Kmf_PLB = parameters[39];
  const double Kmf_SLN = parameters[40];
  const double Kmr_PLBKO = parameters[41];
  const double Kmr_PLB = parameters[42];
  const double Kmr_SLN = parameters[43];
  const double k4 = parameters[44];
  const double kSRleak = parameters[45];

  // Expressions for the hAM_KSMT_nSR component
  const double past = stim_period*std::floor(t/stim_period);
  const double ist = (-past + t <= stim_duration + stim_offset && -past + t
    >= stim_offset ? stim_amplitude : 0.);
  const double cAF_lcell = 1.;
  const double cAF_gCaL = 1.;
  const double cAF_gt = 1.;
  const double cAF_gsus = 1.;
  const double cAF_gK1 = 1.;
  const double cAF_kNaCa = 1.;
  const double cAF_cpumps = 1.;
  const double cAF_PLB = 1.;
  const double cAF_SLN = 1.;
  const double cAF_phos = 1.;
  const double cAF_RyR = 1.;
  const double F = 96487.;
  const double R = 8314.;
  const double T = 306.150000000000;
  const double Ddcell = -1. + 2.*cAF_lcell;
  const double Dvcell = cAF_lcell*(Ddcell*Ddcell);
  const double Vss = 4.99232e-5*Dvcell;
  const double rjunct = 6.5*Ddcell;
  const double lcell = 122.051*cAF_lcell;
  const double dx = 1.625*Ddcell;
  const double Aj_nj = 1.0*M_PI*lcell*rjunct;
  const double xj_nj = 0.01*Ddcell + dx/2.;
  const double xj_nj_Nai = 2.*dx + 0.01*Ddcell;
  const double Vnonjunct1 = 5.0e-7*M_PI*(dx*dx)*lcell;
  const double Vnonjunct2 = 3.*Vnonjunct1;
  const double Vnonjunct3 = 5.*Vnonjunct1;
  const double Vnonjunct4 = 7.*Vnonjunct1;
  const double Vcytosol = Vss + 16.*Vnonjunct1;
  const double VSR1 = 0.0225*Vnonjunct1/Dvcell;
  const double VSR2 = 0.0225*Vnonjunct2/Dvcell;
  const double VSR3 = 0.0225*Vnonjunct3/Dvcell;
  const double VSR4 = 0.0225*Vnonjunct4/Dvcell;
  const double Vnonjunct_Nai = 16.*Vnonjunct1;
  const double gCaL = 15.*cAF_gCaL;
  const double gt = 8.25*cAF_gt;
  const double gsus = 2.25*cAF_gsus;
  const double gKr = 0.5*cAF_gK1;
  const double kNaCa = 0.0084*cAF_kNaCa;
  const double base_phos = 0.1*cAF_phos;
  const double PLB_SERCA_ratio = cAF_PLB;
  const double SLN_SERCA_ratio = cAF_SLN;
  const double SERCAKmf = Kmf_PLB*(1. - base_phos)*PLB_SERCA_ratio +
    Kmf_PLBKO + Kmf_SLN*(1. - base_phos)*SLN_SERCA_ratio;
  const double SERCAKmr = -Kmr_PLB*(1. - base_phos)*PLB_SERCA_ratio -
    Kmr_SLN*(1. - base_phos)*SLN_SERCA_ratio + Kmr_PLBKO;
  const double k3 = k4/(SERCAKmr*SERCAKmr);
  const double k1 = 1000000.*k4;
  const double k2 = (SERCAKmf*SERCAKmf)*k1;
  const double cpumps = 0.03*cAF_cpumps/Dvcell;
  const double ENa = R*T*std::log(Nao/Nass)/F;
  const double EK = R*T*std::log(Ko/Ki)/F;
  const double ECa = R*T*std::log(Cao/Cass)/(2.*F);
  const double INa = Nao*PNa*(F*F)*(INam*INam*INam)*(-1. + std::exp(F*(V -
    ENa)/(R*T)))*(0.9*INah1 + 0.1*INah2)*V/(R*T*(-1. + std::exp(F*V/(R*T))));
  const double INaminf = 1.0/(1. +
    0.0367620699824072*std::exp(-0.121802679658952*V));
  const double INahinf = 1.0/(1. +
    162754.791419004*std::exp(0.188679245283019*V));
  const double INamtau = 2.4e-5 + 4.2e-5*std::exp(-((0.887847222222222 +
    0.0347222222222222*V)*(0.887847222222222 + 0.0347222222222222*V)));
  const double INah1tau = 0.0003 + 0.03/(1. +
    58032.0080361162*std::exp(0.3125*V));
  const double INah2tau = 0.003 + 0.12/(1. +
    58032.0080361162*std::exp(0.3125*V));
  const double ICaLfcainf = 1. - 1./(1. + std::pow(kCa/Cass, kCan));
  const double ICaL = (V - ECa_app)*ICaLd*ICaLf2*ICaLfca*gCaL;
  const double ICaLdinf = 1.0/(1. +
    0.211882344332697*std::exp(-0.172413793103448*V));
  const double ICaLfinf = 0.04 + 1.0/(1. +
    1808.04241445606*std::exp(-0.125*V)) + 0.96/(1. +
    20.815841048468*std::exp(0.119047619047619*V));
  const double ICaLdtau = 0.0005 + 0.00065*std::exp(-((7./6. + V/30.)*(7./6.
    + V/30.)));
  const double ICaLf1tau = -0.00821 + 0.03576/(1. +
    9.36873231613932*std::exp(0.0756568146362647*V)) + 0.04275/(1. +
    1.89702279017873*std::exp(-0.0124048993894433*V)) +
    0.98698*std::exp(-((4.25157034998788 +
    0.1409649899351*V)*(4.25157034998788 + 0.1409649899351*V)));
  const double ICaLf2tau = 0.04 + 1.34*std::exp(-((2.8169014084507 +
    0.0704225352112676*V)*(2.8169014084507 + 0.0704225352112676*V)));
  const double ICaLfcatau = 0.00200000000000000;
  const double It = (-EK + V)*Itr*Its*gt;
  const double Itrinf = 1.0/(1. + std::exp(1./11. - V/11.));
  const double Itsinf = 1.0/(1. +
    33.8432351130073*std::exp(0.0869565217391304*V));
  const double Itrtau = 0.0015 + 0.0035*std::exp(-(V*V)/900.);
  const double Itstau = 0.01414 + 0.025635*std::exp(-((3.30233524526686 +
    0.0629615871356885*V)*(3.30233524526686 + 0.0629615871356885*V)));
  const double Isus = (-EK + V)*Isusr*Isuss*gsus;
  const double Isusrinf = 1.0/(1. +
    0.49774149722499*std::exp(-0.116279069767442*V));
  const double Isussinf = 1.0/(1. + 2.11700001661267*std::exp(V/10.));
  const double Isusrtau = 0.0005 + 0.009/(1. + std::exp(5./12. + V/12.));
  const double Isusstau = 3.05 + 0.59/(1. + std::exp(6. + V/10.));
  const double IKs = gKs*(-EK + V)*IKsn;
  const double IKsninf = 1.0/(1. +
    4.79191026127248*std::exp(-0.078740157480315*V));
  const double IKsntau = 0.7 + 0.4*std::exp(-((-1. + V/20.)*(-1. + V/20.)));
  const double IKrpi = 1.0/(1. + std::exp(55./24. + V/24.));
  const double IKr = (-EK + V)*IKrpa*IKrpi*gKr;
  const double IKrpainf = 1.0/(1. + std::exp(-5./2. - V/6.));
  const double IKrpatau = 0.03118 + 0.21718*std::exp(-((0.907115443521505 +
    0.0450458566821024*V)*(0.907115443521505 + 0.0450458566821024*V)));
  const double IK1 = gK1*std::pow(Ko, 0.4457)*(-EK + V)/(1. + std::exp(F*(5.4 +
    1.5*V - 1.5*EK)/(R*T)));
  const double INab = gNab*(V - ENa);
  const double ICab = gCab*(V - ECa);
  const double INaK = INaKmax*Ko*std::pow(Nass, 1.5)*(150. + V)/((200. +
    V)*(Ko + kNaKK)*(std::pow(Nass, 1.5) + std::pow(kNaKNa, 1.5)));
  const double fCaNCX = 1.;
  const double INaCa = (Cao*(Nass*Nass*Nass)*std::exp(F*gam*V/(R*T)) -
    fCaNCX*(Nao*Nao*Nao)*Cass*std::exp(F*(-1. + gam)*V/(R*T)))*kNaCa/(1. +
    dNaCa*(fCaNCX*(Nao*Nao*Nao)*Cass + Cao*(Nass*Nass*Nass)));
  const double ICaP = ICaPmax*Cass/(Cass + kCaP);
  const double Ifyinf = 1.0/(1. +
    2536.86493784617*std::exp(0.0801266000280443*V));
  const double Ifytau = 1.0/(0.00332*std::exp(-0.0604557273640154*V) +
    23.71839*std::exp(0.0604557273640154*V));
  const double IfNa = gIf*(-0.2677*ENa + 0.2677*V)*Ify;
  const double IfK = gIf*(0.7323*V - 0.7323*EK)*Ify;
  const double If = IfNa + IfK;
  const double betass = 1.0/(1. + KdSLhigh*SLhigh/((KdSLhigh +
    Cass)*(KdSLhigh + Cass)) + KdSLlow*SLlow/((KdSLlow + Cass)*(KdSLlow +
    Cass)) + BCa*KdBCa/((KdBCa + Cass)*(KdBCa + Cass)));
  const double betai1 = 1.0/(1. + BCa*KdBCa/((KdBCa + Cai1)*(KdBCa + Cai1)));
  const double betai2 = 1.0/(1. + BCa*KdBCa/((KdBCa + Cai2)*(KdBCa + Cai2)));
  const double betai3 = 1.0/(1. + BCa*KdBCa/((KdBCa + Cai3)*(KdBCa + Cai3)));
  const double betai4 = 1.0/(1. + BCa*KdBCa/((KdBCa + Cai4)*(KdBCa + Cai4)));
  const double gammai1 = BCa*KdBCa/((KdBCa + Cai1)*(KdBCa + Cai1));
  const double gammai2 = BCa*KdBCa/((KdBCa + Cai2)*(KdBCa + Cai2));
  const double gammai3 = BCa*KdBCa/((KdBCa + Cai3)*(KdBCa + Cai3));
  const double gammai4 = BCa*KdBCa/((KdBCa + Cai4)*(KdBCa + Cai4));
  const double betaSR1 = 1.0/(1. + CSQN*KdCSQN/((CaSR1 + KdCSQN)*(CaSR1 +
    KdCSQN)));
  const double betaSR2 = 1.0/(1. + CSQN*KdCSQN/((CaSR2 + KdCSQN)*(CaSR2 +
    KdCSQN)));
  const double betaSR3 = 1.0/(1. + CSQN*KdCSQN/((CaSR3 + KdCSQN)*(CaSR3 +
    KdCSQN)));
  const double betaSR4 = 1.0/(1. + CSQN*KdCSQN/((CaSR4 + KdCSQN)*(CaSR4 +
    KdCSQN)));
  const double betaNass = 1.0/(1. + BNa*KdBNa/((Nass + KdBNa)*(Nass + KdBNa)));
  const double Jj_nj = 1.0e-6*DCa*(Cass - Cai4)*Aj_nj/xj_nj;
  const double J_SERCASR1 = 2.*(k4*SERCACa1 - (CaSR1*CaSR1)*(cpumps -
    SERCACa1)*k3)*Vnonjunct1;
  const double J_bulkSERCA1 = 2.*((Cai1*Cai1)*(cpumps - SERCACa1)*k1 -
    SERCACa1*k2)*Vnonjunct1;
  const double J_SERCASR2 = 2.*(-(CaSR2*CaSR2)*(cpumps - SERCACa2)*k3 +
    k4*SERCACa2)*Vnonjunct2;
  const double J_bulkSERCA2 = 2.*(-SERCACa2*k2 + (Cai2*Cai2)*(cpumps -
    SERCACa2)*k1)*Vnonjunct2;
  const double J_SERCASR3 = 2.*(-(CaSR3*CaSR3)*(cpumps - SERCACa3)*k3 +
    k4*SERCACa3)*Vnonjunct3;
  const double J_bulkSERCA3 = 2.*((Cai3*Cai3)*(cpumps - SERCACa3)*k1 -
    SERCACa3*k2)*Vnonjunct3;
  const double J_SERCASRss = 2.*(-(CaSR4*CaSR4)*(cpumps - SERCACass)*k3 +
    k4*SERCACass)*Vss;
  const double J_bulkSERCAss = 2.*(-SERCACass*k2 + (Cass*Cass)*(cpumps -
    SERCACass)*k1)*Vss;
  const double RyRtauadapt = 1.;
  const double RyRtauactss = 0.00500000000000000;
  const double RyRtauinactss = 0.0150000000000000;
  const double nuss = 625.*Vss/Dvcell;
  const double RyRSRCass = 1. - 1./(1. + std::exp(-3.0/cAF_RyR + 10.0*CaSR4));
  const double RyRainfss = 0.505 - 0.427/(1. +
    0.0291125663159639*std::exp(12195.1219512195*Cass));
  const double RyRoinfss = 1. - 1./(1. + std::exp(33333.3333333333*Cass -
    7.33333333333333/cAF_RyR - 33.3333333333333*RyRass));
  const double RyRcinfss = 1.0/(1. + 0.135335283236613*std::exp(-100.0*RyRass
    + 100000.0*Cass));
  const double Jrelss = (CaSR4 - Cass)*RyRSRCass*RyRcss*RyRoss*nuss;
  const double RyRtauact = 0.0187500000000000;
  const double RyRtauinact = 0.0875000000000000;
  const double nu1 = Vnonjunct1/Dvcell;
  const double RyRSRCa1 = 1. - 1./(1. + std::exp(10.0*CaSR1 - 3.0/cAF_RyR));
  const double RyRainf1 = 0.505 - 0.427/(1. +
    0.0291125663159639*std::exp(12195.1219512195*Cai1));
  const double RyRoinf1 = 1. - 1./(1. + std::exp(-7.33333333333333/cAF_RyR +
    33333.3333333333*Cai1 - 33.3333333333333*RyRa1));
  const double RyRcinf1 = 1.0/(1. + 0.135335283236613*std::exp(100000.0*Cai1 -
    100.0*RyRa1));
  const double Jrel1 = (CaSR1 - Cai1)*RyRSRCa1*RyRc1*RyRo1*nu1;
  const double nu2 = Vnonjunct2/Dvcell;
  const double RyRSRCa2 = 1. - 1./(1. + std::exp(-3.0/cAF_RyR + 10.0*CaSR2));
  const double RyRainf2 = 0.505 - 0.427/(1. +
    0.0291125663159639*std::exp(12195.1219512195*Cai2));
  const double RyRoinf2 = 1. - 1./(1. + std::exp(-7.33333333333333/cAF_RyR +
    33333.3333333333*Cai2 - 33.3333333333333*RyRa2));
  const double RyRcinf2 = 1.0/(1. + 0.135335283236613*std::exp(-100.0*RyRa2 +
    100000.0*Cai2));
  const double Jrel2 = (CaSR2 - Cai2)*RyRSRCa2*RyRc2*RyRo2*nu2;
  const double nu3 = Vnonjunct3/Dvcell;
  const double RyRSRCa3 = 1. - 1./(1. + std::exp(10.0*CaSR3 - 3.0/cAF_RyR));
  const double RyRainf3 = 0.505 - 0.427/(1. +
    0.0291125663159639*std::exp(12195.1219512195*Cai3));
  const double RyRoinf3 = 1. - 1./(1. + std::exp(-7.33333333333333/cAF_RyR -
    33.3333333333333*RyRa3 + 33333.3333333333*Cai3));
  const double RyRcinf3 = 1.0/(1. + 0.135335283236613*std::exp(-100.0*RyRa3 +
    100000.0*Cai3));
  const double Jrel3 = (CaSR3 - Cai3)*RyRSRCa3*RyRc3*RyRo3*nu3;
  const double JSRCaleak1 = kSRleak*(CaSR1 - Cai1)*Vnonjunct1/Dvcell;
  const double JSRCaleak2 = kSRleak*(CaSR2 - Cai2)*Vnonjunct2/Dvcell;
  const double JSRCaleak3 = kSRleak*(CaSR3 - Cai3)*Vnonjunct3/Dvcell;
  const double JSRCaleakss = kSRleak*(CaSR4 - Cass)*Vss/Dvcell;
  const double JCa1 = JSRCaleak1 + Jrel1 - J_bulkSERCA1;
  const double JCa2 = -J_bulkSERCA2 + Jrel2 + JSRCaleak2;
  const double JCa3 = JSRCaleak3 - J_bulkSERCA3 + Jrel3;
  const double JCa4 = Jj_nj;
  const double JCass = Jrelss - Jj_nj - J_bulkSERCAss + JSRCaleakss;
  const double JSRCa1 = -Jrel1 + J_SERCASR1 - JSRCaleak1;
  const double JSRCa2 = -Jrel2 - JSRCaleak2 + J_SERCASR2;
  const double JSRCa3 = -JSRCaleak3 + J_SERCASR3 - Jrel3;
  const double JSRCa4 = J_SERCASRss - Jrelss - JSRCaleakss;
  const double JNa = 1.0e-6*DNa*(Nass - Nai)*Aj_nj/xj_nj_Nai;
  values[0] = -(0.001*INab + 0.001*INa + 0.001*IK1 + 0.001*If + 0.001*ist +
    0.001*ICaP + 0.001*ICab + 0.001*IKs + 0.001*It + 0.001*INaCa + 0.001*IKr
    + 0.001*INaK + 0.001*ICaL + 0.001*Isus)/Cm;
  values[1] = (0.001*INaminf - 0.001*INam)/INamtau;
  values[2] = (-0.001*INah1 + 0.001*INahinf)/INah1tau;
  values[3] = (-0.001*INah2 + 0.001*INahinf)/INah2tau;
  values[4] = (-0.001*ICaLd + 0.001*ICaLdinf)/ICaLdtau;
  values[5] = (0.001*ICaLfinf - 0.001*ICaLf1)/ICaLf1tau;
  values[6] = (0.001*ICaLfinf - 0.001*ICaLf2)/ICaLf2tau;
  values[7] = (0.001*ICaLfcainf - 0.001*ICaLfca)/ICaLfcatau;
  values[8] = (-0.001*Itr + 0.001*Itrinf)/Itrtau;
  values[9] = (-0.001*Its + 0.001*Itsinf)/Itstau;
  values[10] = (-0.001*Isusr + 0.001*Isusrinf)/Isusrtau;
  values[11] = (-0.001*Isuss + 0.001*Isussinf)/Isusstau;
  values[12] = (-0.001*IKsn + 0.001*IKsninf)/IKsntau;
  values[13] = (-0.001*IKrpa + 0.001*IKrpainf)/IKrpatau;
  values[14] = (0.001*Ifyinf - 0.001*Ify)/Ifytau;
  values[15] = (0.001*RyRoinfss - 0.001*RyRoss)/RyRtauactss;
  values[16] = (-0.001*RyRcss + 0.001*RyRcinfss)/RyRtauinactss;
  values[17] = (-0.001*RyRass + 0.001*RyRainfss)/RyRtauadapt;
  values[18] = (0.001*RyRoinf1 - 0.001*RyRo1)/RyRtauact;
  values[19] = (0.001*RyRcinf1 - 0.001*RyRc1)/RyRtauinact;
  values[20] = (0.001*RyRainf1 - 0.001*RyRa1)/RyRtauadapt;
  values[21] = (-0.001*RyRo2 + 0.001*RyRoinf2)/RyRtauact;
  values[22] = (-0.001*RyRc2 + 0.001*RyRcinf2)/RyRtauinact;
  values[23] = (-0.001*RyRa2 + 0.001*RyRainf2)/RyRtauadapt;
  values[24] = (-0.001*RyRo3 + 0.001*RyRoinf3)/RyRtauact;
  values[25] = (-0.001*RyRc3 + 0.001*RyRcinf3)/RyRtauinact;
  values[26] = (-0.001*RyRa3 + 0.001*RyRainf3)/RyRtauadapt;
  values[27] = (-0.0005*J_SERCASR1 + 0.0005*J_bulkSERCA1)/Vnonjunct1;
  values[28] = (0.0005*J_bulkSERCA2 - 0.0005*J_SERCASR2)/Vnonjunct2;
  values[29] = (0.0005*J_bulkSERCA3 - 0.0005*J_SERCASR3)/Vnonjunct3;
  values[30] = (0.0005*J_bulkSERCAss - 0.0005*J_SERCASRss)/Vss;
  values[31] = 0.001*(-JNa/Vss - (INab + 3.*INaK + INa + IfNa +
    3.*INaCa)/(F*Vss))*betaNass;
  values[32] = 0.001*JNa/Vnonjunct_Nai;
  values[33] = (-0.001*Isus - 0.001*IK1 - 0.001*IKr - 0.001*IKs - 0.001*It -
    0.001*IfK + 0.002*INaK - 0.001*ist)/(F*Vcytosol);
  values[34] = 0.001*((-ICab - ICaL + 2.*INaCa - ICaP)/(2.*F*Vss) +
    JCass/Vss)*betass;
  values[35] = -0.0005*DCaBm*((-Cai1 + Cai2)*(-Cai1 +
    Cai2))*betai1*gammai1/((KdBCa + Cai1)*(dx*dx)) +
    0.001*JCa1*betai1/Vnonjunct1 + 0.0015*(-Cai1 + Cai2)*(DCaBm*gammai1 +
    DCa)*betai1/(dx*dx);
  values[36] = 0.001*(DCa + DCaBm*gammai2)*((-2.*Cai2 + Cai1 + Cai3)/(dx*dx)
    + (-Cai1 + Cai3)/(4.*(dx*dx)))*betai2 - 0.0005*DCaBm*((-Cai1 +
    Cai3)*(-Cai1 + Cai3))*betai2*gammai2/((KdBCa + Cai2)*(dx*dx)) +
    0.001*JCa2*betai2/Vnonjunct2;
  values[37] = 0.001*JCa3*betai3/Vnonjunct3 - 0.0005*DCaBm*((-Cai2 +
    Cai4)*(-Cai2 + Cai4))*betai3*gammai3/((KdBCa + Cai3)*(dx*dx)) +
    0.001*(DCaBm*gammai3 + DCa)*((-2.*Cai3 + Cai2 + Cai4)/(dx*dx) + (-Cai2 +
    Cai4)/(6.*(dx*dx)))*betai3;
  values[38] = 0.001*JCa4*betai4/Vnonjunct4 + 0.001*(DCaBm*gammai4 +
    DCa)*((-Cai3 + Cai4)/(8.*(dx*dx)) + (-Cai4 + Cai3)/(dx*dx))*betai4 -
    0.0005*DCaBm*((-Cai3 + Cai4)*(-Cai3 + Cai4))*betai4*gammai4/((KdBCa +
    Cai4)*(dx*dx));
  values[39] = 0.001*JSRCa1*betaSR1/VSR1 + 0.0015*DCaSR*(CaSR2 -
    CaSR1)*betaSR1/(dx*dx);
  values[40] = 0.001*JSRCa2*betaSR2/VSR2 + 0.001*DCaSR*((CaSR3 -
    CaSR1)/(4.*(dx*dx)) + (CaSR3 + CaSR1 - 2.*CaSR2)/(dx*dx))*betaSR2;
  values[41] = 0.001*JSRCa3*betaSR3/VSR3 + 0.001*DCaSR*((CaSR4 + CaSR2 -
    2.*CaSR3)/(dx*dx) + (CaSR4 - CaSR2)/(6.*(dx*dx)))*betaSR3;
  values[42] = 0.001*DCaSR*((CaSR4 - CaSR3)/(8.*(dx*dx)) + (CaSR3 -
    CaSR4)/(dx*dx))*betaSR4 + 0.001*JSRCa4*betaSR4/VSR4;
}
