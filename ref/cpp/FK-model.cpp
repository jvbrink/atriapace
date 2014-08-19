#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include "nr3.h" //for specially defined variables (doub, vecdoub, matdoub)

#define NDIM 21 //number of parameters
#define EPS 1e-18 
#define NX 100 //number of divisions in cable ("cells")
#define NEXCITE 10 //number of excitations from stimulus current
#define TOTALTIME 30000
#define DT 0.05
#define DX 0.02
#define SPACING_INDEX 10
#define S1limit 5

Doub heav(Doub u); 

using namespace std;
Doub heav(Doub u) //Heaviside function
{
    if (u>0.0)
        return 1.0;
    else
        return 0.0;
}

struct solver_wave
{
    Doub savedt;
    int numt;
    bool success;
    int savepoints;
    VecInt nexcs;
    VecDoub DI;
    VecDoub u,v,w,d,xfi,xso,xsi;
    //VecDoub tso;
    Doub dv,dw,dd,tso;
    VecDoub para,para1,excitimes;
    //MatDoub voltage;
    VecDoub v70;
    VecDoub v30;
    VecDoub inf;
    VecDoub sup;
    
    //----------------KNOW PARAMETERS------------------//
    Doub uo, um, una;
    //---------------UNKNOWN PARAMETERS-------------------//
    Doub uc,uv,uw,ud,tvm,tvp,twm,twp,tsp,tsm,ucsi,xk,td,to,tsoa,tsob,uso,xtso,tsi,D,tvmm;
    
    solver_wave(VecDoub ppara1,Doub ssavedt,VecDoub dDI,VecDoub eexcitimes,VecDoub iinf, VecDoub ssup):para1(ppara1),savedt(ssavedt),DI(dDI),excitimes(eexcitimes),inf(iinf),sup(ssup),para(NDIM,0.0),u(NX,0.0),v(NX,0.0),w(NX,0.0),d(NX,0.0),xfi(NX,0.0),xso(NX,0.0),xsi(NX,0.0),nexcs(NEXCITE,0),v70(100),v30(100)
    {
        //tvmm=10;
        uo=0.0;
        um=1.0;
        una=0.23;
        //t=0.0;
        
        savepoints=int(TOTALTIME/savedt);
        for (int j=0; j<NX; j++) {
            u[j]=0;
            v[j]=1;
            w[j]=1;
            d[j]=0;
        }
        for (int i=0; i<NDIM; i++) {
            para[i]=exp(para1[i])/(exp(para1[i])+1);
            para[i]=inf[i]+para[i]*(sup[i]-inf[i]);
            //printf("para[%i] is %f\n",i,para[i]);
        }
        numt=int(TOTALTIME/DT);
        
        uc=para[0]; //threshold
        uv=para[1]; //fast gate threshold, determines whethere tvm or tvmm is active (chaos 8)
        uw=para[2]; //slow gate threshold
        ud=para[3]; //threshold
        tvm=para[4]; //controls minimum diastolic interval where CV occurs (chaos 8)
        tvp=para[5]; //fast gate closing time
        twm=para[6]; //slow gate opening time (changes APD shape?)
        twp=para[7]; //slow gate closing time (shifts APD up/down?)
        tsp=para[8]; //d-gate variables
        tsm=para[9];
        ucsi=para[10];
        xk=para[11]; //typically around 10
        td=para[12]; //fast current time variable, determines max CV
        to=para[13]; //ungated time constant
        tsoa=para[14]; //curve shape/APD, ungated time, adjusts DI
        tsob=para[15]; //ungated time. Easily adjusts DI, changes APD
        uso=para[16];
        xtso=para[17];
        tsi=para[18]; //slow current time variable, max APD
        D=para[19]; //related to density, mostly changes CV, but can effect everything
        tvmm=para[20]; //controls the steepness of the CV curve (chaos 8)
        //um=para[21]; //related to maximum AP, set to 1
        //una=para[21]; //threshold. Um,una have little effect on APD
        
        for (int i=0; i<NEXCITE; i++) {
            nexcs[i]=int(excitimes[i]/DT);
            //printf("%f, %f \n",nexcs[i]*DT,excitimes[i]);
        }
        
        //voltage.resize(savepoints, NX);
        v70.resize(savepoints);
        v30.resize(savepoints);
        //for (int i=0; i<savepoints; i++)
        //for (int j=0; j<NX; j++)
        //{
        //   voltage[i][j]=0.0;
        //}
        success=true;
    }
    void solving();
};
void solver_wave::solving()
{
    double maxVoltage[NEXCITE];
    int i,j,k;
    int ni=0;
    int ns=0;
    double preV1 =0, preV2=0;
    VecDoub tempu(NX+1,0.0);
    
    int numS1 = 0; //count how many S1 stim was applied (have ~5 inbetween each S2 stim)
    //int S1limit = 1; //how many S1 values to apply before going to S2
    float S1time = 400./DT; //CL(cycle length) between S1 stim
    int lastStimTime = 0;

    
    FILE *fv = fopen("voltage.txt","w+"); //save the voltage to file
    
    for (int init=0;init<NEXCITE;init++){
        maxVoltage[init]=0.0;
    }
    
    for(i=0;i<numt;i++)
    { 
        if(i == S1time + lastStimTime && numS1 < S1limit) {
            numS1++;
            lastStimTime = i;
            u[0] = 0.5;
        }
        else if(i == nexcs[ni] && numS1 >= S1limit) {
            ni++;
            u[0]=0.5;
            lastStimTime = i;
            numS1 = 0;
        }
        
        for(j=0;j<NX;j++)
        {
            tso=tsoa+(tsob-tsoa)*(1+tanh((u[j]-uso)*xtso))/2.0;
            dv=(1-heav(u[j]-una))*(1-v[j])/((1-heav(u[j]-uv))*tvm + tvmm*heav(u[j]-uv)) - heav(u[j]-una)*v[j]/tvp;
            
            dw=(1-heav(u[j]-uw))*(1-w[j])/twm - heav(u[j]-uw)*w[j]/twp;
            
            dd=((1-heav(u[j]-ud))/tsm + heav(u[j]-ud)/tsp)*((1+tanh(xk*(u[j]-ucsi)))/2.0-d[j]);
            v[j]=v[j]+DT*dv;  //fast gate
            w[j]=w[j]+DT*dw;  //slow gate
            d[j]=d[j]+DT*dd;
            
            //----------currents---------//
            xfi[j]=-v[j]*heav(u[j]-una)*(u[j]-una)*(um-u[j])/td;
            xso[j]=(u[j]-uo)*(1-heav(u[j]-uc))/to + heav(u[j]-uc)/tso;
            xsi[j]=-w[j]*d[j]/tsi;
            
        }
        
        //check max voltage (used for finding excitimes)
        if(numS1>0) { //APPLY S1
            if(preV1 > preV2 && u[0] < preV1) { //max slope
                if(preV1 > maxVoltage[ni]) {
                    maxVoltage[ni] = preV1;
                }
            }
        }
        
        if (i==round(ns*savedt/DT)&&ns<savepoints)
        {
            v70[ns] = u[70];
            v30[ns] = u[30];
            
            fprintf(fv, "%f \n", u[70]);
            
            if (ns>0) { 
                if (ns>lastStimTime*DT/savedt+100 && ni<NEXCITE) { 
                    float baseline = 0.1*maxVoltage[ni]; 
                    
                    if(tempu[0]>baseline && u[0]<=baseline && numS1 >= S1limit)
                    {
                        if ((tempu[0]-baseline)>(baseline-u[0]))
                        { 
                            nexcs[ni] = int(ns*savedt/DT) + int((DI[ni])/DT); //APPLY S1
                        }
                        else
                        {
                            nexcs[ni] = int((ns-1)*savedt/DT) + int((DI[ni])/DT); 
                        }
                    }
                }
            }
                                    
            ns=ns+1;
        }
        
        for (k=0; k<NX; k++) {
            tempu[k]=u[k];
        }
        
        for(k=0;k<NX;k++)
        {
            if (k!=0 && k!=NX-1) {
                u[k]=tempu[k]+DT*(D*(tempu[k+1]+tempu[k-1]-2*tempu[k])/pow(DX,2.0)-(xfi[k]+xso[k]+xsi[k]));
            }
            else {
                if (k==0) {
                    preV2 = preV1;
                    preV1 = u[k];
                    u[k]=tempu[k]+DT*(D*(tempu[k+1]-tempu[k])/pow(DX,2.0)-(xfi[k]+xso[k]+xsi[k]));
                }
                else {
                    u[k]=tempu[k]+DT*(D*(tempu[k-1]-tempu[k])/pow(DX,2.0)-(xfi[k]+xso[k]+xsi[k]));
                }
            }
            
        }
    }
    fclose(fv);
}


int main(int argc, char *argv[])
{
    VecDoub para(NDIM),para1(NDIM);
    double runningtime;
	clock_t begin,finish;
    VecDoub inf_limit(NDIM),sup_limit(NDIM); //used in fitting procedure, not really necessary here
    const Doub del = log(1.1); //fitting
	const Doub ftol = 0.0;  //fitting
    Doub savedt = DT;
    
    inf_limit[0]=0;
    inf_limit[1]=EPS;
    inf_limit[2]=EPS;
    inf_limit[3]=EPS;
    inf_limit[4]=EPS;
    inf_limit[5]=EPS;
    inf_limit[6]=EPS;
    inf_limit[7]=EPS;
    inf_limit[8]=EPS;
    inf_limit[9]=EPS;
    inf_limit[10]=EPS;
    inf_limit[11]=EPS;
    inf_limit[12]=EPS;
    inf_limit[13]=EPS;
    inf_limit[14]=EPS;
    inf_limit[15]=EPS;
    inf_limit[16]=EPS;
    inf_limit[17]=EPS;
    inf_limit[18]=EPS;
    inf_limit[19]=0;
    inf_limit[20]=1;
    
    sup_limit[0]=1000;
    sup_limit[1]=1000;
    sup_limit[2]=1000;
    sup_limit[3]=1000;
    sup_limit[4]=1000;
    sup_limit[5]=1000;
    sup_limit[6]=2000;
    sup_limit[7]=1000;
    sup_limit[8]=1000;
    sup_limit[9]=1000;
    sup_limit[10]=1000;
    sup_limit[11]=2000;
    sup_limit[12]=1000;
    sup_limit[13]=1000;
    sup_limit[14]=1000;
    sup_limit[15]=1000;
    sup_limit[16]=1000;
    sup_limit[17]=1000;
    sup_limit[18]=1000;
    sup_limit[19]=1000;
    sup_limit[20]=2000;
    
    ifstream infile1("xfinal.dat",ios::in);
    
    if(!infile1)
    {
        cerr<<"open error!!"<<endl;
        exit(1);
    }
    for(int i=0;i<NDIM;i++)
    {
        infile1>>para[i];
        para[i]=(para[i]-inf_limit[i])/(sup_limit[i]-inf_limit[i]); //normalaze parameters
        para1[i]=log(para[i]/(1-para[i]));
        //printf("%f \n",para[i]);
    }
    infile1.close();
    
    int kmax=2000;
    VecDoub time(kmax,0.0),excitimes(NEXCITE,0.0),data(kmax,0.0),DI(NEXCITE,0.0);
    int length;
    
    for(int i=0;i<NEXCITE;i++) {
        excitimes[i]=TOTALTIME;
    }
    
    //load the experimental data
    int i=0;
    
    //for actual runs
    double maxMAP = -100.;
    double minMAP = 100000.;
    ifstream infile3("MAP.txt",ios::in);
    if(!infile3)
    {
        cerr<<"open error!"<<endl;
        exit(1);
    }
    while(infile3.eof()==0)
    {
        time[i]=i*0.4;
        infile3>>data[i];
        
        if (maxMAP<data[i])
            maxMAP = data[i];
        if (minMAP>data[i])
            minMAP = data[i];
        
        i++;
    }
    length=i-1;
    infile3.close();
    
    //normalize the MAP data
    for (int j=0;j<length;j++) {
        data[j] = (data[j]-minMAP)/(maxMAP-minMAP);
    }
    
    i=0;
    int j=0;
    
    ifstream infile5("DI.txt",ios::in);
    if(!infile5)
    {
        cerr<<"no DI data!"<<endl;
        exit(1);
    }    
    while(infile5.eof()==0)
    {
        infile5>>DI[j];
        j++;
    }
    infile5.close();
    
    solver_wave sol2(para1,savedt,DI,excitimes,inf_limit,sup_limit); //take parameters, get simulation results
    sol2.solving();
    
    return 0;
    
}

