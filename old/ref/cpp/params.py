from math import log, exp
inf = 1E-18
sup = 1000

with open('xfinal.dat', 'r') as infile:
    for line in infile.readlines():
        para = float(line)
        para = (para - inf)/(sup - inf)
        #print para
        para1 = log(para/(1-para))

        para=exp(para1)/(exp(para1)+1);
        para=inf+para*(sup-inf);
        print para

    # for(int i=0;i<NDIM;i++)
    # {
    #     infile1>>para[i];
    #     para[i]=(para[i]-inf_limit[i])/(sup_limit[i]-inf_limit[i]); //normalaze parameters
    #     para1[i]=log(para[i]/(1-para[i]));
    #     //printf("%f \n",para[i]);