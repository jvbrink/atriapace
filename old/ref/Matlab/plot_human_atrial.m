function [t, camean] = plot_human_atrial( filename )

data = load(filename);
 
%**********************************************
% Definition of differential variables
%**********************************************

i_V = 1;
i_lastend = i_V;

% INa
i_start = i_lastend + 1;
i_INam = i_start;
i_INah1 = i_start + 1;
i_INah2 = i_start + 2;
i_lastend = i_INah2;

% ICaL 
i_start = i_lastend + 1;
i_ICaLd = i_start;
i_ICaLf1 = i_start + 1;
i_ICaLf2 = i_start + 2;
i_ICaLfca = i_start + 3;
i_lastend = i_ICaLfca;

% It
i_start = i_lastend + 1;
i_Itr = i_start;
i_Its = i_start + 1;
i_lastend = i_Its;

% Isus
i_start = i_lastend + 1;
i_Isusr = i_start;
i_Isuss = i_start + 1;
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
i_RyRoss = i_start;
i_RyRcss = i_start + 1;
i_RyRass = i_start + 2;
i_RyRo1 = i_start + 3;
i_RyRc1 = i_start + 4;
i_RyRa1 =  i_start + 5;
i_RyRo2 = i_start + 6;
i_RyRc2 = i_start + 7;
i_RyRa2 =  i_start + 8;
i_RyRo3 = i_start + 9;
i_RyRc3 = i_start + 10;
i_RyRa3 =  i_start + 11;
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
i_Nass = i_start;
i_Nai = i_start + 1;
i_Ki = i_start + 2;
i_lastend = i_Ki;


% Cai and CaSR 
i_start = i_lastend + 1;
i_Cass = i_start;
i_Cacenter = i_start + 1;
  
    t = data(:,1);
    dydata = data(:,2:(29 + 6 + 4*2+1));
    
    currentdata =data(1:end,(end-24):(end));
    INa = currentdata(:,1); 
    ICaL = currentdata(:,2);
    It = currentdata(:,3); 
    Isus = currentdata(:,4);
    IKr = currentdata(:,5);
    IKs = currentdata(:,6);
    IK1 = currentdata(:,7);
    INab = currentdata(:,8);
    ICab = currentdata(:,9);
    ICaP = currentdata(:,10);
    INaK = currentdata(:,11);
    INaCa = currentdata(:,12);
    If = currentdata(:,13); 
    Jrelss = currentdata(:,14);
    Jrel1 = currentdata(:,15);
    Jrel2 = currentdata(:,16);
    Jrel3= currentdata(:,17);
    J_bulkSERCA1 = currentdata(:,18);
    J_bulkSERCA2 = currentdata(:,19);
    J_bulkSERCA3 = currentdata(:,20);
    J_bulkSERCAss = currentdata(:,21); 
    JSRCaleak1 = currentdata(:,22);
    JSRCaleak2 = currentdata(:,23); 
    JSRCaleak3 = currentdata(:,24);
    JSRCaleakss = currentdata(:,25);  


  cadata = dydata(1:end,i_Cacenter:i_Cacenter+3);
  caSRdata = dydata(1:end,i_Cacenter+4:i_Cacenter+7);
  cassdata = dydata(1:end,i_Cass);
  x = [0.8125:1.625:5.6875 6.51];
  ca = [cadata cassdata];

  figure
 [ch, ch] = contourf(t,x,ca',100);
  colormap(jet)
     set(ch,'edgecolor','none'); 
  c = colorbar;
  %contourf(t,x,ca,20);
title(filename)

figure
  subplot(311)
  plot(t,cadata(1:end,1),t,cadata(1:end,2),t,cadata(1:end,3),t,cadata(1:end,4),t,cassdata)
  title(['Ca data 1-4 + ss ' filename])
  legend('ca1','ca2','ca3','ca4','cass')
  subplot(312)
  camean = (cadata(1:end,1).*1.625 + cadata(1:end,2).*1.625 + cadata(1:end,3).*1.625 + cadata(1:end,4).*1.625 + cassdata.*0.02)./6.52 ;
  cameannj = mean(cadata');
  
  plot(t,camean,t,cameannj)
  legend('camean','cameannj')
  title('mean ca data 1-4')
  subplot(313)
  plot(t,dydata(1:end,i_SERCACa1),t,dydata(1:end,i_SERCACa2),t,dydata(1:end,i_SERCACa3),t,dydata(1:end,i_SERCACass)); % SERCA Ca bound
  title('SERCA Ca bound 1-3 + ss');
  
figure
  subplot(341)
  plot(t,dydata(1:end,i_V));
  title(['Vm' filename])
  subplot(342)
  plot(t,dydata(1:end,i_Nass), t,dydata(1:end,i_Nai))
  legend('Nass', 'Nai')
    subplot(343)
  plot(t,dydata(1:end,i_Ki))
  legend('Ki')
  subplot(344)
  plot(t,If)
  legend('If')
  subplot(345)
  plot(t,INaCa)
  legend('INCX')
  subplot(346)
  plot(t,INaK)
  legend('INaK')
      subplot(347)
  plot(t,ICaP)
  legend('ICaP')
        subplot(348)
  plot(t,ICab)
  legend('ICab')
          subplot(349)
  plot(t,INab)
  legend('INab')
          subplot(3,4,10)
  plot(t,IK1,t,IKs,t,IKr,t,Isus,t,It)
  legend('IK1','IKs','IKr','Isus','It' )
           subplot(3,4,11)
  plot(t,ICaL)
  legend('ICaL') 
             subplot(3,4,12)
  plot(t,INa)
  legend('INa') 
  
figure

  subplot(421)
  plot(t,caSRdata(1:end,1),t,caSRdata(1:end,2),t,caSRdata(1:end,3),t,caSRdata(1:end,4))
  legend('SRca1','SRca2','SRca3','SRCa4 (rel to ss)');
    title(filename)

  subplot(422)
  plot(t,Jrel1,t,Jrel2,t,Jrel3,t,Jrelss)
  legend('Jrel1','Jrel2','Jrel3','Jrelss')
  subplot(423)
  plot(t,J_bulkSERCA1,t,J_bulkSERCA2,t,J_bulkSERCA3,t,J_bulkSERCAss)
  legend('Serca1','Serca2','Serca3','Sercass')
  subplot(424)
  plot(t,JSRCaleak1,t,JSRCaleak2,t,JSRCaleak3,t,JSRCaleakss)
     legend('SRleak1','SRleak2','SRleak3','SRleakss')
  subplot(425)
  plot(t,dydata(:,i_RyRoss),t,dydata(:,i_RyRcss),t,dydata(:,i_RyRass));
  legend('RyRoss','RyRcss','RyRass')
    subplot(426)
  plot(t,dydata(:,i_RyRo1),t,dydata(:,i_RyRc1),t,dydata(:,i_RyRa1))
  legend('RyRo1','RyRc1','RyRa1')
      subplot(427)
  plot(t,dydata(:,i_RyRo2),t,dydata(:,i_RyRc2),t,dydata(:,i_RyRa2))
  legend('RyRo2','RyRc2','RyRa2')
        subplot(428)
  plot(t,dydata(:,i_RyRo3),t,dydata(:,i_RyRc3),t,dydata(:,i_RyRa3))
  legend('RyRo3','RyRc3','RyRa3')

    


