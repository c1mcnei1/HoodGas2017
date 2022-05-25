% HoodGas GTD data processing 

% load data
% CLM modified for archiving 5/24/22

clear all;

cd('D:\DATA\NOAAShipGTD\DATA\HoodGas2017\ProcessedGTD17th');

load GTDloaded_17th;
GTD.Pgtd=GTD_P;
GTD.Tgtd=GTD_T;
GTD.MT=datenum(TeratermTime)-0.65/60/24; % coorecting for apparent 0.65min clock difference

load 3_3_loaded;
CTD.P=VarName1;
CTD.T=VarName2;
CTD.S=VarName3;
CTD.O2a=VarName4;
CTD.O2b=VarName5;
CTD.MT=datenum('Aug 17 2017')+(VarName11-229-7/24); % * System UTC = 7 hr ahead

save data CTD GTD


clear all;
load data

% interpolate 

[C,IA,IC] = unique(GTD.MT);
CTD.GT=interp1(GTD.MT(IA),GTD.Pgtd(IA),CTD.MT);

[C,IA,IC] = unique(CTD.MT);
GTD.P=interp1(CTD.MT(IA),CTD.P(IA),GTD.MT);
GTD.T=interp1(CTD.MT(IA),CTD.T(IA),GTD.MT);
GTD.S=interp1(CTD.MT(IA),CTD.S(IA),GTD.MT);
GTD.O2a=interp1(CTD.MT(IA),CTD.O2a(IA),GTD.MT);
GTD.O2b=interp1(CTD.MT(IA),CTD.O2b(IA),GTD.MT);
GTD.Min=24.*60.*(GTD.MT-CTD.MT(1));

% figures 

startMT=datenum(datestr('17-Aug-2017 14:00'));
endMT=datenum(datestr('17-Aug-2017 16:30'));
startMin=0;
endMin=150;
g=find(GTD.Min>=startMin & GTD.Min<=endMin);
figure(1);clf;
subplot(2,1,1);grid on;hold on;
 colormap(jet);
 scatter(GTD.Min(g),GTD.P(g),ones(size(GTD.Min(g))).*10,GTD.Min(g));
 set(gca,'ydir','reverse');
 ylabel('CTD P (dbar)');
 title('Hood Canal on RV Robertson on 08/17/2017');
 set(gca,'xlim',[startMin endMin]);
 set(gca,'ylim',[0 26]);
 %set(gca,'xticklabels',[]);
subplot(2,1,2);grid on;hold on;
 scatter(GTD.Min(g),GTD.Pgtd(g),ones(size(GTD.Min(g))).*10,GTD.Min(g));
 set(gca,'xlim',[startMin endMin]);
 ylabel('GTD reading (mbar)');
 xlabel('Time (minutes)');
 set(gca,'ylim',[850 1150]);
 print(1, '-djpeg', '-r300', 'fig1');

 figure(2);clf;
 suptitle('Hood Canal on RV Robertson on 08/17/2017');
 colormap(jet);
 subplot(1,4,1);grid on;hold on;
  scatter(GTD.T(g),GTD.P(g),ones(size(GTD.T(g))).*10,GTD.Min(g));
  scatter(GTD.T(g),GTD.P(g),ones(size(GTD.T(g))).*10,GTD.Min(g));
  set(gca,'ydir','reverse');
  ylabel('CTD P (dbar)');
  xlabel('CTD T (^oC)');
  set(gca,'xlim',[9 21]);
 subplot(1,4,2);grid on;hold on;
  scatter(GTD.S(g),GTD.P(g),ones(size(GTD.S(g))).*10,GTD.Min(g));
  set(gca,'ydir','reverse');
  set(gca,'yticklabels',[]);
  xlabel('CTD S (psu)');
  set(gca,'xlim',[25 30.5]);
subplot(1,4,3);grid on;hold on;
  scatter(GTD.O2a(g),GTD.P(g),ones(size(GTD.O2a(g))).*10,GTD.Min(g));
  set(gca,'ydir','reverse');
  set(gca,'yticklabels',[]);
  xlabel('CTD O_2 (umol/kg)');
  set(gca,'xlim',[50 350]);
 subplot(1,4,4);grid on;hold on;
  scatter(GTD.Pgtd(g),GTD.P(g),ones(size(GTD.P(g))).*10,GTD.Min(g));
  set(gca,'ydir','reverse');
  set(gca,'yticklabels',[]);
  xlabel('GTD P (mbar)');
  set(gca,'xlim',[850 1150]);
print(2, '-djpeg', '-r300', 'fig2');
  
 figure(3);clf;hold on;grid on;
 suptitle('Hood Canal on RV Robertson on 08/17/2017');
  colormap(jet);
  scatter(GTD.S(g),GTD.T(g),ones(size(GTD.S(g))).*10,GTD.Min(g));
  xlabel('S (psu)');
  ylabel('T (^oC)');
  set(gca,'xlim',[25 30.5]);
  set(gca,'ylim',[9 21]);
print(3, '-djpeg', '-r300', 'fig3');

save DATA17 CTD GTD;


cd('D:\DATA\NOAAShipGTD\DATA\HoodGas2017\ProcessedGTD17th');
clear all;
close all;
load DATA17;

data=GTD;
data.Sig0=sw_pden(data.S,data.T,data.P,data.P.*0)-1000;

d=datevec(data.MT);
data.YD=julian([d(:,1) d(:,2) d(:,3) d(:,4) d(:,5) d(:,6)])-julian(2017,0,1,0);

mt=data.MT;
x=data.Min;
z=data.Sig0;
o=data.O2a;
p=data.P;
t=data.T;
s=data.S;
w=psi_pwater(t,s);
y=data.Pgtd-w;
v=data.Tgtd;
g=find(x>10.5 & x<141 & y>800);
N=10;
xx=linspace(min(x(g)),max(x(g)),max(size(x(g))));
yy=interp1(x(g),smooth(y(g),N),xx);
zz=interp1(x(g),smooth(z(g),N),xx);
oo=interp1(x(g),smooth(o(g),N),xx);
pp=interp1(x(g),smooth(p(g),N),xx);
tt=interp1(x(g),smooth(t(g),N),xx);
ss=interp1(x(g),smooth(s(g),N),xx);
ww=interp1(x(g),smooth(w(g),N),xx);
vv=interp1(x(g),smooth(v(g),N),xx);
mtmt=interp1(x(g),smooth(mt(g),N),xx);

x=xx;clear xx;
y=yy;clear yy;
z=zz;clear zz;
o=oo;clear oo;
p=pp;clear pp;
t=tt;clear tt;
s=ss;clear ss;
w=ww;clear ww;
v=vv;clear vv;
mt=mtmt;clear mtmt;

g=find(x>0);

figure(50);clf;hold on;grid on;
plot(x,p,'r');
plot(data.Min,data.P,'.g');
plot(x(g),p(g),'b.');
set(gca,'ydir','reverse');
xlabel('Time (min)');
ylabel('P (dbar)');

figure(51);clf;hold on;grid on;
plot(x,y+w,'ro');
plot(x(g),y(g)+w(g),'b.');
dt=t(g)-v(g);
yc=y(g)./(t(g)+273.15).*(t(g)+273.15+dt);
plot(x(g),yc+w(g),'g.'); % (P1/T1).(T1+dT)=P2
dt=-t(g)+v(g);
yc=y(g)./(t(g)+273.15).*(t(g)+273.15+dt);
plot(x(g),yc+w(g),'c.'); % (P1/T1).(T1+dT)=P2
xlabel('Time (min)');
ylabel('Gas tension (mbar)');
title('Data (red)   Select Data(blue)');

y=yc; % !!!!!!!!!!!!!!!!!!!!!!!!

deltat=mean(diff(x));
Tau=1.8; % 1.9 min
h=(1/Tau)*exp(-(x-x(1))./Tau);
h=h./(sum(h)*deltat);

H=fft(h);
Y=fft(y);
M=Y./(H*deltat);
m=ifft(M);


colormap(jet);
xmin=0;xmax=150;
figure(52);clf;
subplot(3,1,1);hold on;grid on;
title('Hood Canal 17-Aug-2017: CTD (blue) GTD (red) Deconvolved (magenta)');
plot(x(g),p(g),'.b');
set(gca,'ydir','reverse');
set(gca,'ylim',[0 25]);
set(gca,'xlim',[xmin xmax]);
set(gca,'xticklabels',[]);
ylabel('P (dbar)');
subplot(3,1,2);hold on;grid on;
plot(x(g),t(g),'.b');
plot(x(g),v(g),'.r');
set(gca,'ylim',[8 22]);
set(gca,'xlim',[xmin xmax]);
set(gca,'xticklabels',[]);
ylabel('T (^oC)');
subplot(3,1,3);hold on;grid on;
xlabel('Time (min)');
ylabel('Gas tension (mbar)');
% % Argh- that's too noisy.  Try some Tikonov regularization.
Mtik=(conj(H).*Y)./((conj(H).*H+2)*deltat);
mtik=ifft(Mtik);
% [b,a] = butter(2,0.05,'low');           % IIR filter design
% ms = filtfilt(b,a,mtik);                    % zero-phase filtering
% ms=smooth(mtik,100);
% Argh- that's too noisy.  Try some Tikonov regularization.
%[b,a] = butter(2,0.075,'low');           % IIR filter design
%ms = filtfilt(b,a,m);                    % zero-phase filtering
ms=smooth(mtik,145)';
plot(x,ms+w,'m','linewidth',.5);
plot(x(g),ms(g)+w(g),'m.');
plot(x(g),y(g)+w(g),'.r');
set(gca,'ylim',[850 1150]);
set(gca,'xlim',[xmin xmax]);
a=find(diff(ms)>-0.05 & diff(ms)<0.05 );
gg=a;%[find(x>0)];
plot(x(gg),ms(gg)+w(gg),'mo');
print(52, '-djpeg', '-r300', '17th- deconvolved timeseries');


colormap(jet);
figure(53);clf;hold on; grid on;
cmin=0;cmax=150;
plot(y+w,p,'k.');
DotColorZ(ms(gg)+w(gg),p(gg),x(gg),'o',cmin,cmax);
colormap(jet);
H=colorbar;
set(gca,'clim',[cmin,cmax]);
title(H,'Time (min)');
xlabel('Gas tension (mbar)');
ylabel('P (dbar)');
set(gca,'ylim',[0 25]);
set(gca,'xlim',[850 1150]);
set(gca,'ydir','reverse');
title('Hood Canal 17-Aug-2017: deconvolved');
print(53, '-djpeg', '-r300', '17th- deconvolved profile');

% O2 sat
osat=100.*o./O2sol(t,s,p,'umol');
colormap(jet);
figure(54);clf;hold on; grid on;
cmin=20;cmax=135;
plot(y+w,osat,'k.');
colormap(jet);
DotColorZ(ms(gg)+w(gg),osat(gg),x(gg),'o',cmin,cmax);
colormap(jet);
H=colorbar;
set(gca,'clim',[cmin,cmax]);
title(H,'Time (min)');
xlabel('Gas tension (mbar)');
ylabel('O_2sat (%)');
set(gca,'ylim',[20 135]);
set(gca,'xlim',[850 1150]);
title('Hood Canal 17-Aug-2017: deconvolved');
print(54, '-djpeg', '-r300', '17th- deconvolved GT versus O2');

Xn2=0.78084;
Xo2=0.20946;
Xar=1-Xn2-Xo2;

clear GTD_calc*

% x=data.Min;
% z=data.Sig0;
% o=data.O2a;
% p=data.P;
% t=data.T;
% s=data.S;
% w=psi_pwater(t,s);
% y=data.Pgtd-w;
% v=data.Tgtd;
% g=find(x>10.5 & x<141 & y>800);

% ---------- N2 ---------------------------------

for i=1:max(size(gg)),
            MT=mt(gg(i));
            Tt=t(gg(i));
            Ss=s(gg(i));
            TDG=ms(gg(i))+w(gg(i));
            Pp=p(gg(i));
                             %!!! CLM: 10/13/07 potential density
            Dens=sw_pden(Ss,Tt,Pp,Pp.*0)./1000;
            
            C_O2c=o(gg(i)); % have to convert from umol/kg to ml/l
            C_O2c=C_O2c./O2sol(Tt,Ss,Pp,'umol')*O2sol(Tt,Ss,Pp,'ml');
            
            P_H2O = psi_pwater(Tt,Ss);                   % Saturated water vapor pressure in units of 'mbar' over salt water (at 100% RH)  

                  %!!! CLM: 10/13/07 potential temperature
                                             %!!! CLM: 10/13/07  Bunsen at surface pressure
            B_N2  = N2sol(sw_ptmp(Ss,Tt,Pp,Pp.*0),Ss,Pp.*0,'bunsen');
            B_O2  = O2sol(sw_ptmp(Ss,Tt,Pp,Pp.*0),Ss,Pp.*0,'bunsen');
            B_Ar  = Arsol(sw_ptmp(Ss,Tt,Pp,Pp.*0),Ss,Pp.*0,'bunsen');              % Bunsen solubility coefficients for N2, O2 and Ar in units of 'lair/(lwater.atmosphere)'

            P_O2c = C_O2c/B_O2; % Partial pressure of dissolved oxygen in units 'mbar'
            % start Henry's Law correction CLM on 02/24/17
            Vo2=33.1; Phi=1;c2=10.1325; R=82.05; %from Hamme email on : Wednesday, February 3, 2016 9:32 AM
            P_O2=(P_O2c).*exp(((Vo2.*Pp)./(Phi*c2*R*(Tt+273.15)))); % NB: INVERSE - going from potential to in situ
            C_O2 = B_O2 * P_O2c; % overkill, but keeps consistent numerically
            % end Henry's Law correction
            
            P_N2 = (TDG - P_O2 - P_H2O)./(1+((1-Xn2-Xo2)./Xn2));  % from here on, pAr is accounted for as a separate gas, and reduces GT so that pN2 comes or correctly and can then be compared to tropospheric N2 saturation.
            % start Henry's Law correction CLM on 12/22/16
            Vn2=33.1; Phi=1;c2=10.1325; R=82.05; %from Hamme email on : Wednesday, February 3, 2016 9:32 AM
            P_N2c=(P_N2).*exp((-(Vn2.*Pp)./(Phi*c2*R*(Tt+273.15))));
            C_N2 = B_N2 * P_N2c; 
            % end Henry's Law correction
           
            P_O2c_Sat = (1013.25-P_H2O)*Xo2; % Partial pressure of O2 in 1 standard atmosphere of moist air
            P_N2c_Sat = (1013.25-P_H2O)*Xn2; % Partial pressure of N2 in 1 standard atmosphere of moist air; correct, there is NO need to account for argon since its been taken care of already!!!
            P_Arc_Sat = (1013.25-P_H2O)*Xar; % Partial pressure of O2 in 1 standard atmosphere of moist air
                
            GTD_calc_O2sat(i) = 100*P_O2c/P_O2c_Sat; % Percent saturation level of O2 wrt an overlying atmosphere of moist air at 1013.25 mbar total pressure, assuming trace gases are all argon
            GTD_calc_N2sat(i) = 100*P_N2c/P_N2c_Sat; % Percent saturation level of N2 wrt an overlying atmosphere of moist air at 1013.25 mbar total pressure, assuming trace gases are all argon

            P_Arc = P_Arc_Sat*P_N2c/P_N2c_Sat; % same saturation level as N2 
            % start Henry's Law correction CLM on 12/22/16
            Var=33.1; Phi=1;c2=10.1325; R=82.05; %from Hamme email on : Wednesday, February 3, 2016 9:32 AM
            P_Ar=(P_Arc).*exp(((Var.*Pp)./(Phi*c2*R*(Tt+273.15))));
            C_Ar = B_Ar * P_Arc; 
            % end Henry's Law correction

            GTD_calc_Arsat(i) = 100*P_Arc/P_Arc_Sat; % Percent saturation level of N2 wrt an overlying atmosphere of moist air at 1013.25 mbar total pressure, assuming trace gases are all argon

                                                          %!!! CLM: 10/13/07 potential density
            GTD_calc_O2_umol(i)=C_O2./22.3916e-3./Dens; 
            GTD_calc_N2_umol(i)=C_N2./22.403e-3./Dens;
            GTD_calc_Ar_umol(i)=C_Ar./22.3916e-3./Dens;
            
            GTD_calc_O2star_umol(i)=GTD_calc_O2_umol(i)*100./GTD_calc_O2sat(i);
            GTD_calc_N2star_umol(i)=GTD_calc_N2_umol(i)*100./GTD_calc_N2sat(i);
            GTD_calc_Arstar_umol(i)=GTD_calc_Ar_umol(i)*100./GTD_calc_Arsat(i);
            
            GTD_calc_pN2(i)=P_N2;
            GTD_calc_pN2c(i)=P_N2c;
            GTD_calc_pO2(i)=P_O2;
            GTD_calc_pO2c(i)=P_O2c;
            GTD_calc_pAr(i)=P_Ar;
            GTD_calc_pArc(i)=P_Arc;
            
            GTD_calc_TDG(i)=P_N2 + P_O2 + P_Ar + P_H2O;
            GTD_calc_TDGc(i)=P_N2c + P_O2c + P_Arc + P_H2O;
            GTD_calc_dryTDG(i)=P_N2 + P_O2 + P_Ar;
            GTD_calc_dryTDGc(i)=P_N2c + P_O2c + P_Arc;
            
            GTD_calc_P(i)=Pp;
            GTD_calc_T(i)=Tt;
            GTD_calc_S(i)=Ss;
            GTD_calc_Mtime(i)=MT;
end
clear B_Ar B_N2 B_O2 C_N2 C_O2 Dens HH Lat Lon P_H2o P_N2 P_N2_Sat P_N2c P_O2 P_O2_Sat Phi Pp P_H20 c2 a k Vo2 TTL TGD Ss TDG Tt Ut R SS TDG Tt Dr DOFIG DIR P_H2O DR
    
figure(123);clf;hold on;grid on;

title('Hood Canal on 17-Aug-2017: O_2 (red)  N_2 (blue)');
plot([100 100],[0 25],'k--');
y=GTD_calc_P;
x=GTD_calc_N2sat;
%plot(x,y,'.','color',[0 0 0.5])
set(gca,'ydir','reverse');
edges = 0:1:25; % 26 edges, 25 bins
i=discretize(y,edges);
for j=1:max(size(edges))-1,
    b(j)=j-0.5;
    g=find(i==j);
    n(j)=nanmean(x(g));
    e(j)=nanstd(x(g));
end
plot(n,b,'o','markersize',6,'markerfacecolor',[0 0 1]);
errorbar(n,b,e,'o','horizontal','color',[0 0 1]);
y=GTD_calc_P;
x=GTD_calc_O2sat;
%plot(x,y,'.','color',[0.5 0 0])
set(gca,'ydir','reverse');
edges = 0:1:25; % 26 edges, 25 bins
i=discretize(y,edges);
for j=1:max(size(edges))-1,
    b(j)=j-0.5;
    g=find(i==j);
    n(j)=nanmean(x(g));
    e(j)=nanstd(x(g));
end
plot(n,b,'o','markersize',6,'markerfacecolor',[1 0 0]);
errorbar(n,b,e,'o','horizontal','color',[1 0 0]);
ylabel('P (dbar)');
xlabel('Gas Saturation Level (%)');

% export for archiving with NOAA

delete('Deconvolved Data for 08_17_2017.txt');
fileID = fopen('Deconvolved Data for 08_17_2017.txt','w');
fprintf(fileID, 'Date   Time    P   T   S   GT  N2  O2  N2star  O2star\n\n');
fprintf(fileID, '[day-month-year]   [hh:mm:ss]  [dbar]  [oC]    [psu]   [mbar]  [umol/kg]   [umol/kg]   [umol/kg]   [umol/kg] \n\n');
for i=1:max(size(GTD_calc_Mtime))
    i
    fprintf(fileID,'%s %f %f %f %f %f %f %f %f\n',datestr(GTD_calc_Mtime(i)),GTD_calc_P(i),GTD_calc_T(i),GTD_calc_S(i),GTD_calc_TDG(i),GTD_calc_N2_umol(i),GTD_calc_O2_umol(i),GTD_calc_N2star_umol(i),GTD_calc_O2star_umol(i));
end
fclose(fileID);
figure(99);clf;hold on;
plot(GTD.MT,GTD.Pgtd,'.b');
plot(GTD_calc_Mtime,GTD_calc_TDG,'or');
datetick


