
%%
%%%%load('/Users/Matt/Dropbox/Elastic Filament Velocimetry/Mafu Numerical/ExpDataCompiled.mat')
%%%%load('/Users/Matt/Dropbox/Elastic Filament Velocimetry/EFV Confocal Data/ExperimentDataAfter.mat')

%load('C:\Users\Matt\Dropbox\Elastic Filament Velocimetry\Experimental Data\For real data now\750x8both.mat')
L = 840E-6;             %Length of wire in m
L0 = 750E-6;            %Length of rectangle in m
th = 150e-9;             %Thickness of wire in m
w = 6.5E-6;              % Width of wire in m
%th =w;
rho =20e3;             %Density of Platinum in kg/m^3
air = 1
GF=2.4           % Gain Factor
if(air)
    rhof = 1.2;             %density of fluid
    mu = 1.8E-5;          % Viscosity of water in Pa sec
    U = 30;                  %Fluid Velocity in m/s
else
    rhof = 1000;             %density of fluid
    mu = 1E-3;          % Viscosity of water in Pa sec
    U = 0.7;                  %Fluid Velocity in m/s
end

A= th*w;                 %Cross Sectional Area of the wire
us =logspace(-8,log10(U),1000);
E = 168e9;              %Youngs Modulus in GPA
I = th.^3.*w./12;        %Second Moment


ReNu = rhof*U*w./mu
shft =0
dragcoef = @(Re) (1.18+6.8./(Re+shft).^(0.89)+1.96./(Re+shft).^(0.5)-0.0004*(Re+shft)./(1+3.64E-7*(Re+shft).^2)).*(Re+shft)/2; %DRAG COEFFICIENT FOR A CYLINDER
%dragcoef = @(Re) 0.5.*Re + 5.*Re.^(1/3); %DRAG COEFFICIENT FOR A CYLINDER
%dragcoef = @(Re) 8; %DRAG COEFFICIENT FOR A CYLINDER


q = dragcoef(rhof*us*w./mu).*us.*mu;            %Load per unit span
sp0 = 0;
Hu = 3.*q.*L./(E*A);
Q = @(HU, SP0) (HU+sqrt(HU.^2-SP0.^3)).^(1/3);
DELT = @(Hu,sp0) 2.^(2/3)./8.*L.*(sp0./Q(Hu,sp0)+Q(Hu,sp0));
EPS = @(Hu,sp0) 8/3.*DELT(Hu,sp0).^2./L^2 ;



clf

%TENSION IS Nwgative
N = 1;
Nint = floor(32/N);
cs = parula;
ns = -linspace(3,3.6,N)

for i = 1:N
    sp = 10^(ns(i));
    bar = sqrt((sp).^3*32).^(1/3);
    %p = plot(Hu,(DELT(Hu,-10^(ns(i))*2*4^(1/3))-DELT(0,-10^(ns(i))*2*4^(1/3))),'-',Hu,(DELT(Hu,10^(ns(i))*2*4^(1/3))-DELT(0,10^(ns(i)))*2*4^(1/3)),'-.','LineWidth',2);
    p = plot(Hu,(EPS(Hu,-10^(ns(i))*2*4^(1/3)))-(EPS(0,-10^(ns(i))*2*4^(1/3))),'k--',Hu,(EPS(Hu,10^(ns(i))*2*4^(1/3)))-(EPS(0,10^(ns(i))*2*4^(1/3))),'k:','LineWidth',1.5);
end
close all
w25 =2.5E-6;
L375 = 526E-6;
L200 =200E-6;
ms =12


figure('Position',[ 700 300 800 600])
hold on
for i = 1:N
    sp = 10^(ns(i));
    bar = sqrt((sp).^3*32).^(1/3);
    %p = plot(Hu,(DELT(Hu,-10^(ns(i))*2*4^(1/3))-DELT(0,-10^(ns(i))*2*4^(1/3))),'-',Hu,(DELT(Hu,10^(ns(i))*2*4^(1/3))-DELT(0,10^(ns(i)))*2*4^(1/3)),'-.','LineWidth',2);
    p = plot(us,((EPS(Hu,10^(ns(i))*2*4^(1/3)))-(EPS(0,10^(ns(i))*2*4^(1/3))))*GF,'k--','LineWidth',1.5);
end
a = plot(us,EPS(Hu,0)*GF,'k-');
%d = plot(AnumT(:,1),8/3*((AnumT(:,2)-7)./840).^2-8/3*((12-7)/840).^2,'ko','MarkerSize',ms,'MarkerFaceColor',[0.8 0.8 0.8 ]);
%d2 = plot(AnumT(:,1),AnumT(:,3),'kd','MarkerSize',ms,'MarkerFaceColor',[0.8 0.8 0.8]);
if(air)
    d = errorbar(1.1*Vair,(dRRair+dRRair(end-1)),errair,'ks','MarkerSize',ms,'MarkerFaceColor',[0.8 0.8 0.8 ]);
    axis([0 30 -2E-4 10E-4])
    grid on
xlabel('x')
ylabel('y')
set(gca,'FontSize',20)

print(gcf,'-depsc','velair')
else
    d2 = errorbar(Vwater*2,(dRRwater-dRRwater(end-1)),errwater,'ks','MarkerSize',ms,'MarkerFaceColor',[0.3 0.3 0.3]);
    axis([0 0.6 -2E-4 10E-4])
    grid on
xlabel('x')
ylabel('y')
set(gca,'FontSize',20)

print(gcf,'-depsc','velwater')
end
hold off


%%
figure('Position',[ 100 300 800 600])
hold on
for i = 1:N
    sp = 10^(ns(i));
    bar = sqrt((sp).^3*32).^(1/3);
    %p = plot(Hu,(DELT(Hu,-10^(ns(i))*2*4^(1/3))-DELT(0,-10^(ns(i))*2*4^(1/3))),'-',Hu,(DELT(Hu,10^(ns(i))*2*4^(1/3))-DELT(0,10^(ns(i)))*2*4^(1/3)),'-.','LineWidth',2);
    p = plot(Hu,(EPS(Hu,10^(ns(i))*2*4^(1/3)))-(EPS(0,10^(ns(i))*2*4^(1/3))),'k--','LineWidth',1.5);
end

grid on
a = plot(Hu,EPS(Hu,0),'k-');
d = errorbar(real((1.1.*AnumT(:,1).*(1.8E-5).*L0.*dragcoef(1.2.*1.1.*AnumT(:,1).*w./(1.8E-5)).*3./(E.*A))),AnumT(:,3),AnumT(:,4),'kd','MarkerSize',ms,'MarkerFaceColor',[0.8 0.8 0.8]);
hold off
%xlabel('$\frac{C_D U \mu L}{E  w t}$','Interpreter','latex');
%ylabel('$\varepsilon$','Interpreter','latex');
xlabel('x')
ylabel('y')
set(gca,'FontSize',20)
axis([0 1E-4 -0.6E-4 3.5E-4])
%legend([a;p;d],{'Theory, \sigma^+_0 = 0 ','Theory w/ Pre-tension,\sigma^+_0 = -10^{-4}','Theory w/ Pre-deflection,\sigma^+_0 = 10^{-4}','750 \times 6.5 micron Confocal Data in Air'},'location','BestOutside')

print(gcf,'-depsc','confocaldata')
%%
figure('Position',[ 300 300 800 600])
hold on
for i = 1:N
    sp = 10^(ns(i));
    bar = sqrt((sp).^3*32).^(1/3);
    %p = plot(Hu,(DELT(Hu,-10^(ns(i))*2*4^(1/3))-DELT(0,-10^(ns(i))*2*4^(1/3))),'-',Hu,(DELT(Hu,10^(ns(i))*2*4^(1/3))-DELT(0,10^(ns(i)))*2*4^(1/3)),'-.','LineWidth',2);
    p = plot(Hu,(EPS(Hu,10^(ns(i))*2*4^(1/3)))-(EPS(0,10^(ns(i))*2*4^(1/3))),'k--','LineWidth',1.5);
end

a = plot(Hu,EPS(Hu,0),'k-');
b =errorbar((1.1*Vair.*muair.*L0.*dragcoef(1.2*1.1*Vair*w./muair).*3./(E.*A)),(dRRair+dRRair(end-1))/GF,errair/GF,'ks','MarkerSize',ms,'MarkerFaceColor',[0.8 0.8 0.8]);
c = errorbar((2*Vwater.*muwater.*L0.*dragcoef(2*1000*Vwater*w./muwater).*3./(E.*A)),(dRRwater-dRRwater(end-1))/GF,errwater/GF,'ks','MarkerSize',ms,'MarkerFaceColor',[0.3 0.3 0.3]);
%b =plot((1.1*Vair.*muair.*L0.*dragcoef(1.2*1.1*Vair*w./muair).*3./(E.*A)),(dRRair+dRRair(end-1))/GF,...
%    'ks','MarkerSize',ms,'MarkerFaceColor',[0.8 0.8 0.8]);
%c = plot((2*Vwater.*muwater.*L0.*dragcoef(2*1000*Vwater*w./muwater).*3./(E.*A)),(dRRwater-dRRwater(end-1))/GF,...
%    'ks','MarkerSize',ms,'MarkerFaceColor',[0.3 0.3 0.3]);
%d = errorbar(real((1.1.*AnumT(:,1).*(1.8E-5).*L0.*dragcoef(1.2.*1.1.*AnumT(:,1).*w./(1.8E-5)).*3./(E.*A))),AnumT(:,3),AnumT(:,4),'kd','MarkerSize',ms,'MarkerFaceColor',[0.8 0.8 0.8]);
%d2 = plot((1.1.*AnumT(:,1).*(1.8E-5).*L0.*dragcoef(1.2.*1.1.*AnumT(:,1).*w./(1.8E-5)).*3./(E.*A)),8/3*((AnumT(:,2)-7)./840).^2-8/3*((12-7)/840).^2,'ko','MarkerSize',19,'MarkerFaceColor',[0.8 0.8 0.8]);
%e = plot((VwaterT.*muwaterT.*L200.*dragcoef(1000*VwaterT*w25./muwaterT)).*3/(E*(w25*th)),(dRRwaterT-dRRwaterT(end-1))/GF,'ko','MarkerSize',19,'MarkerFaceColor',[1 0 0]);
hold off
xlabel('x')
ylabel('y')
%xlabel('$\frac{C_D U \mu L}{E  w t}$','Interpreter','latex');
%ylabel('$\varepsilon$','Interpreter','latex');
set(gca,'FontSize',20)
axis([0 1E-4 -0.6E-4 3.5E-4])
grid on
%{
legend([a;p;b;c],{'Theory, \sigma^+_0 = 0 ','Theory w/ Pre-tension,\sigma^+_0 = -10^{-4}',...
    'Theory w/ Pre-deflection,\sigma^+_0 = 10^{-4}','750 \times 6.5 micron Resistance Data in Air',...
    '750 \times 6.5 micron Resistance Data in Water',},'location','BestOutside')
%}

print(gcf,'-depsc','750by8both')
%%
figure('Position',[ 500 300 800 600])
hold on
for i = 1:N
    sp = 10^(ns(i));
    bar = sqrt((sp).^3*32).^(1/3);
    %p = plot(Hu,(DELT(Hu,-10^(ns(i))*2*4^(1/3))-DELT(0,-10^(ns(i))*2*4^(1/3))),'-',Hu,(DELT(Hu,10^(ns(i))*2*4^(1/3))-DELT(0,10^(ns(i)))*2*4^(1/3)),'-.','LineWidth',2);
    p = plot(Hu,(EPS(Hu,10^(ns(i))*2*4^(1/3)))-(EPS(0,10^(ns(i))*2*4^(1/3))),'k--','LineWidth',1.5);
end
a = plot(Hu,EPS(Hu,0),'k-');
b =plot((1.1*Vair.*muair.*L0.*dragcoef(1.2*1.1*Vair*w./muair).*3./(E.*A)),(dRRair+dRRair(end-1))/GF,'ks','MarkerSize',ms,'MarkerFaceColor',[0.8 0.8 0.8]);
c = plot((2*Vwater.*muwater.*L0.*dragcoef(2*1000*Vwater*w./muwater).*3./(E.*A)),(dRRwater-dRRwater(end-1))/GF,'ks','MarkerSize',ms,'MarkerFaceColor',[0.3 0.3 0.3]);
d = plot(real((1.1.*AnumT(:,1).*(1.8E-5).*L0.*dragcoef(1.2.*1.1.*AnumT(:,1).*w./(1.8E-5)).*3./(E.*A))),AnumT(:,3),'kd','MarkerSize',ms,'MarkerFaceColor',[0.8 0.8 0.8]);
f = plot((2*Vwater3754.*muwater3754*L375.*dragcoef(2*1000*Vwater3754*w25./muwater3754)*3)/(E*(w25*th)),(dRRwater3754-dRRwater3754(end-1))/GF,'kv','MarkerSize',ms,'MarkerFaceColor',[0.3 0.3 0.3]);
g = plot((2*Vwater3758.*muwater3758*L375.*dragcoef(2*1000*Vwater3758*w./muwater3758)*3)/(E*(w*th)),(dRRwater3758-dRRwater3758(end-1))/GF,'k^','MarkerSize',ms,'MarkerFaceColor',[0.3 0.3 0.3]);
%xlabel('$\frac{C_D U \mu L}{E  w t}$','Interpreter','latex');
%ylabel('$\varepsilon$','Interpreter','latex');
xlabel('x')
ylabel('y')
set(gca,'FontSize',20)
grid on

axis([0 1E-4 -0.6E-4 3.5E-4])
%{
legend([a;p;b;c;d;f;g],{'Theory, \sigma^+_0 = 0 ','Theory w/ Pre-tension,\sigma^+_0 = -10^{-4}',...
    'Theory w/ Pre-deflection,\sigma^+_0 = 10^{-4}','750 \times 6.5 micron Resistance Data in Air',...
    '750 \times 6.5 micron Resistance Data in Water', '750 \times 6.5 micron Confocal Data in Air'...
    '375 \times 2.5 micron Resistance Data in Water','375 \times 6.5 micron Resistance Data in Water'},'location','BestOutside')
%}
print(gcf,'-depsc','alldata')

%%

figure('Position',[ 700 300 1100 600])
hold on
for i = 1:N
    sp = 10^(ns(i));
    bar = sqrt((sp).^3*32).^(1/3);
    %p = plot(Hu,(DELT(Hu,-10^(ns(i))*2*4^(1/3))-DELT(0,-10^(ns(i))*2*4^(1/3))),'-',Hu,(DELT(Hu,10^(ns(i))*2*4^(1/3))-DELT(0,10^(ns(i)))*2*4^(1/3)),'-.','LineWidth',2);
    p = plot(us,(EPS(Hu,-10^(ns(i))*2*4^(1/3)))-(EPS(0,-10^(ns(i))*2*4^(1/3))),'k--',us,(EPS(Hu,10^(ns(i))*2*4^(1/3)))-(EPS(0,10^(ns(i))*2*4^(1/3))),'k:','LineWidth',1.5);
end
a = plot(us,EPS(Hu,0),'k-');
%d = plot(AnumT(:,1),8/3*((AnumT(:,2)-7)./840).^2-8/3*((12-7)/840).^2,'ko','MarkerSize',ms,'MarkerFaceColor',[0.8 0.8 0.8 ]);
%d2 = plot(AnumT(:,1),AnumT(:,3),'kd','MarkerSize',ms,'MarkerFaceColor',[0.8 0.8 0.8]);
d = plot(AnumT(:,1),8/3*((AnumT(:,2)-7)./840).^2-8/3*((12-7)/840).^2,'ko','MarkerSize',ms,'MarkerFaceColor',[0.8 0.8 0.8 ]);
d2 = plot(AnumT(:,1),AnumT(:,3),'kd','MarkerSize',ms,'MarkerFaceColor',[0.8 0.8 0.8]);
hold off
grid on
axis([0 25 -0.6E-4 3.5E-4])
set(gca,'FontSize',20)


