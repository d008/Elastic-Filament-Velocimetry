

L = 750E-6;             %Length of wire in m
L0 = 750E-6;            %Length of rectangle in m
th = 150e-9;             %Thickness of wire in m
w = 6.5E-6;              % Width of wire in m
rho =20e3;             %Density of Platinum in kg/m^3
V = @(t) 1/2*sin(t/20)+1/2
air = 1
if(air)
    rhof = 1.2;             %density of fluid
    mu = 1.8E-5;          % Viscosity of water in Pa sec
    U = 2;                 %Fluid Velocity in m/s
else
    rhof = 1000;             %density of fluid
    mu = 1E-3;          % Viscosity of water in Pa sec
    U = 0.4;                  %Fluid Velocity in m/s
end
A= th*w;                 %Cross Sectional Area of the wire

E = 168e9;              %Youngs Modulus in GPA
I = th.^3.*w./12;        %Second Moment
ReNu = rhof*U*w./mu     %Reynolds Number

%dragcoef2 = @(Re) 4*pi/(log(2*L/w)-0.5); %DRAG COEFFICIENT FOR A CYLINDER
dragcoef = @(Re) 0.5.*Re + 5.*Re.^(1/3); %DRAG COEFFICIENT FOR A CYLINDER
dragcoef = @(Re) (1.18+6.8./(Re).^(0.89)+1.96./(Re).^(0.5)-0.0004*(Re)./(1+3.64E-7*(Re).^2)).*(Re)/2; %DRAG COEFFICIENT FOR A CYLINDER
dragcoefN = @(Re)dragcoef(Re)*2/Re;
%dragcoef = @(Re) 10
%dragcoef2 = @(Re) 4*pi/(log(2*L/w)-0.5); %DRAG COEFFICIENT FOR A CYLINDER

Cd = dragcoef(ReNu)
q = Cd*U*mu;            %Load per unit span

Hu = 3.*q.*L./(E*A);
Q = @(HU, SP0) (HU+sqrt(HU.^2-SP0.^3)).^(1/3);
DELT = @(Hu,sp0) 2.^(2/3)./8.*L.*(sp0./Q(Hu,sp0)+Q(Hu,sp0));
EPS = @(Hu,sp0) 8/3.*DELT(Hu,sp0).^2./L^2 ;

%Scaling Parameters of deflection and timeclear aasda
delta= (Cd*(L/2)^4*U*mu/(E*A))^(1/3)
time = sqrt(rho*(L/2)^4/(E*delta^2));

%SET flexural part: 1 for include, 0 to exclude
flex = true;

%Scaling Parameters of individual terms

%Elastic-Axial Force
R1 = E*delta^2*time^2/(rho*(L/2)^4);
%Flexural Rigidity
R2 = E*I*time^2/(rho*A*(L/2)^4) * flex;
%Forcings
R3 = Cd*mu*U*time^2/(rho*A*delta);
%Damping
R4 = delta/(U*time);
%Frequency
omeger = 1/time;
%Damping term 2
damp =R3*R4;

%Linearized Parameters
omega0 = 8*(delta)/L^2*sqrt(E/(3*rho));
zeta = sqrt(3/(E*rho))*Cd*L^2*mu/(16*A*delta);
freqResp = (4/(omega0*zeta))^(-1);

%Domain
N =20;
x = 0:N;
x = cos(x/(N)*pi)';
W = x; W(abs(x)<=L0/L) = 1;W(abs(x)>L0/L) = (abs(x(abs(x)>L0/L)*L/L0)-1)*10+1;

%Domain Considerations
u = 0*x;
v = x*0;
urecord = u;
trecord =0;eps = 0;

%Derivative Matrices
D = zeros(N+1,N+1);
c = @(j) 1 + (j==1) +(j==(N+1));
for j = 1:N+1
    for k = 1:N+1
        if(j~=k)
            D(j,k) = c(j)*(-1)^(j-1+k-1)/(c(k)*(x(j)-x(k)));
        elseif(j==k & j~=1 & j~= N+1)
            D(j,k) = -x(j)/(2*(1-x(j)^2));
        elseif(j==1 & k==1)
            D(j,k) = (2*(N^2)+1)/6;
        elseif(j==N+1 & k==N+1)
            D(j,k) = -(2*(N^2)+1)/6;
        end
    end
end
D;
D4 = D^4;
D2 = D^2;
Dl = D(1,:)';
aa = -Dl(3:N-1)'./(Dl(2)+Dl(N));


%Prediction from low order model
%LowOrderModel = [1/(2*3^(1/3))*(q*L/(E*A))^(2/3)*10^4,(3/64*q*L^4/(E*A))^(1/3)*10^6]

%% Rk4
T=1/(omega0*zeta)*omeger;
eps0=10^(-36);
T0 = eps0*((L/2)^2/delta^2); %Tension is negative

sigma=@(t,u,v) ((-trapz(x*L/2,sqrt(1+(D*u*delta/L*2).^2))./(L)-1)*(L/2)^2/delta^2)-T0;
fu = @(t,u,v) v;
fv= @(t,u,v) (-R2*D2*(D2*(u).*W)+sigma(t,u,v)*R1*D2*u./W+R3*(1-R4*v).*dragcoef(abs((1-v).*U.*rhof.*W.*w./mu))./Cd);


delta2 = real(DELT(0,eps0*2*4^(1/3)))
u = delta2/(delta)*(1-(x).^2);
urecord = u;

eps = (sigma(0,u,v))/((L/2)^2/delta^2);


LowOrderModel =real([(EPS(Hu,eps0*2*4^(1/3))-EPS(0,eps0*2*4^(1/3)))*10^4,...
    (DELT(Hu,eps0*2*4^(1/3))-DELT(0,eps0*2*4^(1/3)))*10^6])

dt =2*(1-cos(1/N*pi))/4;
t=0;
figure(1)
tic

%4 order



ax = gca;
ax.NextPlot = 'replaceChildren';
h = waitbar(0,'Please wait...');
%%
for t=0:dt:T
    if(isnan(eps(end)))
        break
    end
    waitbar(t/T,h,sprintf('Strain = %0.2f, t = %0.2f',(eps(end)-eps(1))*10^4,t));
    u0 = u;
    v0 = v;
    %K1
    k1u = fu(t,u,v);
    k1v = fv(t,u,v);
    u(ind) = u0(ind)+k1u(ind)*dt/2;
    v(ind) = v0(ind)+k1v(ind)*dt/2;
    %u(2) = u(3)/4;
    %u(N) = u(N-1)/4;
    if(flex ==true)
        u(2) = aa*u(3:N-1);
        u(N) =u(2);
    end
    %k2
    k2u = fu(t,u,v);
    k2v = fv(t,u,v);
    u(ind) = u0(ind)+k2u(ind)*dt/2;
    v(ind) = v0(ind)+k2v(ind)*dt/2;
    if(flex ==true)
        u(2) = aa*u(3:N-1);
        u(N) =u(2);
    end
    
    %k3
    k3v = fv(t,u,v);
    k3u = fu(t,u,v);
    u(ind) = u0(ind)+k3u(ind)*dt;
    v(ind) = v0(ind)+k3v(ind)*dt;
    if(flex ==true)
        u(2) = aa*u(3:N-1);
        u(N) =u(2);
    end
    
    %k4
    k4u = fu(t,u,v);
    k4v = fv(t,u,v);
    u(ind) = u0(ind)+dt/6*(k1u(ind)+k2u(ind)*2+k3u(ind)*2+k4u(ind));
    v(ind) = v0(ind)+dt/6*(k1v(ind)+k2v(ind)*2+k3v(ind)*2+k4v(ind));
    if(flex ==true)
        u(2) = aa*u(3:N-1);
        u(N) =u(2);
    end
    eps=[eps;(sigma(t,u,v))/((L/2)^2/delta^2)];
    [t,eps(end)*10^4,8/3*delta^2/(L)^2*max(u).^2*10^4];
    urecord = [urecord,u];
    trecord = [trecord,t*time];
    sigma(t,u,v);
    
end
close(h)
urecord = urecord*delta;
x = x*L/2;

%Prediction from simulations order model
Simulation = [eps(end)*10^4,max(u*delta)*10^6]
toc

plot(x,u*delta,'-')
hold on
plot(x,sqrt((0.3467*(q*L/(E*A))^(2/3))*3/8*L^2).*(1-(2*x./L).^2))
hold off
%%
espComp = real([(eps(end))-(eps(1)),EPS(Hu,eps0*2*4^(1/3))-EPS(0,eps0*2*4^(1/3))]*10^4)
deltComp = real([max(u)*delta-delta2,DELT(Hu,eps0*2*4^(1/3))-DELT(0,eps0*2*4^(1/3))]*10^6)

vidObj = VideoWriter('temp.avi');
open(vidObj);
skip =100;
for i = 1:length(trecord)/skip-1
    figure(1)
    subplot(2,1,1)
    %plot(x*10^6,(delta2)*(1-(x/L*2).^2)*10^6,'--')
    plot(x*10^6,urecord(:,i*skip)*10^6,'-*')
    hold on
    plot(x*10^6,sqrt((0.3467*(q*L/(E*A))^(2/3))*3/8*L^2).*(1-(2*x./L).^2)*10^6,':')
    hold off
    xlabel('x (microns)')
    ylabel('\delta (microns)')
    title(strcat('t=',int2str(trecord(i*skip)*10^6),'microsec'))
    %axis([-L/2*10^6 L/2*10^6 0 1.5*delta*10^6*max(u)])
    set(gca,'FontSize',24)
    subplot(2,1,2)
    plot(trecord,eps)
    ylabel('Strain')
    xlabel('Time (sec)')
    set(gca,'FontSize',24)
    hold on
    plot(trecord(i*skip),eps(i*skip),'o')
    hold off
    drawnow
    pause(0.001)
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
end
close(vidObj);
