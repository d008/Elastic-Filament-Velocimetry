%% Execute Header
function [trecord, eps] = EFVMain()
EFV_Simulation_Header
%% Velocity
fluid = air;      %set fluid
material = constantan;      %set wire material
wire = r750x4;    %set the wire geometry
U = 10;            %set the velocity

%SET flexural part: 1 for include, 0 to exclude
flex = true;
record = false;
eps0=10^(-4);       %set pretension
N = 30;             %

L = wire.L; L0 = wire.L0;th =wire.th;w= wire.w;A=wire.A;I=wire.I;
rho_s = material.rho;E=material.E;
mu = fluid.mu; rho_f = fluid.rho;
%% Domain Parameters
%Make Grid
[x,D,aa,ind] = chebyGridMaker(N,flex);

%Further grid consideration
D4 = D^4;
D2 = D^2;

%Data storage and allocation
u = 0*x;v = x*0;
urecord = u;trecord =0;eps = 0;

%Width Profile
W = x; W(abs(x)<=L0/L) = 1;W(abs(x)>L0/L) = (abs(x(abs(x)>L0/L)*L/L0)-1)*10+1;

%% Approximate Steady State Solutions
Cd = cdV(w*U*rho_f/mu);     %Steady state Cd
q = Cd*U*mu;            %Load per unit span
Cm = 2;

Hu = 3.*q.*L./(E*A);
Q = @(HU, SP0) (HU+sqrt(HU.^2-SP0.^3)).^(1/3);
DELT = @(Hu,sp0) 2.^(2/3)./8.*L.*(sp0./Q(Hu,sp0)+Q(Hu,sp0));
EPS = @(Hu,sp0) 8/3.*DELT(Hu,sp0).^2./L^2 ;

%Scaling Parameters of deflection and timeclear aasda
delta= (Cd*(L/2)^4*U*mu/(E*A))^(1/3);
time = sqrt(rho_s*(L/2)^4/(E*delta^2));

%Scaling Parameters of individual terms
%Elastic-Axial Force
R1 = E*delta^2*time^2/((rho_s+rho_f*Cm)*(L/2)^4);
%Flexural Rigidity
R2 = E*I*time^2/((rho_s+rho_f*Cm)*A*(L/2)^4) * flex;
%Forcings
R3 = Cd*mu*U*time^2/((rho_s+rho_f*Cm)*A*delta);
%Damping
R4 = delta/(U*time);
%Added Mass
R5 = rho_s*Cm*U*time/((rho_s+rho_f*Cm)*delta)*0;

%%
%Linearized Parameters
omega0 = 8*(delta)/L^2*sqrt(E/(3*rho_s));
zeta = sqrt(3/(E*rho_s))*Cd*L^2*mu/(16*A*delta);
freqResp = (4/(omega0*zeta))^(-1);

T=2/zeta;
dt =(1-cos(1/N*pi));

T0 = eps0*((L/2)^2/delta^2); %Tension is negative

%Computational Functions
sigma=@(t,u,v) ((-trapz(x*L/2,sqrt(1+(D*u*delta/L*2).^2))./(L)-1)*(L/2)^2/delta^2)-T0;
V =@(t) t*0+1;%(-cos(t)+1)/2+(-cos(t/3)+1)/2+1;%(square((t+dt)/10)+1)/2%
fu = @(t,u,v) v;
fv= @(t,u,v) (-R2*D2*(D2*(u).*W)+...
    sigma(t,u,v)*R1*D2*u./W+...
    R3*(V(t)-R4*v).*cdV(abs((V(t)-v).*U.*rho_f.*W.*w./mu))./Cd)+...
    R5*(V(t+dt)-V(t-dt))/(dt*2);

%Set the inital deflection
delta2 = real(DELT(0,eps0*2*4^(1/3)));
u = delta2/(delta)*(1-(x).^2);
urecord = u;

eps = (sigma(0,u,v))/((L/2)^2/delta^2);
%% RK4 Iteration
ax = gca;
ax.NextPlot = 'replaceChildren';
h = waitbar(0,'Please wait...');
tic
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


model = [real(DELT(Hu,eps0*2*4^(1/3))-DELT(0,eps0*2*4^(1/3)))*10^6,real(EPS(Hu,eps0*2*4^(1/3))-EPS(0,eps0*2*4^(1/3)))*10^4]
simluation = [(max(u)*delta-delta2)*10^6,(eps(end)-eps(1))*10^4]
toc

%% Post Process
if record
    vidObj = VideoWriter('temp.avi');
    open(vidObj);
end
test = real(EPS(cdV(w.*U.*V(trecord./time).*rho_f./mu).*U.*V(trecord./time)*mu*L*3/(E*A),eps0*2*4^(1/3))-...
    EPS(0,eps0*2*4^(1/3)));
skip =10;
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
    plot(trecord,test,'-')
    hold off
    drawnow
    pause(0.001)
    if record
        currFrame = getframe(gcf);
        writeVideo(vidObj,currFrame);
    end
end
if record
    close(vidObj);
end