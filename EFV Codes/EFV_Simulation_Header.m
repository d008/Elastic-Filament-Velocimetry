clc
%% Create Common Wire Structs

%Common Fluids
water = struct('is','fluid','rho',1000,'mu',9E-4);
ethyleneglyocl = struct('is','fluid','rho',1113,'mu',1.61E-2);
air = struct('is','fluid','rho',1.2,'mu',1.8E-5);

%Common Solids
pt = struct('is','solid','rho',20e3,'E',168E9,'TCR',3.9E-3,'rho_R',10.6E-8);
ni = struct('is','solid','rho',8.9e3,'E',200E9,'TCR',6E-3,'rho_R',7E-8);
constantan = struct('is','solid','rho',8.9e3,'E',162E9,'TCR',8E-6,'rho_R',5e-7);
th =150e-9;

%Common Wire Geometries
for l = [60,375,535,750]
    for w = [2,4,8]
        %Geomtry w/o triangles
        v = genvarname(sprintf('r%dx%d',l,w));
        ans  =struct('L',l/(10^6),'L0',l/(10^6),'th',th,'w',(w-1.5)/(10^6),...
            'A',th*(w-1.5)/(10^6),'I',th.^3*(w-1.5)/(10^6)/12);
        eval([ v '=ans;']);
        %Geomtry w triangles
        v = genvarname(sprintf('efv%dx%d',l,w));
        ans  =struct('L',l/(10^6),'L0',l*1.12/(10^6),'th',150e-9,'w',(w-1.5)/(10^6),...
            'A',th*(w-1.5)/(10^6),'I',th.^3*(w-1.5)/(10^6)/12);
        eval([ v '=ans;']);
        %Geomtry w stubs
        v = genvarname(sprintf('n%dx%d',l,w));
        ans  =struct('L',2.6666*l/(10^6),'L0',l/(10^6),'th',150e-9,'w',(w)/(10^6),...
            'A',th*(w)/(10^6),'I',th.^3*(w)/(10^6)/12);
        eval([ v '=ans;']);
    end
end

% Drag Coefficient Power fit- Viscous
cdV = @(Re) (1.18.*(Re)/2+6.8.*(Re).^(0.11)/2+1.96.*(Re).^(0.5)/2-...
    0.0004*(Re).^2./(1+3.64E-7*(Re).^2)/2); 
% Drag Coefficient Power fit- Intertial
cdI = @(Re) cdV(Re)*2/Re;


