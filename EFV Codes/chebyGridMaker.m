function [x,D,aa,ind] = chebyGridMaker(N,flex)
%chebyGridMaker - Generates grid with N grid points w/ chebyshev spacing
%   Takes a value N an returns the grid and derivatives for the chebyshev
%   spaced grid.
%   Returns - grid, x
%           - width profile, W
%           - first Derivative, D
%           - Boundary derivative -aa

%Make a grid with cosine spacing
x = 0:N;
x = cos(x/(N)*pi)';

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

Dl = D(1,:)';
aa = -Dl(3:N-1)'./(Dl(2)+Dl(N));

ind = 2:N;
if(flex ==true)
    ind = 3:N-1;
end

end

