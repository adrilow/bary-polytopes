function [ b ] = maxEntrCoords( omega, v )
% maxEntrCoords Calculates maximum entropy coordinates for a given Polytope
%
% Uses Algorithm by Hormann & Sukumar (2008), currently with constant
% priors, so it ressembles the coordinates by Arroyo & Ortiz (2006), which
% satisfy partition of unity and linear precission, but the lagrange
% constraint only for strictly convex polytopes.
%
% Inputs:
% omega : [v1, ..., vn] (vi [w1, ..., wd] d : dim of omega) vertices of
% the polytope omega.
% v : [a1, ..., ad] point in omega.
%
% Outputs:
% b : [b1v, ..., bnv] barycentric coordinates for v in omega.

b = zeros(length(omega),length(omega(1)));



%1)
vtilde = arrayfun(@(vi) vi-v, omega);
m = ones(length(omega)); % Using constant priors, this has to be updated

%2)
k = 0;
lambda = 0;
epsilon = 10e-10;

while True
    %3)
    f,g,H = brownfgh(lambda);
    
    %4)
    deltalambda = -(1/H) * g;
    
    %5)
    alpha = 1; %for the moment
    lambda = lambda + (alpha * deltalambda);
    
    %6)
    if (norm(gradient(F(lambda))) <= epsilon)
        break;
    else
        k = k+1
    end
end

%7)
for i = 1:length(omega)
    b(i) = Zi(lambda,i) / Z(lambda);
end

% Helper functions

    % Z is the partition function
    function [zi] = Zi(lambda,j)
        zi = m(i)*exp((-lambda)*vtilde(j));
    end

    function [z] = Z(lambda)
        z = 0;
        for j = 1:length(omega)
            z = z + zi(lambda,j);
        end
    end

    % F 
    function [f] = F(lambda)
        f = ln(Z(lambda));
    end


end

