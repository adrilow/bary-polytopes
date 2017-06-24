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

if(length(v(1,:))==1) %If the point is given as column vector
    v = v.';          %Transform to line vector
end
b = zeros(length(omega),length(omega(1)));



%1)
vtilde = zeros(length(omega),length(omega(1,:)));
for i = 1:length(omega)
    vtildei = vtilde(i,:);
    omegai = omega(i,:);
    
    vtilde(i,:) = omega(i,:)-v;
end
m = ones(1,length(omega)); % Using constant priors, this has to be updated

%2)
k = 0;
lambda = 0;
epsilon = 10e-10;

while true
    %3)
    %In the paper they talk about the gradient and the Hessian, but as
    %these are scalar functions, these are just the first and second
    %derivative at lambda.
    [g,~] = derivest(@F,lambda); 
    [H,~] = derivest(@F,lambda,'DerivativeOrder',2);
    
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
        vtildej = vtilde(j,:);
        zi = m(j)*exp((-lambda)*norm(vtildej));
    end

    function [z] = Z(lambda)
        z = 0;
        for j = 1:length(omega)
            z = z + Zi(lambda,j);
        end
    end

    % F 
    function [f] = F(lambda)
        zlambda = Z(lambda);
        f = log(zlambda);
    end


end

