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
    vtilde(i,:) = omega(i,:)-v;
end
m = ones(1,length(omega)); % Using constant priors, this has to be updated

%2)
k = 0;
lambda = zeros(1,length(omega(1,:)));
epsilon = 10e-5;

while true
    %3)
    
    %Symbolic version of the F function
    [H,~] = hessian(@F,lambda);
    [g,~,~] = gradest(@F,lambda);
    
    
    %4)
    deltalambda = -( (H^(-1)) *g.' );
    
    %5)
    alpha = 1; %for the moment
    lambda = lambda + (alpha * deltalambda.');
    
    %6)
    %The algorithm checks for the next g = f', but I think it's not
    %necessary.
    %[g,~] = derivest(@F,lambda)
    convcheck = norm(g); %check for convergence
    if (convcheck <= epsilon)
        break;
    else
        k = k+1;
    end
end

%disp(strcat('Iterations performed: ',int2str(k+1)))

%7)
for i = 1:length(omega)
    b(i) = Zi(lambda,i) / Z(lambda);
end

% Helper functions 

    % Z is the partition function
    function [zi] = Zi(lambda,j)
        zi = m(j)*exp(dot(-lambda,vtilde(j,:)));
    end

    function [z] = Z(lambda)
        z = 0;
        for j = 1:length(omega)
            z = z+ Zi(lambda,j);
        end
    end

    % F 
    function [f] = F(lambda)
        f = log(Z(lambda));
    end


end

