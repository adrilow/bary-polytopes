function [ b, K ] = maxEntrCoordsSym( omega, v )
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

%Retrieve or generate Gradient and Hessian of F for the current dim and
%number of points of the polytope.
gradStr = sprintf('Grad_%d_dim_%d_n',length(omega(1,:)),length(omega));
hesseStr = sprintf('Hesse_%d_dim_%d_n',length(omega(1,:)),length(omega));
if(exist(gradStr,'file')==0 || exist(hesseStr,'file')==0)
    generateGradHess(length(omega(1,:)),length(omega));
end
g = str2func(gradStr);
H = str2func(hesseStr);
selector = ones(1,length(omega)); %this will change, when n_max>n

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

    [gNum,HNum,~] = gradandhessian(@F,lambda)
    g_k = g(selector,lambda.',omega.')
    H_k = H(selector,lambda.',omega.')
    
    %4)
    deltalambda = -( ((H_k)^-1) * g_k );
    
    %5)
    alpha = 1; %for the moment
    lambda = lambda + (alpha * deltalambda.');
    
    %6)
    convcheck = norm(g_k); %check for convergence
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
K = k;

% Helper functions 

    % Z is the partition function
    function [zi] = Zi(lambda,j)
        zi = m(j)*exp(-(lambda*vtilde(j,:).'));
    end

    function [z] = Z(lambda)
        %It calls Z_i on every vertex
        res = arrayfun(@(x) Zi(lambda,x),1:length(omega));
        %And returns the result
        z = sum(res);

    end

    % F 
    function [f] = F(lambda)
        f = log(Z(lambda));
    end



end
