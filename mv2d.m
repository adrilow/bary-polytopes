function [phi, dphi] = mv2d(v,x)
%
% Evaluate MV basis functions and their gradients
% in a convex polygon
%
% Inputs:
% v : [x1 y1; x2 y2; ...; xn yn], the vertices, ccw
% x : [x(1) x(2)], point inside polygon
% Outputs:
% phi : basis functions = [phi_1; ...; phi_n]
% dphi : gradients of basis functions = [dphi_1; ...; dphi_n]

n = size(v,1);
w = zeros(n,1);
R = zeros(n,2);
phi = zeros(n,1);
dphi = zeros(n,2);

r = zeros(n,1);
e = zeros(n,2);

for i = 1:n
	d = v(i,:) - x;
	r(i) = norm(d);
	e(i,:) = d / r(i);
end

mod1 = inline('mod(i-1,n)+1','i','n');

t = zeros(n,1);
c2 = zeros(n,2);

for i = 1:n
	ip1 = mod1(i+1,n);
	cosa = dot(e(i,:),e(ip1,:));
	sina = det([e(i,:);e(ip1,:)]);
	t(i) = (1-cosa)/sina;
	c = e(i,:) / r(i) - e(ip1,:) / r(ip1);
	c2(i,:) = [-c(2) c(1)] / sina;
end

for i = 1:n
	im1 = mod1(i-1,n);
	w(i) = (t(im1) + t(i))/r(i);
	R(i,:) = (t(im1)*c2(im1,:)+t(i)*c2(i,:))/...
	(t(im1)+t(i))+e(i,:)/r(i);
end

phi = w/sum(w);

phiR = phi' * R;
for k = 1:2
	dphi(:,k) = phi .* (R(:,k) - phiR(:,k));
end