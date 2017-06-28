function [] = testMaxEntrCoord( dim )
%TestMaxEntrCoord tests the coordinates generated by maxEntrCoords() and
%proofs whether they are:
% Positive
% Linearly Precise (with tolerance of 10e-10)
% A partition of unity
% By generating random dim-Dimensional polytopes and testing points in the
% inside.
%
% Inputs:
% Dimension of the polyopes.

%Generate 1000 random Polytopes of dimension <= dim
disp(strcat('Initializing test with dim:',int2str(dim)))
for i = 1:1
    
    disp(strcat(int2str(i),'-th Iteration'))
    
    points = randi([dim+1,dim*10]); % number of points
    setindim = zeros(points,dim); % points in the 100^d Hypercube with an edge in [0,...,0]
    for d = 1:points
        setindim(d,:) = 100.*rand(1,dim);
    end
    P = Polyhedron(setindim); %generate Polyhedron from the given points
    U = PolyUnion(P); %The converHull method is defined on Polyunion
    PC = U.convexHull; %This will be the convex polytope on which the maxEntrCoords will be calculated
    PC.minVRep(); %We want its (non-redundant) vertices -> strictly convex
    omega = PC.V; %The Vertices of the Polyhedron
    
    disp('Convex Polytope tested:')
    disp(omega)
    
    
    strict = false;
    
    while ~strict
        vSol = PC.interiorPoint; % Determine an arbitrary interior point
        strict = vSol.isStrict; % Which is strictly in the relative interior
    end
    
    
    v = vSol.x; % Take the coordinates of the point.
    
    disp('On the point:')
    disp(v)
    
    b = maxEntrCoords(omega,v) %Calculate b1,...,bn for the Point v in the Polytope omega.
    
    %Test of Positivity
    disp('Testing positivity...')
    for j = 1:length(b)
        if(b(j) < 0)
            disp('Failed!')
            return
        end
    end
    disp('Passed')
    
    %Test of Linear Precission
    disp('Testing linear precission')
    vIs = zeros(1,dim);
    for j = 1:length(b)
        vIs = vIs + b(j)*omega(j,:);
    end
    
    dist = norm(vIs.'-v);
    disp(strcat('Should: ',mat2str(v.'),' Is: ',mat2str(vIs), ' Distance: ',mat2str(dist)))
    if(dist > 10e-10) %note that epsilon is only defined for part. of unity, but this way some error due to rounding et al. is allowed
        disp('Failed')
        return
    end
    disp('Passed')
    
    %Test of Partition of unity
    disp('Testing Partition of unity')
    total = 0;
    for j = 1:length(b)
        total = total + b(j);
    end
    disp(total)
    if(abs(1-total)>10e-10)
        disp('Failed')
        return
    end
    disp('Passed')
    disp('All tests passed')
    
end


end

