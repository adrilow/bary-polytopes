function [] = profile_max_entr_coords() 

%profile on -history
for dim= 2:2% 2:20 
    results = zeros(2,(100*dim-(dim+1)));
    it = 1;
    for points = dim+1:dim*5 %really dim*100
        disp(strcat('Initializing test with dim:',int2str(dim)))
        runtime_average = 0;
        tic;
        
        runs = 2; %for testing, this number has to be bigger
        for i = 1:runs %runs, for computing average
    
            %disp(strcat(int2str(i),'-th Iteration'))
    
            setindim = zeros(points,dim); % points in the 100^d Hypercube with an edge in [0,...,0]
            for d = 1:points
                setindim(d,:) = 100.*rand(1,dim);
            end
            P = Polyhedron(setindim); %generate Polyhedron from the given points
            U = PolyUnion(P); %The converHull method is defined on Polyunion
            PC = U.convexHull; %This will be the convex polytope on which the maxEntrCoords will be calculated
            PC.minVRep(); %We want its (non-redundant) vertices -> strictly convex
            omega = PC.V; %The Vertices of the Polyhedron
    
            %disp('Convex Polytope tested:')
            %disp(omega)
     
            strict = false;
    
            while ~strict
                vSol = PC.interiorPoint; % Determine an arbitrary interior point
                strict = vSol.isStrict; % Which is strictly in the relative interior
            end
            
            v = vSol.x; % Take the coordinates of the point.
    
            %disp('On the point:')
            %disp(v)
    
            tstart = tic;
            b = maxEntrCoords(omega,v); %Calculate b1,...,bn for the Point v in the Polytope omega.
            runtime_average = runtime_average + toc(tstart)/runs;
    
            %disp(fprintf(strcat('b_i: \n',mat2str(b))))
    
        end
        results(1,it) = points;
        results(2,it) = runtime_average;
        it = it+1;
    end
    %csvwrite(strcat(int2str(dim),'-Dimensions.csv'),results);
end
%p = profile('info')
%profsave(p);
%profile off
%DataField = fieldnames(p);
%jsonencode(struct2cell(p.(DataField{1})));
end