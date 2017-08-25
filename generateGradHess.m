function [] = generateGradHess (d,n_max)

    %d= 5; %Dimension of Space
    %n_max = 20; %Maximum vertices of Polytope

    V_TILDE = sym('v_tilde',[d,n_max]) %Will then be used to select n [1..n..1,0..n_max-n..0]
    syms n k

    lambda = sym('lambda',[d,1]); 
    a = sym('a',[1,n_max]);				%Switch, 1 for element of Vtilde, 0 for Dummy (can also be used as prior function)

    F=log(a*exp(-(transpose(V_TILDE)*lambda)));

    grad_F = transpose(jacobian(F,lambda));
    hesse_F = jacobian(grad_F,lambda);
    
    %Out of higher dimensional grads/hesses one can get all lower
    %dimensional versions for the same n_max by taking the first lowdim
    %(or lowdim x lowdim) elements of the grad/hess.
    
    matlabFunction( grad_F, 'File', sprintf('Grad_%d_dim_%d_n.m',d,n_max), 'vars',{a,lambda,V_TILDE});
    matlabFunction( hesse_F, 'File', sprintf('Hesse_%d_dim_%d_n.m',d,n_max), 'vars',{a,lambda,V_TILDE});
end