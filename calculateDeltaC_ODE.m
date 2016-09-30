function [ dC_vec ] = calculateDeltaC_ODE( x_in, D, Hlocal, bicg_maxiter, bicg_tolerance)
   
global    Bzero_GLOBAL;

    C_vec = x_in;

    C = transpose(reshape(C_vec, D, D));
    
    [ C, Eeigvals, L, R ] = calculateEigenvectors_inv(C);
    
    Rtrunc = R(2:end, :);
    Ltrunc = L(2:end, :);
    
    %Some useful contractions
 
    %%vl_C_Cbar %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vl_C_Cbar = bsxfun(@times, transpose(L(1,:)), C);
    vl_C_Cbar = transpose(vl_C_Cbar)*conj(C); 
    
    %%C_Cbar_vr %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    C_Cbar_vr = bsxfun(@times, C, R(1,:));
    C_Cbar_vr = C_Cbar_vr*C';
    
    %%%%%%%%% W1
    O1 = Hlocal{4}.*C_Cbar_vr;
    O2 = C*O1*C';
    OR = eye(D).*O2;
    
    O1 = eye(D).*C_Cbar_vr;
    O2 = C*O1*C';
    OR = OR + Hlocal{1}.*O2;
    
    O1 = Hlocal{3}.*C_Cbar_vr;
    O2 = C*O1*C';
    OR = OR + Hlocal{2}.*O2;
    
    O3 = bsxfun(@times, transpose(L(1,:)), C)*OR;
   
    O4 = O3.*(1./C);
    W1 = conj(Rtrunc)*O4*transpose(conj(L));
    
    %%%%%%%%% W3 
    O1 = vl_C_Cbar.*Hlocal{1};
    O2 = transpose(C)*O1*conj(C);
    OL = O2.*eye(D);
    
    O1 = vl_C_Cbar.*eye(D);
    O2 = transpose(C)*O1*conj(C);
    OL = OL + O2.*Hlocal{4};
    
    O1 = vl_C_Cbar.*Hlocal{2};
    O2 = transpose(C)*O1*conj(C);
    OL = OL + O2.*Hlocal{3};
    
    O3 = transpose(OL)*bsxfun(@times, C, R(1,:));
    
    O4 = O3.*(1./C);
    W3 = conj(Rtrunc)*O4*transpose(conj(L));
    
    %%%%%%%%%% W2
    L2 = vl_C_Cbar.*Hlocal{1};  R2 = eye(D).*C_Cbar_vr;
    
    T1 = transpose(L2)*C*R2;
    T2 = T1.*(1./C);
    
    W2 = conj(Rtrunc)*T2*transpose(conj(L));
    
    L2 = vl_C_Cbar.*eye(D);  R2 = Hlocal{4}.*C_Cbar_vr;
    
    T1 = transpose(L2)*C*R2;
    T2 = T1.*(1./C);
    
    W2 = W2+ conj(Rtrunc)*T2*transpose(conj(L));
    
    L2 = vl_C_Cbar.*Hlocal{2};  R2 = Hlocal{3}.*C_Cbar_vr;
    
    T1 = transpose(L2)*C*R2;
    T2 = T1.*(1./C);
    
    W2 = W2+ conj(Rtrunc)*T2*transpose(conj(L));
    
    %%%%%%%%%%%  WF  non-local term    
    Rtrunc_tilde = bsxfun(@times, (1./(1- Eeigvals(2:end))), Rtrunc);
    
    %%OL is recycled from W3...
    F = Rtrunc_tilde*diag(transpose(C)*OL*conj(C));
    T2 = transpose(F)*Ltrunc;
    T3 = bsxfun(@times, transpose(T2), C);
    T3 = bsxfun(@times, T3, R(1,:));
    T4 = T3.*(1./C);
    
    WF = conj(Rtrunc)*T4*transpose(conj(L));
      
    %%%%%%%%%%%%%%%%%
    W = W1+W2+W3+WF;
    W_vec = reshape(transpose(W), D*(D-1), 1);
    
   % B = G_inv_action( W_vec, C, L, R, Rtrunc, bicg_maxiter, bicg_tolerance, Bzero_GLOBAL );
    
    B = G_inv_action_precond( W_vec, C, L, R, Ltrunc, Rtrunc, bicg_maxiter, bicg_tolerance, Bzero_GLOBAL );
    
    Bzero_GLOBAL = B;

    Bmat = transpose(reshape(B, D, D-1));
    deltaC = (transpose(Rtrunc)*Bmat*L).*(1./conj(C));
    
    dC_vec  =1i*reshape(transpose(deltaC), D^2, 1);
    

end

