function [ W_vec ] = calculateDeltaC( C, L, R, Ltrunc, Rtrunc, Eeigvals, Hlocal)
    
    [~, D] = size(C);
    
    %Some useful pre-contractions
 
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
    
    
    %left gauge condition
 %   counter_tmp=1;
 %   for kk=1:(D-1)
 %       for mm=1:D
 %           dc{counter_tmp} = kron(transpose(R(kk+1,:)), L(mm,:))./conj(C);
 %           counter_tmp = counter_tmp+1;
 %       end
 %   end
    
    %hope we are not dividing by \approx zero
   % fprintf('Min norm squared of element of C is %d\n', min(reshape(transpose(C.*conj(C)), 4, 1)));

    %G = zeros(D*(D-1), D*(D-1));
    %for kk=1:D*(D-1)
    %    for mm = 1:D*(D-1)
    %        G(kk,mm) = L(1, :)*(dc{mm}.*conj(dc{kk}))*transpose(R(1,:));
    %    end
    %end
    
 
    %B = G\W_vec;
    
  %  maxiter = 2000; tolerance = 5e-9; vec_zero = cpxrand(D*(D-1), 1);
  %  B = G_inv_action( W_vec, C, L, R, Rtrunc, maxiter, tolerance, vec_zero );

 %   deltaC=zeros(D);
 %   for kk=1:D*(D-1)
 %       deltaC = deltaC + B(kk)*dc{kk};
 %   end
    


end

