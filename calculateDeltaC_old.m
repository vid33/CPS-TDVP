function [ deltaC, B, W, G] = calculateDeltaC_old( C, H, vl, vr, ul, ur, lm )

    %initialization
    vl_C_Cbar = zeros(2,2);
    C_Cbar_vr = zeros(2,2);
    C_Cbar_ur = zeros(2,2);
    vl_C_dCbar = zeros(2,2);
    C_dCbar_vr = zeros(2,2);
    G = zeros(2,2);

    %LGC
    dca = kron(ur, vl)./conj(C);
    dcb = kron(ur, ul)./conj(C);
   
    %hope we are not dividing by \approx zero
   % fprintf('Min norm squared of element of C is %d\n', min(reshape(transpose(C.*conj(C)), 4, 1)));

    G(1,1) = vl*(dca.*conj(dca))*vr;
    G(1,2) = vl*(dcb.*conj(dca))*vr;
    G(2,1) = vl*(dca.*conj(dcb))*vr; 
    G(2,2) = vl*(dcb.*conj(dcb))*vr;
    
    %BL, BM, BR  
    
    %%vl_C_Cbar %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    for kk=1:2
        vl_C_Cbar(kk, :) = vl(kk)*C(kk, :); 
    end
    vl_C_Cbar = transpose(vl_C_Cbar)*conj(C);
    vl_C_Cbar_vec = reshape(transpose(vl_C_Cbar), 1, 4);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
     
    %%C_Cbar_vr %%%%%%%%%%%%%%%%%%%%%%%%
    for kk=1:2
        C_Cbar_vr(:, kk) = C(:, kk)*vr(kk);
    end
    C_Cbar_vr = C_Cbar_vr*C';
    C_Cbar_vr_vec = reshape(transpose(C_Cbar_vr), 4, 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    %%C_Cbar_ur %%%%%%%%%%%%%%%%%%%%%%%%
    for kk=1:2
        C_Cbar_ur(:, kk) = C(:, kk)*ur(kk);
    end
    C_Cbar_ur = C_Cbar_ur*C';
    C_Cbar_ur_vec = reshape(transpose(C_Cbar_ur), 4, 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    
    %Cbar_H_C%%%%%%%%%%%%%%%%%%    
    Cbar_H_C = H;
    Cbar_covec = reshape(transpose(conj(C)), 1, 4);
    C_vec = reshape(transpose(C), 4, 1);
    for kk=1:4
       Cbar_H_C(kk, :) = Cbar_covec(kk)*Cbar_H_C(kk,:); 
    end
    for kk=1:4
        Cbar_H_C(:, kk) = C_vec(kk)*Cbar_H_C(:, kk);
    end
    %rearrange indices
    Cbar_H_C = reshape(Cbar_H_C, [2, 2, 2, 2]);
    tmp = Contract({Cbar_H_C},{ [-1,-2, -3, -4] } );
    tmp = Contract({tmp},{ [-3,-1, -4, -2] } );
    Cbar_H_C = reshape(tmp, 4, 4); 
    clearvars tmp;
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    %BL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for kk=1:2
        vl_C_dCbar(kk, :) = vl(kk)*C(kk, :); 
    end
    vl_C_dCbar_a = transpose(vl_C_dCbar)*conj(dca);
    vl_C_dCbar_b = transpose(vl_C_dCbar)*conj(dcb);
    vl_C_dCbar_a_vec = reshape(transpose(vl_C_dCbar_a), 1, 4);
    vl_C_dCbar_b_vec = reshape(transpose(vl_C_dCbar_b), 1, 4);
    
    BLa = vl_C_dCbar_a_vec*Cbar_H_C*C_Cbar_vr_vec;
    BLb = vl_C_dCbar_b_vec*Cbar_H_C*C_Cbar_vr_vec;
    
    %BR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for kk=1:2
        C_dCbar_vr(:, kk) = C(:, kk)*vr(kk);
    end

    C_dCbar_vr_a = C_dCbar_vr*dca';
    C_dCbar_vr_b = C_dCbar_vr*dcb';
    C_dCbar_vr_a_vec = reshape(transpose(C_dCbar_vr_a), 4, 1);
    C_dCbar_vr_b_vec = reshape(transpose(C_dCbar_vr_b), 4, 1);
    
    BRa = vl_C_Cbar_vec*Cbar_H_C*C_dCbar_vr_a_vec;
    BRb = vl_C_Cbar_vec*Cbar_H_C*C_dCbar_vr_b_vec;
    
    %BC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    dCbar_H_C_a = H;
    dCbar_H_C_b = H;
    dCbar_covec_a = reshape(transpose(conj(dca)), 1, 4);
    dCbar_covec_b = reshape(transpose(conj(dcb)), 1, 4);
    C_vec = reshape(transpose(C), 4, 1);
    for kk=1:4
       dCbar_H_C_a(kk, :) = dCbar_covec_a(kk)*dCbar_H_C_a(kk,:); 
       dCbar_H_C_b(kk, :) = dCbar_covec_b(kk)*dCbar_H_C_b(kk,:); 
    end
    for kk=1:4
        dCbar_H_C_a(:, kk) = C_vec(kk)*dCbar_H_C_a(:, kk);
        dCbar_H_C_b(:, kk) = C_vec(kk)*dCbar_H_C_b(:, kk);
    end
    %rearrange indices
    dCbar_H_C_a = reshape(dCbar_H_C_a, [2, 2, 2, 2]);
    tmp = Contract({dCbar_H_C_a},{ [-1,-2, -3, -4] } );
    tmp = Contract({tmp},{ [-3,-1, -4, -2] } );
    dCbar_H_C_a = reshape(tmp, 4, 4);
    
    dCbar_H_C_b = reshape(dCbar_H_C_b, [2, 2, 2, 2]);
    tmp = Contract({dCbar_H_C_b},{ [-1,-2, -3, -4] } );
    tmp = Contract({tmp},{ [-3,-1, -4, -2] } );
    dCbar_H_C_b = reshape(tmp, 4, 4 );   
    clearvars tmp;
    
    BCa = vl_C_Cbar_vec*dCbar_H_C_a*C_Cbar_vr_vec;
    BCb = vl_C_Cbar_vec*dCbar_H_C_b*C_Cbar_vr_vec;
    

    %Non-local term%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    Fa = (1/(1-lm))*vl_C_Cbar_vec*Cbar_H_C*C_Cbar_ur_vec*(ul*(C.*conj(dca))*vr);
    Fb = (1/(1-lm))*vl_C_Cbar_vec*Cbar_H_C*C_Cbar_ur_vec*(ul*(C.*conj(dcb))*vr); %this is zero dcb also obeys r gauge condition
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    W = zeros(2,1);
    
    W(1) = Fa + BLa + BCa + BRa;
    W(2) = Fb + BLb + BCb + BRb;

    B = G\W;

    deltaC = B(1)*dca  + B(2)*dcb;
    

end

