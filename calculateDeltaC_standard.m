function [ deltaC, B, W, G, BL, BC, BR, F] = calculateDeltaC_standard( C, L, R, Eeigvals, J, h )
    
    [~, D] = size(C);
    siteOverlap = log2(D);
    
    C_Cbar_ur = cell(D-1, 1);
    C_Cbar_ur_vec = cell(D-1, 1);
    vl_C_dCbar = cell(D*(D-1), 1);
    vl_C_dCbar_vec = cell(D*(D-1), 1);
    C_dCbar_vr = cell(D*(D-1), 1);
    C_dCbar_vr_vec = cell(D*(D-1), 1);
    dc = cell((D-1)*D, 1);
    
    %left gauge condition
    counter_tmp=1;
    for kk=1:(D-1)
        for mm=1:D
            dc{counter_tmp} = kron(transpose(R(kk+1,:)), L(mm,:))./conj(C);
            counter_tmp = counter_tmp+1;
        end
    end
    
    %hope we are not dividing by \approx zero
   % fprintf('Min norm squared of element of C is %d\n', min(reshape(transpose(C.*conj(C)), 4, 1)));

    G = zeros(D*(D-1), D*(D-1));
    for kk=1:D*(D-1)
        for mm = 1:D*(D-1)
            G(kk,mm) = L(1,:)*(dc{mm}.*conj(dc{kk}))*transpose(R(1,:));
        end
    end
    
    %BL, BM, BR  
    
    %%vl_C_Cbar %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    vl_C_Cbar = diag(L(1,:))*C;
    vl_C_Cbar = transpose(vl_C_Cbar)*conj(C);
    vl_C_Cbar_vec = reshape(transpose(vl_C_Cbar), 1, D^2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
     
    %%C_Cbar_vr %%%%%%%%%%%%%%%%%%%%%%%%
  
    C_Cbar_vr = C*diag(transpose(R(1,:)));
    C_Cbar_vr = C_Cbar_vr*C';
    C_Cbar_vr_vec = reshape(transpose(C_Cbar_vr), D^2, 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    
    %%C_Cbar_ur %%%%%%%%%%%%%%%%%%%%%%%%
    for mm =1:(D-1)
        C_Cbar_ur{mm} = C*diag(transpose(R(mm+1,:)));    
        C_Cbar_ur{mm} = C_Cbar_ur{mm}*C';
        C_Cbar_ur_vec{mm} = reshape(transpose(C_Cbar_ur{mm}), D^2, 1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    
    %C_H_Cbar%%%%%%%%%%%%%%%%%%    
    [ C_H_Cbar  ] = calculate_H_effective( C, C, siteOverlap, D, J, h );
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    %BL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    vl_C_dCbar_tmp = diag(L(1,:))*C;
    
    for kk=1:D*(D-1)
        vl_C_dCbar{kk} = transpose(vl_C_dCbar_tmp)*conj(dc{kk});
        vl_C_dCbar_vec{kk} = reshape(transpose(vl_C_dCbar{kk}), 1, D^2);
    end
    
    BL = zeros(D*(D-1), 1);  
    for kk=1:D*(D-1)
        BL(kk) = vl_C_dCbar_vec{kk}*C_H_Cbar*C_Cbar_vr_vec;
    end
    
    %BR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    C_dCbar_vr_tmp = C*diag(transpose(R(1,:)));

    for kk=1:D*(D-1)
        C_dCbar_vr{kk} = C_dCbar_vr_tmp*(dc{kk}');
        C_dCbar_vr_vec{kk} = reshape(transpose(C_dCbar_vr{kk}), D^2, 1);
    end
    
    BR = zeros(D*(D-1), 1);  
    for kk=1:D*(D-1)
        BR(kk) = vl_C_Cbar_vec*C_H_Cbar*C_dCbar_vr_vec{kk};
    end
    %BRb = vl_C_Cbar_vec*Cbar_H_C*C_dCbar_vr_b_vec;
    
    %BC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    C_H_dCbar = cell(D*(D-1), 1);

    for mm=1:D*(D-1)
        C_H_dCbar{mm} = calculate_H_effective( C, dc{mm}, siteOverlap, D, J, h );
    end
    
    BC = zeros(D*(D-1), 1); 
    
    for kk=1:D*(D-1)
        BC(kk) = vl_C_Cbar_vec*C_H_dCbar{kk}*C_Cbar_vr_vec;
    end

    %Non-local term%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    F = zeros(D*(D-1), 1); 
    for kk=1:D*(D-1)
        for mm=1:(D-1)
            F(kk) = F(kk)+ (1/(1-Eeigvals(mm+1)))*vl_C_Cbar_vec*C_H_Cbar*C_Cbar_ur_vec{mm}*(L(mm+1,:)*(C.*conj(dc{kk}))*transpose(R(1,:)));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    W = F +BL +BR +BC;
    
    F = transpose(reshape(F, D, D-1));
    BL = transpose(reshape(BL, D, D-1));
    BR = transpose(reshape(BR, D, D-1));
    BC = transpose(reshape(BC, D, D-1));
    
    B = G\W;

    deltaC=zeros(D);
    for kk=1:D*(D-1)
        deltaC = deltaC + B(kk)*dc{kk};
    end

end

