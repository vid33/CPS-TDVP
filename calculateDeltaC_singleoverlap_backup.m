function [ deltaC, B, W, G] = calculateDeltaC( C, L, R, Eeigvals, J, h )

%old function [ deltaC, B, W, G ] = calculateDeltaC( C, H, vl, vr, ul, ur, lm )
    %initialization
    
    [~, D] = size(C);
    
    C_Cbar_ur = cell(D-1, 1);
    C_Cbar_ur_vec = cell(D-1, 1);
    vl_C_dCbar = cell(D*(D-1), 1);
    vl_C_dCbar_vec = cell(D*(D-1), 1);
    C_dCbar_vr = cell(D*(D-1), 1);
    C_dCbar_vr_vec = cell(D*(D-1), 1);

    G = zeros(D*(D-1), D*(D-1));

    dc = cell((D-1)*D, 1);
    
    %left gauge condition
    for kk=1:(D-1)
        for mm=1:D
            dc{mm*kk} = kron(R{kk+1}, L{mm})./conj(C);
        end
    end
    
    
    %hope we are not dividing by \approx zero
   % fprintf('Min norm squared of element of C is %d\n', min(reshape(transpose(C.*conj(C)), 4, 1)));

    for kk=1:D*(D-1)
        for mm = 1:D*(D-1)
            G(kk,mm) = L{1}*(dc{mm}.*conj(dc{kk}))*R{1};
        end
    end
    
    %BL, BM, BR  
    
    %%vl_C_Cbar %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    vl_C_Cbar = zeros(D);
    for kk=1:D
        vl_C_Cbar(kk, :) = L{1}(kk)*C(kk, :); 
    end
    vl_C_Cbar = transpose(vl_C_Cbar)*conj(C);
    vl_C_Cbar_vec = reshape(transpose(vl_C_Cbar), 1, D^2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
     
    %%C_Cbar_vr %%%%%%%%%%%%%%%%%%%%%%%%
    C_Cbar_vr = zeros(D);
    for kk=1:D
        C_Cbar_vr(:, kk) = C(:, kk)*R{1}(kk);
    end
    C_Cbar_vr = C_Cbar_vr*C';
    C_Cbar_vr_vec = reshape(transpose(C_Cbar_vr), D^2, 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    %%C_Cbar_ur %%%%%%%%%%%%%%%%%%%%%%%%
    for mm =1:(D-1)
        C_Cbar_ur{mm} = zeros(D);
        for kk=1:D
            C_Cbar_ur{mm}(:, kk) = C(:, kk)*R{mm+1}(kk);
        end
    
        C_Cbar_ur{mm} = C_Cbar_ur{mm}*C';
        C_Cbar_ur_vec{mm} = reshape(transpose(C_Cbar_ur{mm}), D^2, 1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    
    %C_H_Cbar%%%%%%%%%%%%%%%%%%    

    sx = [0,1;1,0]; 
    sx_vec_l = reshape(transpose(sx), 1, D^2); sx_vec_r = reshape(transpose(sx), D^2, 1);
    sz = [1,0;0,-1]; 
    sz_vec_l = reshape(transpose(sz), 1, D^2); sz_vec_r = reshape(transpose(sz), D^2, 1);
    eye_vec_l = reshape(eye(D), 1, D^2); eye_vec_r = reshape(eye(D), D^2, 1);
    
    
    C_kron_Cbar = kron(C, conj(C));
    C_H_Cbar_1 = zeros(D^2);
    C_H_Cbar_2l = zeros(D^2);
    C_H_Cbar_2r = zeros(D^2);
    
    for kk=1:D^2
        C_H_Cbar_1(kk, :) = sx_vec_l(kk)*C_kron_Cbar(kk, :); 
    end
    
    for kk=1:D^2
        C_H_Cbar_1(:, kk) = C_H_Cbar_1(:, kk)*sx_vec_r(kk);
    end

    
    for kk=1:D^2
        C_H_Cbar_2l(kk, :) = sz_vec_l(kk)*C_kron_Cbar(kk, :); 
    end
    for kk=1:D^2
        C_H_Cbar_2l(:, kk) = C_H_Cbar_2l(:, kk)*eye_vec_r(kk);
    end
    
    for kk=1:D^2
        C_H_Cbar_2r(:, kk) = C_kron_Cbar(:, kk)*sz_vec_r(kk);
    end
    for kk=1:D^2
        C_H_Cbar_2r(kk, :) = eye_vec_l(kk)*C_H_Cbar_2r(kk, :); 
    end
    
    C_H_Cbar = -J*C_H_Cbar_1 + (h/2)*(C_H_Cbar_2l + C_H_Cbar_2r);
    
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    %BL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    vl_C_dCbar_tmp = zeros(D);
    for kk=1:D
        vl_C_dCbar_tmp(kk, :) = L{1}(kk)*C(kk, :); 
    end
    for kk=1:D*(D-1)
        vl_C_dCbar{kk} = transpose(vl_C_dCbar_tmp)*conj(dc{kk});
        vl_C_dCbar_vec{kk} = reshape(transpose(vl_C_dCbar{kk}), 1, D^2);
    end
    
    BL = zeros(D*(D-1), 1);  
    for kk=1:D*(D-1)
        BL(kk) = vl_C_dCbar_vec{kk}*C_H_Cbar*C_Cbar_vr_vec;
    end
    %BLa = vl_C_dCbar_a_vec*Cbar_H_C*C_Cbar_vr_vec;
    
    %BR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C_dCbar_vr_tmp = zeros(D);
    for kk=1:D
        C_dCbar_vr_tmp(:, kk) = C(:, kk)*R{1}(kk);
    end

    for kk=1:D*(D-1)
        C_dCbar_vr{kk} = C_dCbar_vr_tmp*((dc{kk})');
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
    
        C_kron_dCbar = kron(C, conj(dc{mm}));
        C_H_dCbar_1 = zeros(D^2, D^2);
        C_H_dCbar_2l = zeros(D^2, D^2);
        C_H_dCbar_2r = zeros(D^2, D^2);

        for kk=1:D^2
            C_H_dCbar_1(kk, :) = sx_vec_l(kk)*C_kron_dCbar(kk, :); 
        end

        for kk=1:D^2
            C_H_dCbar_1(:, kk) = C_H_dCbar_1(:, kk)*sx_vec_r(kk);
        end


        for kk=1:D^2
            C_H_dCbar_2l(kk, :) = sz_vec_l(kk)*C_kron_dCbar(kk, :); 
        end
        for kk=1:D^2
            C_H_dCbar_2l(:, kk) = C_H_dCbar_2l(:, kk)*eye_vec_r(kk);
        end

        for kk=1:D^2
            C_H_dCbar_2r(:, kk) = C_kron_dCbar(:, kk)*sz_vec_r(kk);
        end
        for kk=1:D^2
            C_H_dCbar_2r(kk, :) = eye_vec_l(kk)*C_H_dCbar_2r(kk, :); 
        end

        C_H_dCbar{mm} = -J*C_H_dCbar_1 + (h/2)*(C_H_dCbar_2l + C_H_dCbar_2r);
    
    
    end
    
    BC = zeros(D*(D-1), 1); 
    
    for kk=1:D(D-1)
        BC(kk) = vl_C_Cbar_vec*C_H_dCbar{kk}*C_Cbar_vr_vec;
    end
    

    %Non-local term%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    F = zeros(D*(D-1), 1);
    
    for kk=1:D*(D-1)
        for mm=1:(D-1)
            F(kk) = (1/(1-Eeigvals(mm+1)))*vl_C_Cbar_vec*C_H_Cbar*C_Cbar_ur_vec{mm}*(L{mm+1}*(C.*conj(dc{kk}))*R{1});
        end
    end
        %   Fa = (1/(1-lm))*vl_C_Cbar_vec*Cbar_H_C*C_Cbar_ur_vec*(ul*(C.*conj(dca))*vr);
 %   Fb = (1/(1-lm))*vl_C_Cbar_vec*Cbar_H_C*C_Cbar_ur_vec*(ul*(C.*conj(dcb))*vr); %this is zero dcb also obeys r gauge condition
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    W = F +BL +BR +BC;
    
 %   W = zeros(2,1);
    
 %   W(1) = Fa + BLa + BCa + BRa;
 %   W(2) = Fb + BLb + BCb + BRb;

    B = G\W;

    deltaC=zeros(D);
    for kk=1:D*(D-1)
        deltaC = deltaC + B(kk)*dc{kk};
    end

end

