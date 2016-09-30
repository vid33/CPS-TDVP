function W_out=W_action_RQ2(x_in, mass, potential, interaction, D, REAL_TIME)

        Rvec = x_in(1:D^2); Qvec = x_in(D^2+1:end);    

        R = transpose(reshape(Rvec, D, D));

        Q = transpose(reshape(Qvec, D, D));
        
        %impose matl = eye(D) exactly
        Ktmp = (1i/2)*(Q - Q');
        Ktmp = (1/2)*(Ktmp+Ktmp');
        Q = -(1/2)*(R')*R - 1i*Ktmp;
               
        Rkin = Q*R - R*Q;       
        
        %function [ matr ] = lrEigenvector_fast_load( R, Q, D)
       
        [ matr ] = eigenvector_r_fast_load( R, Q, D);
        
             
        %[F, ~] = calculateF( R, Rkin, Q, matr, potential, interaction, F_zero, D );
        
       % function [ F ] = calculateF_load( R, Rkin, Q, matr, potential, interaction, D )

        [ F ] = calculateF_load( R, Rkin, Q, matr, potential, interaction, D );

        matF = transpose(reshape(F, D, D));
        
    
        if (REAL_TIME == false) %imaginary time evolution
            Rstar = (matr*(R'))/matr;
            
            W_out_R = -1*( transpose(matF)*R-R*transpose(matF)...
                +(1/(2*mass))*((Rkin*Q - Q*Rkin) + (Rkin*R - R*Rkin)*Rstar...
                + (R*R' - R'*R)*Rkin)...
                +potential*R + interaction*(R*R*Rstar + R'*R*R) );
        
            W_out_Q = -Rstar_l*W_out_R;
        
        elseif (REAL_TIME == true) %real time evolution
            
            Rstar = (matr*(R'))/matr;
            
            W_out_R = 1i*( transpose(matF)*R-R*transpose(matF)...
                +(1/(2*mass))*((Rkin*Q - Q*Rkin) + (Rkin*R - R*Rkin)*Rstar...
                + (R*R' - R'*R)*Rkin)...
                +potential*R + interaction*(R*R*Rstar + R'*R*R) );
        
            W_out_Q = -R'*W_out_R;
            
        end
        
        W_out_R = reshape(transpose(W_out_R), D^2, 1);
        W_out_Q = reshape(transpose(W_out_Q), D^2, 1);
        
        W_out = [W_out_R ; W_out_Q];
       
end

