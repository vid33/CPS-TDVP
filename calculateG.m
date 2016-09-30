function [ G ] = calculateG( C, R, L )

    [~, D] = size(C);

    dc = cell((D-1)*D, 1);

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


end

