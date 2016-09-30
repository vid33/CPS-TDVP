%for overlap 2 decompose into 2 mps tensors A1 and A2

clear;

A1 = zeros(2, 2, 2, 2, 2);
A2 = zeros(2, 2, 2, 2, 2);

for k=1:2
    A1(k, k, k, k, k) = 1;
    A2(k, k, k, k, k) = 1;
end


A1(2, 1, 2, 2, 1) = 1;
A1(1, 2, 1, 1, 2) = 1;

A2(2, 1, 1, 2, 1) = 1;
A2(1, 2, 2, 1, 2) = 1;


Ctensor = Contract({A1, A2}, {[-1, -2, -3, 1, 2], [1, 2, -4, -5, -6]});

C = reshape(Ctensor,  4, 4, 4);


EE_CPS = Contract({C, conj(C)}, {[-1, 1, -3], [-2, 1, -4]});
EE_CPS = reshape(EE_CPS, 4^2, 4^2);

eig(EE_CPS)

%"small" transfer matrices

EE1 = Contract({A1, conj(A1)}, {[-1, -2,  1, -5, -6], [-3, -4, 1, -7, -8]});
EE1 = reshape(EE1, 4^2, 4^2);


break;

for kk=1:D
    for mm=1:D
        for nn=1:D
            if kk==mm && kk==nn
                bigDelta(kk,mm,nn)=1;
            end
        end
    end
end

bigC = Contract({C, bigDelta}, {[-1, 1], [1, -2, -3]});

bigCbar = Contract({conj(C), bigDelta}, {[-1, 1], [1 -2, -3]});

bigEE = Contract({bigC, bigCbar}, {[-1, 1, -3], [-2, 1, -4]});


bigEEmat = transpose(reshape(bigEE, D^2, D^2));



bigA = Contract({A, bigDelta}, {[-1, 1], [1, -2, -3]});

bigB = Contract({B, bigDelta}, {[-1, 1], [1, -2, -3]});


bigEEmix = Contract({bigA, bigB}, {[-1, 1, -3], [-2, 1, -4]});


bigEEmix_mat = transpose(reshape(bigEEmix, D^2, D^2));


EEmix = A.*B;


