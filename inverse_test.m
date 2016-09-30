%efficient inverse of product of two matrices inv(G) G = A*B;

clear;

D1=6; D2 = 4; Ddiff = D1-D2;

A = rand(D2, D1); B = rand(D1, D2);

G = A*B;

[U, D, V] = svd(A);


Dtmp = D(:, 1:end-Ddiff);

Dinv = diag(1./diag(Dtmp));

Dinv = [Dinv; zeros(Ddiff, D2)];

P = Dinv*D; %projector

M = P*V'*B;

Mtrunc = M(1:D2, :);

Mtrunc_inv = inv(Mtrunc);

Minv = zeros(D2, D1);
Minv(1:D2, 1:D2)  = Mtrunc_inv;

Ginv = Mtrunc_inv*Dinv(1:D2, :)*U';

