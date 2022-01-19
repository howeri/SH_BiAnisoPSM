function [X_hat, b_hat, res] = OMP_modified(A, b, k, l)  % A*X=b, |X_group|_0=k, l=number of columns in 1 group
assert(l>1);
Lambda = []; %Lambda matraix
r = b;       %residual
d = [];      %index of group in A

A_group = reshape(A, [size(A, 1), l, size(A, 2)/l]); % Devide A into groups
A_group_orth = zeros(size(A_group));
for i=1:size(A, 2)/l
    orthonormalBasis_i = orth(A_group(:,:,i));
    A_group_orth(:, 1:size(orthonormalBasis_i,2), i) = orthonormalBasis_i; % When A_group is not indepednet, extra columns are zeros.
end

for i = 1:k
    inner = sum(A_group_orth.*conj(r), 1);
    align = squeeze(vecnorm(inner));    
    [~, ind] = max(align);
    d = [d ind];
    Lambda = [Lambda A_group(:,:,ind)];
    Lambda_orth = orth(Lambda);
    b_hat_temp = Lambda_orth*(Lambda_orth'*b);  %takes care of complex numbers
    r = b - b_hat_temp;
end

X_hat_nnzOnly = reshape(pinv(Lambda)*b, [l, k]);
X_hat_full = zeros(l,size(A, 2)/l);
X_hat_full(:,d) = X_hat_nnzOnly;
X_hat = reshape(X_hat_full, [size(A, 2), 1]);
b_hat = A*X_hat;
res = norm(b - b_hat);

end

