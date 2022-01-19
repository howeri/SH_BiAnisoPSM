function [a, b] = NorIterAlg(e, Y_obs, Y_exc, r, K)
% Solves bilinear equations a'(Qm)b=rm
SH_num = size(Y_obs, 1);
len_sol = K*SH_num;
a = ones(len_sol,1)*sqrt(1/len_sol);
iterNum = 30;
disp(['Intinal cost = ', num2str(sum(r.^2)/2/length(r))]);
for n =1:iterNum
    Q_linear = zeros(length(r), len_sol);
    for k =1:K
        a_k = a((k-1)*SH_num+1:k*SH_num);
        e_k = reshape(e(k,:), [1,1,size(e, 2)]);
        Q_linear(:,(k-1)*SH_num+1:k*SH_num) = squeeze(sum(a_k.*(Y_exc.*Y_obs.*e_k), 1))';
    end
    b = pinv(Q_linear)*r;
    
    Q_linear = zeros(length(r), len_sol);
    for k =1:K
        b_k = b((k-1)*SH_num+1:k*SH_num);
        e_k = reshape(e(k,:), [1,1,size(e, 2)]);
        Q_linear(:,(k-1)*SH_num+1:k*SH_num) = squeeze(sum(b_k'.*(Y_exc.*Y_obs.*e_k), 2))';
    end
    a = pinv(Q_linear)*r;
    
    ind = find(a, 1);
    xi = sign(a(ind));
    a = xi*a/norm(a);
    b = xi*b*norm(a);
    
    r_modeled_temp = zeros(length(r), len_sol);
    for k =1:K
        a_k = a((k-1)*SH_num+1:k*SH_num);
        e_k = reshape(e(k,:), [1,1,size(e, 2)]);
        r_modeled_temp(:,(k-1)*SH_num+1:k*SH_num) = squeeze(sum(a_k.*(Y_exc.*Y_obs.*e_k), 1))';
    end
    r_modeled = r_modeled_temp*b;
    cost = sum((r_modeled - r).^2)/2/length(r);
    disp(['Iter', num2str(n), ' cost = ', num2str(cost)])
end
end