function [c, L, LLR_1] = rmdec_rpa(L, m, r, Nmax, theta)
% function [c, L, LLR_1] = rmdec_rpa(L, m, r, Nmax, theta)
%RMDEC_RPA decodes a Reed-Muller codeword using the RPA algorithm.
% 
% Source: "Recursive projection-aggregation decoding of Reed-Muller codes", Min Ye Emmanuel Abbe, 2020
% 
% Inputs:
%   L - input LLR vector (channel output)
%   m - RM code degree
%   r - RM code order
%   Nmax - maximal number of iterations
% 
% Outputs:
%   c - decoded word
%   L - LLR values of decoded word
%   LLR_1 - LLR values after the first iteration of the algorithm
% 
%   Written by: Asaf Goren, Asaf Arad, 2021
%   Updated by: Yuval Ben-Hur, 29/11/2021
% 


    LLR_1 = L; % initialize
    
    % map all the binary field
    Z = de2bi(2^m-1:-1:0,'right-msb').';
    [~,AscendingOrder] = sort(sum(Z,1), 'ascend');
    Z = Z(:,AscendingOrder);
    E_FHT = Z';

    E = de2bi(0:2^m-1,'left-msb');
    if r == 1 % stop condition - FHT decoding
        c = FHTdecoder(L, m, E_FHT);
        return;
    end

    for j = 1:Nmax % loop over iterations
       y_hat_mat = zeros(2^m-1,2^m);
       for i = 1:2^m - 1 % loop over all vectors in E (except 0 vector)
           bi = E(i+1,:);
           B = [E(1,:); bi]; % sub-group for the coset
           % define flags for building the cosets
           T = [b2d(flip(E,2)), b2d(flip((mod(E + bi,2)),2))];
           [T_unique, ia, ic] = unique(sort(T,2), 'rows');
           L_B_T = log(exp(L(1+T_unique(:,1)) + L(1+T_unique(:,2))) + 1) - ...
               log(exp(L(1+T_unique(:,1))) + exp(L(1+T_unique(:,2))));
           
           [y_hat_B, ~] = rmdec_rpa(L_B_T, m-1, r-1, Nmax, theta);
           y_hat_mat(i,:) = y_hat_B(ic);
       end
       
       L_hat = aggregation(L, m, y_hat_mat,E);
       
       if j == 1
           LLR_1 = L_hat; % LLR after the first iteration
       end
       if max(abs(L_hat-L) ./ abs(L)) <= theta
           break;
       end
       L = L_hat;
    end
    
    c = L < 0;
end

function c = FHTdecoder(L, m, E)
    Z = E';
    n = 2^m;
    FHT = fwht(L,n,'hadamard');
    [~,max_ind] = max(abs(FHT));
    u = de2bi(max_ind-1,m,'left-msb'); % FHT Algorithm - left MSB 
    
    Z_1 = Z(:,sum(Z,1)==1);
    c = Eval(Z_1(:,logical(u)));
    if FHT(max_ind) < 0
        c = not(c);
    end
end

function y = b2d(x)
% Convert a binary array to a decimal number
% 
% Similar to bin2dec but works with arrays instead of strings and is found to be 
% rather faster

    z = 2.^(size(x,2)-1:-1:0);
    y = sum(x.*z,2);
end

function L_hat = aggregation(L,m,y_hat_mat,E)

    % initialize cummLLR(z)
    cumLLR = zeros(1,2^m);
    n = 2^m;
    L_hat = zeros(1,2^m);

    y_hat_mat_padded = [zeros(1,2^m) ; y_hat_mat]; % padding with first row of zeros
    for j = 1:n
        z = E(j,:);
        cumLLR(j) = sum((1-2*y_hat_mat_padded(2:end,j))' .* L(bi2de(bitxor(z,E(2:end,:)))+1));
        L_hat(j) = cumLLR(j) / (n-1);
    end
end

