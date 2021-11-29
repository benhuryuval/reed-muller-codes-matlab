function [c, v] = rmlistdec_rpa(L, m, r, Nmax, theta, t)
% function [c, v] = rmlistdec_rpa(L, m, r, Nmax, theta, t)
%RMLISTDEC_RPA decodes a Reed-Muller list of codeword using the list decoder variant of the RPA algorithm.
% 
% Source: "Recursive projection-aggregation decoding of Reed-Muller codes", Min Ye Emmanuel Abbe, 2020
% 
% Inputs:
%   L - input LLR vector (channel output)
%   m - RM code degree
%   r - RM code order
%   Nmax - maximal number of iterations
%   t - log2 of list size
% 
% Outputs:
%   c - list of decoded words
%   v - LLR values of decoded words
% 
%   Written by: Asaf Goren, Asaf Arad, 2021
%   Updated by: Yuval Ben-Hur, 29/11/2021
% 

    c_list = zeros(2^t,2^m);
    list_size = 2^t;
    L_tilde = L;
    [~, z] = mink(abs(L),t);

    L_max = 2*max(abs(L));
    for i = 1:list_size % for each u in {L_max,-L_max}^t
        u = ((d2b(i-1,t) * 2) - 1) * L_max;
        L(z) = u;
        [c_hat_u, ~] = rmdec_rpa(L, m, r, Nmax, theta);
        [c_hat_u,~] = rmdec_reed(c_hat_u,r,m);
        c_list(i,:) = c_hat_u; % concat to the decoded list
    end
    c_grades = sum( (-1).^c_list .* L_tilde,2);
    [~, I] = sort(c_grades,'descend');

    c = c_list(I,:);
    v = var(c_grades);
end

function y = d2b(x, m)

% Convert a decimanl number into a binary array
% 
% Similar to dec2bin but yields a numerical array instead of a string and is found to
% be rather faster

    c = floor(log2(x)) + 1;
    if x == 0
        y = 0;
        y = [zeros(1,m-length(y)), y]; % pad with zeros to be length m
        return;
    end
    y(c) = 0; % Initialize output array
    for i = 1:c
        r = floor(x / 2);
        y(c+1-i) = x - 2*r;
        x = r;
    end
    % padding with zeros to length of m
    y = [zeros(1,m-length(y)), y];

end
