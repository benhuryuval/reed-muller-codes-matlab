function [codeword,G,monomG] = rmenc(message,r,m,G,monomG)
% function [codeword,G,monomG] = rmenc(message,r,m,G,monomG)
% 
% RMENC encodes a binary message to an RM(r,m) codeword using a generator matrix. The generator matrix can
% be supplied as input, or calculated inside the function.
% 
% Source:   [1] Reed-Muller Codes: Theory and Algorithms by Emmanuel Abbe, Amir
%               Shpilka and Min Ye, 2020.
%           [2] Reed-Muller Codes, Sebastian Raaphorst, Carleton University, 2003.
% 
% Written by Yuval Ben-Hur, 26/02/2020
% 

if nargin<5
    monomG = zeros(1,m); % first row is 1, hence corres. monom is [0,0,...,0]
    G = ones(1,2^m); % first row is Eval(1)
    
    Z = de2bi(0 : 1 : 2^m-1,'right-msb').';
    [~,AscendingOrder] = sort(sum(Z,1),'ascend');
    Z = Z(:,AscendingOrder);
    
    for rr = 1 : 1 : r
        for jj = fliplr(find(sum(Z,1)==rr))
            monomG = [monomG; Z(:,jj).'];
            G = [G; Eval(Z(:,jj))];
        end
    end
    G = fliplr(G);
    monomG = fliplr(monomG);
end

codeword = mod(message*G,2);
