function codeword = rmenc_v2(message,r,m)
% function codeword = rmenc_v2(message,r,m)
% 
% RMENC encodes a binary message to an RM(r,m) codeword by polynomial evaluation.
% 
% Source:   [1] Reed-Muller Codes: Theory and Algorithms by Emmanuel Abbe, Amir
%               Shpilka and Min Ye, 2020.
% 
% Written by Yuval Ben-Hur, 26/02/2020
% Updated by Yuval Ben-Hur, 03/03/2020
% 

Z = de2bi(2^m-1:-1:0,'right-msb').';
[~,AscendingOrder] = sort(sum(Z,1),'ascend');
Z = Z(:,AscendingOrder);

Z_r = Z(:,sum(Z,1)<=r);

codeword = flip(Eval(Z_r(:,logical(message))));