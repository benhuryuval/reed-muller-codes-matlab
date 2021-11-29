function [decodedCodeword,decodedMessage] = rmdec_dumer(LLR,r,m)
% function [message,decodedCodeword,messageLR] = rmdec_dumer(LLR,r,m)
% 
% RMDEC_DUMER decodes an RM(r,m) codeword using Dumer's recursive decoding algorithm.
% 
% Input:
%   LLR -   Log-Likelihood Ratio of channel output.
%   r -     Reed-Muller code order.
%   m -     Reed-Muller code length exponent.
%   mask -  2^m binary vector indicating the locations of frozen bits (0). Default: ones(2^m,1).
% 
% Output:
%   decodedCodeword -   decoded codeword.
%   decodedMessage -    decoded informatio message.
% 
% Source:   [1] Reed-Muller Codes: Theory and Algorithms by Emmanuel Abbe, Amir
%               Shpilka and Min Ye, 2020.
% 
% Written by Yuval Ben-Hur, 26/02/2020
% Updated by Yuval Ben-Hur, 13/05/2020 (r==0, not tested yet)
% 

% fprintf('Decoding RM(%d,%d).\n',r,m);

if nargin<4
    mask = ones(2^m,1);
end

LLR = LLR(:).';

if r>1 && r<m
    Z = de2bi(0 : 1 : 2^m-1,'right-msb').';
    
    % divide LLR according to (u,u+v) partition
    L_zm0 = LLR(Z(m,:)==0); % L^{z_m=0}
    L_zm1 = LLR(Z(m,:)==1); % L^{z_m=1}
    L_slash_zm = log(exp(L_zm0+L_zm1)+1) - log(exp(L_zm0)+exp(L_zm1)); % L^{\z_m}(f) *channel dependent*
    
    % recursively call decoder for subcode v
    [c_hat_slash_zm,av] = rmdec_dumer(L_slash_zm,r-1,m-1);
    
    % Calculate \tilde{L}^{z_m=0} and \tilde{L}^{\z_m}
    L_gal_zm0 = L_zm0 + (-1).^c_hat_slash_zm .* L_zm1; % *channel dependent?*
    L_gal_slash_zm = L_gal_zm0(:,Z(m,:)==0);
    
    % recursively call decoder for subcode u
    [c_hat_zm0,au] = rmdec_dumer(L_gal_slash_zm,r,m-1);
    
    decodedCodeword = [c_hat_zm0, mod(c_hat_zm0+c_hat_slash_zm,2)];
    decodedMessage = rmdemap(decodedCodeword,r,m);
%     decodedMessage_ = [av,au];
    
elseif r==1 % 1st order RM code
    
    [decodedCodeword,decodedMessage] = rmdec_fht(LLR,r,m);
    
elseif r==m % full code

    decodedCodeword = double(LLR<0); % ML decoder
    decodedMessage = rmdemap(decodedCodeword,r,m);
    
elseif r==0 % repetition code
    
    decodedMessage = (1-sign(sum(LLR)))/2;
    decodedCodeword = decodedMessage .* ones(1,2^m);
    
else % decoder should never reach here

    error('Decoder error.');

end
