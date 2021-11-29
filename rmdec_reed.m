function [decodedCodeword,decodedMessage] = rmdec_reed(y,r,m)
% function [decodedCodeword,decodedMessage] = rmdec_reed(codeword,r,m)
% 
% RMDEC_REED decodes an RM(r,m) codeword using Reed's majority-logic algorithm.
% 
% Source:   [1] Reed-Muller Codes: Theory and Algorithms by Emmanuel Abbe, Amir
%               Shpilka and Min Ye, 2020.
%           [2] Reed-Muller Codes, Sebastian Raaphorst, Carleton University, 2003.
% 
% Written by Yuval Ben-Hur, 03/03/2020
% 

if ~isequal(size(y),[2^m,1]) && ~isequal(size(y),[1,2^m])
    error('Invalid input size.');
end
if ~all(y==1 | y==0)
    error('Only binary RM supported.'); 
end

y = y(:).';

Z = de2bi(2^m-1:-1:0,'right-msb').';
[~,AscendingOrder] = sort(sum(Z,1),'ascend');
Z_ = Z(:,AscendingOrder);
Z_r = Z_(:,sum(Z_,1)<=r);

kRM = size(Z_r,2);
decodedMessage = zeros(1,kRM); % initialize decoded message
for ii = kRM : -1 : 1
    % define monom
    A = find(Z_r(:,ii)); % monom variables
    t = numel(A); % monom rank
    % define cosets
    logic = [];
    if t==0 % A is the empty set {}
        z = 1 : 2^m;
        logic = y(z);
    elseif t==m % A is the full set {1,...,m}
        logic = mod(sum(y),2);
    else
        cosets = de2bi(0:1:2^(m-t)-1,m-t,'right-msb').';
        % sum over each coset
        for b = cosets
            z = all(Z(setdiff(1:m,A),:) == b,1);
            logic = [logic mod(sum(y(z)),2)];
        end
    end
    decodedMessage(ii) = mode(logic); % decoded bit is majority of logic
    % if next rank is different, reduce monomials of previous rank
    if ii==1 || numel(find(Z_r(:,ii-1)))<t
        y = mod(y - rmenc_v2(decodedMessage.*(sum(Z_r,1)==t),r,m),2);
    end
end

decodedCodeword = rmenc_v2(decodedMessage,r,m);
