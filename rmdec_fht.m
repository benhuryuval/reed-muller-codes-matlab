function [decodedCodeword,decodedMessage] = rmdec_fht(L,r,m)
% function [decodedCodeword,decodedMessage] = rmdec_fht(L,r,m)
% 
% RMDEC_FHT decodes an RM(1,m) codeword using Fast Hadamard Transform.
% 
% Source: Reed-Muller Codes: Theory and Algorithms by Emmanuel Abbe, Amir
%   Shpilka and Min Ye, 2020.


if r~=1, error('FHT decodes only 1st-order RM codes'); end

L = L(:).';
hatLu = fwht(L,[],'hadamard');
[~,ustar_idx] = max(abs(hatLu));
ustar = de2bi(ustar_idx-1,m,'right-msb').';
poly = eye(m);
if hatLu(ustar_idx)>0
    if all(ustar==0)
        decodedCodeword = zeros(1,2^m);
    else
        decodedCodeword = Eval(poly(:,ustar==1));
    end
else
    if all(ustar==0)
        decodedCodeword = ones(1,2^m);
    else
        decodedCodeword = Eval([zeros(m,1), poly(:,ustar==1)]);
    end
end

% find information message
if nargout==2
    k = 0; for rr = 0:1:r, k = k+sum(nchoosek(m,rr)); end
    Z = de2bi(0 : 1 : 2^k-1,'right-msb');
    decodedMessage = Z(all(rmenc(Z,r,m)==decodedCodeword,2),:);
end

% % consider trying: 
% decodedMessage = [double(hatLu(ustar_idx)<=0), ustar(:).'];