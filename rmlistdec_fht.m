function [decodedCodewordList,decodedMessageList] = rmlistdec_fht(LLR,r,m)
% function [decodedCodeword,decodedMessage] = rmdec_reed(L,r,m,G)
% 
% RMDEC_FHT decodes an RM(1,m) codeword using Fast Hadamard Transform.
% 
% Source: Reed-Muller Codes: Theory and Algorithms by Emmanuel Abbe, Amir
%   Shpilka and Min Ye, 2020.


if r~=1, error('FHT decodes only 1st-order RM codes'); end

LLR = LLR(:).';
hatLu = fwht(LLR,[],'hadamard');
[~,ListOrder] = sort(abs(hatLu),'descend');
uStarList = de2bi(ListOrder-1,m,'right-msb').';

ListSize = 2^(1+m);
decodedCodewordList = zeros(ListSize,2^m);
decodedMessageList = zeros(ListSize,1+m);

poly = eye(m);
for ll = 1 : size(ListOrder,2)
    ustar = uStarList(:,ll);
    if hatLu(ListOrder(ll))>0
        if all(ustar==0)
            decodedCodewordList(ll,:) = zeros(1,2^m);
        else
            decodedCodewordList(ll,:) = Eval(poly(:,ustar==1));
        end
    else
        if all(ustar==0)
            decodedCodewordList(ll,:) = ones(1,2^m);
        else
            decodedCodewordList(ll,:) = Eval([zeros(m,1), poly(:,ustar==1)]);
        end
    end
    decodedCodewordList(size(ListOrder,2) + ll,:) = mod(1+decodedCodewordList(ll,:),2);

    % find information message
    if nargout==2
        k = 0; for rr = 0:1:r, k = k+sum(nchoosek(m,rr)); end
        Z = de2bi(0 : 1 : 2^k-1,'right-msb');
        decodedMessageList(ll,:) = Z(all(rmenc(Z,r,m)==decodedCodewordList(ll,:),2),:);
    end

end

% for ll = 1 : size(ListOrder,2)
%     ustar = uStarList(:,ll);
%     if hatLu(ListOrder(ll))<=0
%         if all(ustar==0)
%             decodedCodewordList(size(ListOrder,2) + ll,:) = zeros(1,2^m);
%         else
%             decodedCodewordList(size(ListOrder,2) + ll,:) = Eval(poly(:,ustar==1));
%         end
%     else
%         if all(ustar==0)
%             decodedCodewordList(size(ListOrder,2) + ll,:) = ones(1,2^m);
%         else
%             decodedCodewordList(size(ListOrder,2) + ll,:) = Eval([zeros(m,1), poly(:,ustar==1)]);
%         end
%     end
% 
%     % find information message
%     if nargout==2
%         k = 0; for rr = 0:1:r, k = k+sum(nchoosek(m,rr)); end
%         Z = de2bi(0 : 1 : 2^k-1,'right-msb');
%         decodedMessageList(size(ListOrder,2) + ll,:) = Z(all(rmenc(Z,r,m)==decodedCodewordList(ll,:),2),:);
%     end
% end

