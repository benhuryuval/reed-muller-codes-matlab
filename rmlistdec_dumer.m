function [FinalList] = rmlistdec_dumer(LLR,r,m,L)
% function [List] = rmlistdec_dumer(LLR,r,m,L)
% 
% RMLISTDEC_DUMER decodes a list of codewords based on Dumer's recursive-lists algorithm.
% 
% Input:
%   LLR - Input word bits log-likelihood ratio (also denoted y_in)
%   r - Reed-Muller code order
%   m - Reed-Muller code length exponent
%   L - integer <2^k defining the maximal list size, where
%           k = sum_{j=0,...,r} {nchoosek(m,j)}
% 
% Output:
%   List - List of decoded codewords ordered by descending reliability.
%          Type is struct with fields:
%               a - information
%               rho - cost (or reliability)
%               y - codeword
% 
% Sources: 
%       [1] Reed-Muller Codes: Theory and Algorithms by Emmanuel Abbe, Amir
%           Shpilka and Min Ye, 2020.
%       [2] Soft decision decoding of Reed-Muller codes: recursive lists by 
%           Ilya Dumer and Kirill Shabunov, 2017.
%       [3] ECCLab\dtrm_glp_bg - Simulation of recursive list decoding for
%           RM codes, 2020.
% 
% Written by Yuval Ben-Hur, 25/09/2020
%     

    % define persistent (static) variable path counter (s)
    global node_counter
    if isempty(node_counter), node_counter=0; end
    
    % initialize lists
    Lists = struct( 'ListInfo',                         ... % list of information bits properties
                                struct( 'a', 0,         ... % information bits
                                        'lind2xl',0,    ... % list index corresponding to information bits
                                        'parent',-1,    ... % list item parents
                                        'rho', 0 )      ... % confidence of bit sequence
                                        ,               ...
                    'ListItems',    repmat(             ... % lists of LLRs for codes 1,...,m
                                        struct('yList',[], 'ListLen',0) ...
                                    ,[1,m])             ...
                    );
    
    % calculate y_v = y' + y''
    Lists.ListItems(m-1).yList( Lists.ListItems(m-1).ListLen+1 ) = left(LLR) .* right(LLR);
    Lists.a = 0;
    Lists.ListItems(m-1).ListLen = Lists.ListItems(m-1).ListLen + 1;
    
    % decode v
    Lists = rmlistdec_dumer_inner(Lists,r-1,m-1,L);
    
    % calculate y_u = (y'+v) + y''
    for idx = 1 : Lists.ListItems(m-1).ListLen
        y1 = left(LLR) .* Lists.ListItems(m-1).yList(idx);
        Lists.ListItems(m-1).yList(idx) = (y1 + right(LLR)) ./ (1 + y1.*right(LLR)) ;
    end
    
    % decode u
    Lists = rmlistdec_dumer_inner(Lists,r,m-1,L);
    
    % return sorted list
    [~,order] = sort([Lists.ListInfo.rho],'ascend');
    FinalList = struct( 'a', [Lists.ListInfo.a],...
                        'rho', [Lists.ListInfo.rho],...
                        'ListItems', Lists.ListItems(m).yList,...
                        'Order', order );
    
    clear rmlistdec_dumer % discard persistant variables

% - * - * - Auxilliary functions - * - * -
function [Lists] = rmlistdec_dumer_inner(Lists,r,m,L)
    
    if r==m
        if node_counter==0
            rmm_skip(Lists,m);
        else
            if m==1
                rm11_branch(Lists,node_counter);
            else
                rmm_branch(Lists,m,node_counter);
            end
        end
        node_counter = node_counter+1;
        return;
    end
    
    if r==1
        if node_counter==0
            rm1_skip(Lists,m);
        else
            rm1_branch(Lists,m,node_counter);
        end
        node_counter = node_counter+1;
    else
        
        % calculate y_v = y' + y''
        for idx = 1 : size(List,1)
            List(idx).v = left(List(idx).y) + right(List(idx).y);
        end
        
        % Decode v
        List = rmlistdec_dumer_inner(List,r-1,m-1,L);
        
        % Calculate y_u = y_1 xor v + y_2.
        for idx = 1 : size(List,1)
            List(idx).u = (List(idx).v + left(List(idx))) + right(List(idx));
        end
        
        % Decode u
        List = rmlistdec_dumer_inner(List,r,m-1,L);

        % y_dec = (u xor v | u)
        % TODO
        
    end
    
function [Lists] = rm11_branch(Lists,L)
    n = 2; % code length
    binmask = de2bi(0:1:2^n-1); % binary mask for bit flipping
    % iterate list and decode information bits
    for idx = 1 : length(Lists(1).ListItems)
        ytmp = abs(Lists(1).ListItems(idx).y);
        xtmp = double(Lists(1).ListItems(idx).y <= 0); % 1 if y<=0, 0 if y>0
        % add all possible flips of xtmp to list
        for jdx = 1 : size(binmask,1)
            Lists.a(length(Lists.a)+1) = [Lists.a(length(Lists.a)+1), xtmp+binmask(jdx,:)];
            Lists.lind2xl(length(Lists.rho)+1) = idx;
            Lists.parent(length(Lists.rho)+1) = idx;
            Lists.rho(length(Lists.rho)+1) = Lists.rho(length(Lists.rho)+1) + sum(lnp(ytmp,binmask(jdx,:)));
        end
    end
    % sort and cut list (originally, partition and cut)
    [~,order] = sort([Lists.ListInfo.rho],'ascend');
    Lists.ListInfo = Lists.ListInfo(order(1:L));
    
function z = lnp(y,b)
% function z = lnp(y,b)
% LNP Calculates ln(Prob{y=b|y})
    z = log(1 + sign(0.5-b).*y);

% functions that slice a vector y to left and right parts
function yy = left(y)
    yy = y(1 : length(y)/2);

function yy = right(y)
    yy = y(length(y)/2+1 : length(y));


function b = partition(a,k)
    st = 0;
    for ii = 1 : length(a)
        if a(ii) >= a(end)
            tmp = a(ii); a(ii) = a(st); a(st) = tmp; % swap a(ii) and a(st)
            st = st + 1;
        end
    end
    if k==st, return; end
    if st > k,  b = partition(a(1:st),k);
    else,       b = partition(a(st+1:end),k-st);
    end
    
% TODO:
function [Lists] = rm1_branch(Lists,m)
    for idx = 1 : length(Lists(m).ListItems)
        v = left(Lists(m).ListItems(idx).y) + right(Lists(m).ListItems(idx).y);
        if sum(log(1 + v)) > sum(log(1 - v)) % s0>s1
            Lists(m).ListItems(idx).a = [Lists(m).ListItems(idx).y ]
        end
    end
    
function [Lists] = rmm_skip(Lists,m)
    for idx = 1 : length(Lists(m).ListItems)
        Lists(m).ListItems(idx).y = left(Lists(m).ListItems(idx).y) + right(Lists(m).ListItems(idx).y);
    end


