function [code,information,scores] = rmlistdec_dumerexe(in,r,m,mode,L)
% function [code,information,scores] = rmlistdec_dumerexe(L,r,m)
% 
% RMLISTDEC_DUMEREXE encodes and/or decodes based on Dumer's recursive list-decoding algorithm.
% It uses a C implementation by Kirill Shabunov downloaded from GitHub.
% Source address: https://github.com/kshabunov/ecclab
% 
% Instructions for using EXE:
%   Encode - File.exe -e r m list_size mask info_bits(unspaced)
%   Decode - File.exe -d r m list_size mask code(spaced floating point values)
% 
% Input:
%   in -    information bits for encoding (in case of mode='enc') or 
%           channel bits for decoding (in case of mode='dec')
%   r -     Reed-Muller code order
%   m -     Reed-Muller code length exponent
%   mode -  string indicating the requested operation. set 'enc' for
%           encoding and 'dec' for decoding.
%   L -     List size (unecessary for mode="enc")
% 
% Output:
%   out -   information bits for encoding (in case of mode='enc') or 
%           list of decoded codewords in descending score (in case of mode='dec').
% 
% Sources: 
%       [1] Reed-Muller Codes: Theory and Algorithms by Emmanuel Abbe, Amir
%           Shpilka and Min Ye, 2020.
%       [2] Soft decision decoding of Reed-Muller codes: recursive lists by 
%           Ilya Dumer and Kirill Shabunov, 2017
% 
% Written by Yuval Ben-Hur, 02/12/2020
% with the help of Neriya Golan.
% 
% Description of tests: 
%       The function was tested on a number of codewords without channel
%       errors. It was also tested when a random number of channel errors, 
%       smaller than dmin/2, was applied.
%       Tests run:
%       (2,3), (2,4), (3,4) - all codewords in code without channel errors and with channel errors.
%       (2,5), (3,5), (4,5), (2,6), (3,6), (4,6), (5,6), (2,7), (3,7), (4,7), (5,7) - 20 codewords from code with channel errors. decoder failed for (6,7).
%       (2,8), (3,8), (4,8), (5,8) - 20 codewords from code with channel errors. decoder failed for (6,8) (7,8).
%       (2,9), (3,9), (4,9), (5,9), (6,9) - 20 codewords from code with channel errors. decoder failed for (7,9) (8,9).
%       (2,10), (3,10), (4,10), (5,10) - 20 codewords from code with channel errors. decoder failed for (6,10), (7,10) (8,10), (9,10).
% 
    
    % check input r,m
    if r>m, error('Doesnt work for specified r, m.'); end
    switch m
        case {2,3,4,5,6}
            
        case {7,8,9}
            if r>=6, error('Doesnt work for specified r, m.'); end            
        otherwise
            error('Doesnt work for specified r, m.');
    end
    
    if nargin<5, L=1; end

    exepath = 'C:\Users\Yuval\Google Drive\PhD\RetentionPaper\attraction-repo\Matlab\functions\reedmuller\RMEncodeDecode.exe';
    if ~isfile(exepath)
        error('Can''t find RMEncodeDecode.exe.');
    end
    exepathStr = ['"' exepath '"'];
    
    maskStr = sprintf('%d',ones(1,numnodes(r,m)));
    rStr = sprintf('%d',r);
    mStr = sprintf('%d',m);
    LStr = sprintf('%d',L);
    switch mode
        case "enc"
            % prepare input bits
            inStr = sprintf('%d',in);
            % encode
            [status, cmdout] = system(strcat(exepathStr," -e ",rStr," ",mStr," ",LStr," ",maskStr," ",inStr));
            code = double(str2num(cmdout)>0);
            information = in;
            scores = [];
        case "dec"
            % check input
            if any(size(in)~=[1,2^m]), error('invalid input size.'); end
            % prepare input bits
%             inStr = sprintf('%f',in);
            inStr = sprintf('%1.1f ',in); inStr(end) = [];
            % decode
            [status1, cmdout1] = system(strcat(exepathStr," -d ",rStr," ",mStr," ",LStr," ",maskStr," ",inStr));
            [scores, information, code] = parseDecOutput(cmdout1);
    end

end

function [scores, information, code] = parseDecOutput(str)
    s = split(str, newline);
    sz = (length(s) - 1) / 3;
    k = length(str2num(s{2}));
    n = length(str2num(s{3}));
    scores = zeros(sz, 1);
    information = zeros(sz, k);
    code = zeros(sz, n);
    for i = 1:sz
        scores(i) = str2double(s{1+3*(i-1)});
        information(i,:) = str2num(s{2 + 3*(i-1)});
        code(i, :) = str2num(s{3 + 3*(i-1)});
    end
end

function n = numnodes(r,m)
    if r==m || r==0
        n=1;
    else
        n = numnodes(r-1,m-1)+numnodes(r,m-1);
    end
end
