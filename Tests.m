% Test script for RM functions
clear
AvgRunTime = struct('rmenc_old',0,'rmenc',0,'rmenc_v2',0,'rmdemap',0,'rmdec_fht',0,...
    'rmlistdec_fht',0,'rmdec_reed',0,'rmdec_dumer',0,'rmdumerexe_enc',0,'rmdumerexe_dec',0,...
    'rmdec_rpa',0,'rmdec_listrpa',0);

ErrorsFlag = 2; % 0 - no channel \ 1 - binary errors \ 2 - erasure channel \ 3 - errors and erasures
MaxIterNum = 1e3;

TestEncoders = 0;
TestDemap = 0;
TestFHT = 0;
TestListDecFHT = 0;
TestML = 0;
TestReed = 0;
TestDumer = 0;
TestListDecDumer = 0;
TestRPA = 1;
TestListRPA = 1;

for m = 3 : 7
    for r = 1 : 2%m
        fprintf('Working on RM(%d,%d)...',r,m);
        
        kRM = 0; for rr = 0:1:r, kRM=kRM+nchoosek(m,rr); end
        nRM = 2.^m;
        
        for msgdec = 0 : 1 : 2^kRM-1
            message = de2bi(msgdec,kRM,'left-msb');
            if msgdec > MaxIterNum
                warning(sprintf('Reached %d codewords. Proceeding to next prarameters.',MaxIterNum));
                break;
            end
            
            if TestListDecDumer && r>=2 && r<m
                tic
                Codeword_dumerexe = rmlistdec_dumerexe(message,r,m,"enc");
                AvgRunTime.rmdumerexe_enc = AvgRunTime.rmdumerexe_enc + toc;
                if ErrorsFlag
                    Codeword = Codeword_dumerexe;
                end
            else
                tic
                Codeword = rmenc_v2(message,r,m);
                AvgRunTime.rmenc_v2 = AvgRunTime.rmenc_v2 + toc;
            end
            
            if TestEncoders
                tic
                [Codeword_old,G_old,monomG_old] = rmenc_old(message,r,m);
                AvgRunTime.rmenc_old = AvgRunTime.rmenc_old + toc;
                
                tic
                [Codeword,G,monomG] = rmenc(message,r,m);
                AvgRunTime.rmenc = AvgRunTime.rmenc + toc;
                
                tic
                Codeword_v2 = rmenc_v2(message,r,m);
                AvgRunTime.rmenc_v2 = AvgRunTime.rmenc_v2 + toc;
                
                if ~isequal(Codeword,Codeword_old,Codeword_v2), error('encoder failed.'); end
            end
                        
            if TestDemap
                tic
                decodedMessage = rmdemap(Codeword,r,m);
                AvgRunTime.rmdemap = AvgRunTime.rmdemap + toc;
                if ~isequal(message,decodedMessage), error('demapper failed.'); end
            end
            
            ChannelIn = Codeword;
            switch ErrorsFlag
                case 0 
                    ChannelOut = ChannelIn;
                case 1
                    dRM = 2^(m-r);
                    if dRM==1, t = 0; else, t = randi([0 dRM/2-1],[1 1]); end
                    tIdxs = randi([1 nRM],[1 t]);
                    ChannelOut = ChannelIn;
                    ChannelOut(tIdxs) = 1-ChannelOut(tIdxs);
                case 2
                    dRM = 2^(m-r);
                    if dRM==1, t = 0; else, t = randi([0 dRM-1],[1 1]); end
                    tIdxs = randi([1 nRM],[1 t]);
                    ChannelOut = ChannelIn;
                    ChannelOut(tIdxs) = 1/2;
                case 3
                    dRM = 2^(m-r);
                    if dRM==1, tErr = 0; else, tErr = randi([0 dRM/2-1],[1 1]); end
                    tErrIdxs = randi([1 nRM],[1 tErr]);
                    if dRM==1, tErsr = 0; else, tErsr = randi([0 dRM-1-2*tErr],[1 1]); end
                    tErsrIdxs = randi([1 nRM],[1 tErsr]);
                    ChannelOut = ChannelIn;
                    ChannelOut(tErrIdxs) = 1-ChannelOut(tErrIdxs);
                    ChannelOut(tErsrIdxs) = 1/2;
            end
            ChannelOut = 1-2*ChannelOut;
            
            if r==1
                if TestFHT
                    tic
                    [decodedCodeword,decodedMessage] = rmdec_fht(ChannelOut,r,m);
                    AvgRunTime.rmdec_fht = AvgRunTime.rmdec_fht + toc;
                    if ~isequal(decodedCodeword,Codeword), error('FHT decoder failed (codeword).'); end
                    if ~isequal(decodedMessage,message), error('FHT decoder failed (message).'); end
                end
                
                if TestListDecFHT
                    tic
                    [decodedCodewordList,decodedMessageList] = rmlistdec_fht(ChannelOut,r,m);
                    AvgRunTime.rmlistdec_fht = AvgRunTime.rmlistdec_fht + toc;
                    ListLocation = find(all(decodedCodewordList==Codeword,2));
                    if ~isequal(ListLocation,1), error('FHT list-decoder failed (codeword).'); end
                    if ~isequal(decodedMessageList(1,:),message), error('FHT decoder failed (message).'); end
                end
            end
            
            if r==m
                if TestML
                    decodedCodeword = double(ChannelOut<0); % ML decoder
                    decodedMessage = rmdemap(decodedCodeword,r,m);
                    if ~isequal(decodedCodeword,Codeword), error('ML decoder failed (codeword).'); end
                    if ~isequal(decodedMessage,message), error('ML decoder failed (message).'); end
                end
            end
                        
            if TestReed
                tic
                [decodedCodeword,decodedMessage] = rmdec_reed((1-ChannelOut)/2,r,m);
                AvgRunTime.rmdec_reed = AvgRunTime.rmdec_reed + toc;
                if ~isequal(decodedCodeword,Codeword)
                    fprintf('Reed''s decoder failed (codeword).\n'); 
                    fprintf('codeword: %s.\n',num2str(Codeword));
                    fprintf('decoded : %s.\n',num2str(decodedCodeword));
                    tIdxs
                end
                if ~isequal(decodedMessage,message)
                    fprintf('Reed''s decoder failed (r=%d,m=%d,message).\n',r,m); 
                    fprintf('message: %s.\n',num2str(message));
                    fprintf('decoded: %s.\n',num2str(decodedMessage));
                end
            end
            
            if TestDumer
                tic
                [decodedCodeword,decodedMessage] = rmdec_dumer(ChannelOut,r,m,mask); % <- [av,au] doewnt work, used demap instead
                AvgRunTime.rmdec_dumer = AvgRunTime.rmdec_dumer + toc;
                if ~isequal(decodedCodeword,Codeword), error('Dumer''s decoder failed (codeword).'); end
                if ~isequal(decodedMessage,message), error('Dumer''s decoder failed (message).'); end
            end
            
            if TestListDecDumer && r>=2 && r<m
                L = min([32,2^kRM]);
                tic;
                if ErrorsFlag
                    [code,information,scores] = rmlistdec_dumerexe(-1*ChannelOut,r,m,"dec",L);
                else
                    [code,information,scores] = rmlistdec_dumerexe(-1+2*Codeword_dumerexe,r,m,"dec",L);
                end
                decodedCodeword = (code(1,:)+1)/2;
                decodedMessage = information(1,:);
                AvgRunTime.rmdumerexe_dec = AvgRunTime.rmdumerexe_dec + toc;
                if ~isequal(decodedCodeword,Codeword_dumerexe), error('Dumer''s list-decoder failed (codeword).'); end
                if ~isequal(decodedMessage,message), error('Dumer''s list-decoder failed (message).'); end
            end
            
            if TestRPA
                Nmax = 10; theta = 0.05;
                tic
                [decodedCodeword, ~, ~] = rmdec_rpa(ChannelOut, m, r, Nmax, theta);
                AvgRunTime.rmdec_rpa = AvgRunTime.rmdec_rpa + toc;
                if ~isequal(decodedCodeword,Codeword), error('RPA decoder failed (codeword).'); end
            end
            
            if TestListRPA
                Nmax = 10; theta = 0.05; t = 3;
                tic
                [decodedCodeword, ~] = rmlistdec_rpa(ChannelOut, m, r, Nmax, theta, t);
                decodedCodeword = decodedCodeword(1,:);
                AvgRunTime.rmdec_listrpa = AvgRunTime.rmdec_listrpa + toc;
                if ~isequal(decodedCodeword,Codeword), error('List-RPA decoder failed (codeword).'); end
            end
            
        end
        
        fprintf('\t Done.\n');
    end
end
