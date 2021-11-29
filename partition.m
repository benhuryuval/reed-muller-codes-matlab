function b = partition(s,k)
    
    len = length(s);
    
    st = 1;
    for ii = 1 : len
        if s(ii)<s(len), continue; end
        tmp = s(ii); s(ii) = s(st); s(st) = tmp; % swap a(ii) and a(st)
        st = st + 1;
    end
    tmp = s(len); s(len) = s(st); s(st) = tmp; % swap a(end) and a(st)
    if k==st, b=s; return; end
    if st > k,  b = partition(s(1:st),k);
    else,       b = partition(s(st+1:end),k-st);
    end