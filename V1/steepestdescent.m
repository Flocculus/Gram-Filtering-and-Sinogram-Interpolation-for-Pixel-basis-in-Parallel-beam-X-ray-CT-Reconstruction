function x = steepestdescent(fA, b,x, iMax)
%start

tmp = fliplr(convxh(fliplr(x),fA,[],true));




r = b - tmp;
delta = r(:)'*r(:);
delta0 = delta;
for i = 1:iMax
    
    if delta <= eps^2*delta0
        %fprintf("Done");
        break;
    end
    q = fliplr(convxh(fliplr(r),fA,[],true));
    alpha = delta/(r(:)'*q(:));
    x = x + alpha * r;
    r = r - alpha * q;
    delta = r(:)'*r(:);
    
end
