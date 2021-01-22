function [out] = CRT_matrix(x,r)
%Mingyuan Zhou, Duke ECE, May, 2012
%x~NegBin(r,p)
%Finding the latent L
if length(r)==1
    x  = x(:);
    [xx,ii,jj] = unique(x);
    L = zeros(size(x));
    Lsum = 0;
    if ~isempty(x)
        for i=1:length(xx)
            y = xx(i);
            if y>0
                L(jj==i)  = sum( bsxfun(@le, rand(nnz(jj==i),y), r./(r+(1:y)-1) ), 2);
            end
        end
        Lsum = sum(L);
    end
else
    [h w] = size(x);
    x = x(:);
    [xx,ii,jj] = unique(x);
    r = r(:);       
    % Lsum=0;
    L = zeros(size(x));
    for i=1:length(xx)
        y = xx(i);
        L(jj==i) = sum( bsxfun(@le, rand(nnz(jj==i),y), 1./(1+bsxfun(@rdivide,(0:y-1), r(jj==i)) )  ), 2);
    end
    % Lsum = sum(L);
    
    out = reshape(L, [h w]);
end