%%sE-NMF
%code by Dongyuan Li 2019/12/23 
function [Bt,Ft]=senmf(A,k,iter)
[n,~]=size(A);
B=rand(n,d);
F=rand(d,n);

for o=1:iter
B=B.*((2*A*F'+4*alpha*B*B'Bt)./(2*Bt*Ft*Ft'+4*alpha*Bt*Bt'*Bt+eps));
Ft=Ft.*((Bt'*A)./(Bt'*Bt*Ft+eps));
end

