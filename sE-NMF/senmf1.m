%%sE-NMF
%code by Dongyuan Li 2019/12/23 
function [Bt,Ft,error]=senmf1(A,k,iter)

[n,~]=size(A);
Bt=abs(rand(n,k));
Ft=abs(rand(k,n));
error=zeros(iter,1);
for o=1:iter
Bt=Bt.*((A*Ft')./(Bt*(Ft*Ft')+eps));
Ft=Ft.*((Bt'*A)./(Bt'*Bt*Ft+eps));
error(o,1)=norm((A-Bt*Ft),'fro');
end

