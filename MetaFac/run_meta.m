%{
A=zeros(393,393,10);
for o=1:10
k=importdata(['/home/godl/IJCAI_CODING/Dataset/GN/syn-var-z=5/z_5/synvar_5.t0',num2str(o),'.edges']);
G=graph(k(:,1),k(:,2));
A(:,:,o)=adjacency(G);
end
clear G k o;




%%%获取真实标签的方法:
for i=1:10
p=dlmread(['/home/godl/IJCAI_CODING/Dataset/GN/syn-var-z=5/z_5/synvar_5.t0',num2str(i),'.comm']);
[m,n]=size(p);
%p=p+1;

mm=zeros(393,1);

for o=1:m
for k=1:n
if p(o,k)~=0
t=p(o,k);
mm(t,1)=o;
end
end
end
%mm=mm';
%mm(mm==0)=1;
eval(['real',num2str(i),'=','mm',';']);
end 
clear i k m mm n o t p;


A1=A(:,:,1);
A2=A(:,:,2);
A3=A(:,:,3);
A4=A(:,:,4);
A5=A(:,:,5);
A6=A(:,:,6);
A7=A(:,:,7);
A8=A(:,:,8);
A9=A(:,:,9);
A10=A(:,:,10);
%}
tic
A1=tensor(A1);
A2=tensor(A2);
A3=tensor(A3);
A4=tensor(A4);
A5=tensor(A5);
A6=tensor(A6);
A7=tensor(A7);
A8=tensor(A8);
A9=tensor(A9);
A10=tensor(A10);



clear RG Vdims
K = 20; 
w=1;
prior={};
fprintf('\n => generate 3-order data tensor for R1');
RG{1} = {A1,1:2,w};
Vdims = [size(A1)];
Vdims,RG
fprintf('\n => run MetaFac for R1');
[S1,U1,iters,cvtime]=metafac(RG, Vdims, K, prior);

disp('-----1------')
fprintf('\n => generate 3-order data tensor for R2');
RG{2} = {A2,[2,3],w};
Vdims = [Vdims size(A2,2)];
fprintf('\n => run MetaFac for R2');
[S2,U2,iters,cvtime]=metafac(RG, Vdims, K, prior);

disp('-----2------')
fprintf('\n => generate 3-order data tensor for R3');
RG{3} = {A3,[3,4],w};
Vdims = [Vdims size(A3,2)];
fprintf('\n => run MetaFac for R3');
[S3,U3,iters,cvtime]=metafac(RG, Vdims, K, prior);


disp('-----3------')
fprintf('\n => generate 3-order data tensor for R4'); 
RG{4} = {A4,[4,5],w};
Vdims = [Vdims size(A4,2)];
[S4,U4,iters,cvtime]=metafac(RG, Vdims, K, prior);

disp('-----4------')
fprintf('\n => generate 3-order data tensor for R5');
RG{5} = {A5,[5,6],w};
Vdims = [Vdims size(A5,2)];
[S5,U5,iters,cvtime]=metafac(RG, Vdims, K, prior);


disp('-----5------')
fprintf('\n => generate 3-order data tensor for R6');
RG{6} = {A6,[6,7],w};
Vdims = [Vdims size(A6,2)];
[S6,U6,iters,cvtime]=metafac(RG, Vdims, K, prior);


disp('-----6------')
fprintf('\n => generate 3-order data tensor for R7');
RG{7} = {A7,[7,8],w};
Vdims = [Vdims size(A7,2)];
[S7,U7,iters,cvtime]=metafac(RG, Vdims, K, prior);

disp('-----7------')
fprintf('\n => generate 3-order data tensor for R8');
RG{8} = {A8,[8,9],w};
Vdims = [Vdims size(A8,2)];
[S8,U8,iters,cvtime]=metafac(RG, Vdims, K, prior);


disp('-----8------')
fprintf('\n => generate 3-order data tensor for R9');
RG{9} = {A9,[9,10],w};
Vdims = [Vdims size(A9,2)];
[S9,U9,iters,cvtime]=metafac(RG, Vdims, K, prior);


disp('-----9-----')
fprintf('\n => generate 3-order data tensor for R10');
RG{10} = {A10,[10,11],w};
Vdims = [Vdims size(A10,2)];
[S10,U10,iters,cvtime]=metafac(RG, Vdims, K, prior);



[a1,b1]=max(U1{1,2}');
[a2,b2]=max(U2{1,3}');
[a3,b3]=max(U3{1,4}');
[a4,b4]=max(U4{1,5}');
[a5,b5]=max(U5{1,6}');
[a6,b6]=max(U6{1,7}');
[a7,b7]=max(U7{1,8}');
[a8,b8]=max(U8{1,9}');
[a9,b9]=max(U9{1,10}');
[a10,b10]=max(U10{1,11}');



 Z=zeros(10,1005);
Z(1,:)=b1;
Z(2,:)=b2;
Z(3,:)=b3;
Z(4,:)=b4;
Z(5,:)=b5;
Z(6,:)=b6;
Z(7,:)=b7;
Z(8,:)=b8;
Z(9,:)=b9;
Z(10,:)=b10;


NMI(real1,b1)
NMI(real2,b2)
NMI(real3,b3)
NMI(real4,b4)
NMI(real5,b5)
NMI(real6,b6)
NMI(real7,b7)
NMI(real8,b8)
NMI(real9,b9)
NMI(real10,b10)
toc
%{
for i=1:10
clear a;
l=2;
mm=max(Z(i,:));
for o=1:mm
k=find(Z(i,:)==o);
[m,n]=size(k);
l=max(l,n);
end

a=zeros(mm,l);
for o=1:mm
k=find(Z(i,:)==o);
[m,n]=size(k);
a(o,1:n)=k;
end
dlmwrite((['/home/godl/IJCAI_CODING/Algorithm/PisCES-master/NF11/aa',num2str(i)]),a,'precision',5,'delimiter',' ');
end
%}