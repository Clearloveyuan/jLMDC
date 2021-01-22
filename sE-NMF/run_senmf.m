tic
disp('group1')
zzz=1;
A=zeros(10000,10000,10);
for o=1:10
k=importdata(['/home/godl/IJCAI_CODING/Dataset/Green/birthdeath/birthdeath.t0',num2str(o),'.edges']);
G=graph(k(:,1),k(:,2));
A(:,:,o)=adjacency(G);
end
clear G k o;

for i=1:10
p=dlmread(['/home/godl/IJCAI_CODING/Dataset/Green/birthdeath/birthdeath.t0',num2str(i),'.comm']);
[m,n]=size(p);
%p=p+1;
mm=zeros(10000,1);

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


[B1,F1]=senmf1(A1,50);
clear A1 F1 
[B2,F2]=senmf1(A2,50);
clear A2 F2 
[B3,F3]=senmf1(A3,50);
clear A3 F3 
[B4,F4]=senmf1(A4,50);
clear A4 F4 
[B5,F5]=senmf1(A5,50);
clear A5 F5 
[B6,F6]=senmf1(A6,50);
clear A6 F6 
[B7,F7]=senmf1(A7,50);
clear A7 F7 
[B8,F8]=senmf1(A8,50);
clear A8 F8 
[B9,F9]=senmf1(A9,50);
clear A9 F9 
[B10,F10]=senmf1(A10,50);
clear A10 F10
%[a1,b1]=kmeans(B1,4,'Replicates',100);
%[a2,b2]=kmeans(B2,4,'Replicates',100);
%[a3,b3]=kmeans(B3,4,'Replicates',100);
%[a4,b4]=kmeans(B4,4,'Replicates',100);
%[a5,b5]=kmeans(B5,4,'Replicates',100);
%[a6,b6]=kmeans(B6,4,'Replicates',100);
%[a7,b7]=kmeans(B7,4,'Replicates',100);
%[a8,b8]=kmeans(B8,4,'Replicates',100);
%[a9,b9]=kmeans(B9,4,'Replicates',100);
%[a10,b10]=kmeans(B10,4,'Replicates',100);





[a1,b1]=kmeans(B1,4,'Replicates',30);
[a2,b2]=kmeans(B2,4,'Replicates',30);
[a3,b3]=kmeans(B3,4,'Replicates',30);
[a4,b4]=kmeans(B4,4,'Replicates',30);
[a5,b5]=kmeans(B5,4,'Replicates',30);
[a6,b6]=kmeans(B6,4,'Replicates',30);
[a7,b7]=kmeans(B7,4,'Replicates',30);
[a8,b8]=kmeans(B8,4,'Replicates',30);
[a9,b9]=kmeans(B9,4,'Replicates',30);
[a10,b10]=kmeans(B10,4,'Replicates',30);






NMI(real1,a1)
NMI(real2,a2)
NMI(real3,a3)
NMI(real4,a4)
NMI(real5,a5)
NMI(real6,a6)
NMI(real7,a7)
NMI(real8,a8)
NMI(real9,a9)
NMI(real10,a10)
save(['/home/godl/IJCAI_CODING/Dataset/Green/se_NMF_f1/group',num2str(zzz),'.mat']);
clear 




disp('group2')
zzz=2;
A=zeros(10000,10000,10);
for o=1:10
kk=o+10;
k=importdata(['/home/godl/IJCAI_CODING/Dataset/Green/birthdeath/birthdeath.t0',num2str(kk),'.edges']);
G=graph(k(:,1),k(:,2));
A(:,:,o)=adjacency(G);
end
clear G k o;

for i=1:10
kkk=i+10;
p=dlmread(['/home/godl/IJCAI_CODING/Dataset/Green/birthdeath/birthdeath.t0',num2str(kkk),'.comm']);
[m,n]=size(p);
%p=p+1;
mm=zeros(10000,1);

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


[B1,F1]=senmf1(A1,50);
clear A1 F1 
[B2,F2]=senmf1(A2,50);
clear A2 F2 
[B3,F3]=senmf1(A3,50);
clear A3 F3 
[B4,F4]=senmf1(A4,50);
clear A4 F4 
[B5,F5]=senmf1(A5,50);
clear A5 F5 
[B6,F6]=senmf1(A6,50);
clear A6 F6 
[B7,F7]=senmf1(A7,50);
clear A7 F7 
[B8,F8]=senmf1(A8,50);
clear A8 F8 
[B9,F9]=senmf1(A9,50);
clear A9 F9 
[B10,F10]=senmf1(A10,50);
clear A10 F10
%[a1,b1]=kmeans(B1,4,'Replicates',100);
%[a2,b2]=kmeans(B2,4,'Replicates',100);
%[a3,b3]=kmeans(B3,4,'Replicates',100);
%[a4,b4]=kmeans(B4,4,'Replicates',100);
%[a5,b5]=kmeans(B5,4,'Replicates',100);
%[a6,b6]=kmeans(B6,4,'Replicates',100);
%[a7,b7]=kmeans(B7,4,'Replicates',100);
%[a8,b8]=kmeans(B8,4,'Replicates',100);
%[a9,b9]=kmeans(B9,4,'Replicates',100);
%[a10,b10]=kmeans(B10,4,'Replicates',100);





[a1,b1]=kmeans(B1,4,'Replicates',30);
[a2,b2]=kmeans(B2,4,'Replicates',30);
[a3,b3]=kmeans(B3,4,'Replicates',30);
[a4,b4]=kmeans(B4,4,'Replicates',30);
[a5,b5]=kmeans(B5,4,'Replicates',30);
[a6,b6]=kmeans(B6,4,'Replicates',30);
[a7,b7]=kmeans(B7,4,'Replicates',30);
[a8,b8]=kmeans(B8,4,'Replicates',30);
[a9,b9]=kmeans(B9,4,'Replicates',30);
[a10,b10]=kmeans(B10,4,'Replicates',30);






NMI(real1,a1)
NMI(real2,a2)
NMI(real3,a3)
NMI(real4,a4)
NMI(real5,a5)
NMI(real6,a6)
NMI(real7,a7)
NMI(real8,a8)
NMI(real9,a9)
NMI(real10,a10)
save(['/home/godl/IJCAI_CODING/Dataset/Green/se_NMF_f1/group',num2str(zzz),'.mat']);
clear 






disp('group3')
zzz=3;
A=zeros(10000,10000,10);
for o=1:10
    kk=o+20;
k=importdata(['/home/godl/IJCAI_CODING/Dataset/Green/birthdeath/birthdeath.t0',num2str(kk),'.edges']);
G=graph(k(:,1),k(:,2));
A(:,:,o)=adjacency(G);
end
clear G k o;

for i=1:10
    kkk=i+20;
p=dlmread(['/home/godl/IJCAI_CODING/Dataset/Green/birthdeath/birthdeath.t0',num2str(kkk),'.comm']);
[m,n]=size(p);
%p=p+1;
mm=zeros(10000,1);

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


[B1,F1]=senmf1(A1,50);
clear A1 F1 
[B2,F2]=senmf1(A2,50);
clear A2 F2 
[B3,F3]=senmf1(A3,50);
clear A3 F3 
[B4,F4]=senmf1(A4,50);
clear A4 F4 
[B5,F5]=senmf1(A5,50);
clear A5 F5 
[B6,F6]=senmf1(A6,50);
clear A6 F6 
[B7,F7]=senmf1(A7,50);
clear A7 F7 
[B8,F8]=senmf1(A8,50);
clear A8 F8 
[B9,F9]=senmf1(A9,50);
clear A9 F9 
[B10,F10]=senmf1(A10,50);
clear A10 F10

%[a1,b1]=kmeans(B1,4,'Replicates',100);
%[a2,b2]=kmeans(B2,4,'Replicates',100);
%[a3,b3]=kmeans(B3,4,'Replicates',100);
%[a4,b4]=kmeans(B4,4,'Replicates',100);
%[a5,b5]=kmeans(B5,4,'Replicates',100);
%[a6,b6]=kmeans(B6,4,'Replicates',100);
%[a7,b7]=kmeans(B7,4,'Replicates',100);
%[a8,b8]=kmeans(B8,4,'Replicates',100);
%[a9,b9]=kmeans(B9,4,'Replicates',100);
%[a10,b10]=kmeans(B10,4,'Replicates',100);





[a1,b1]=kmeans(B1,4,'Replicates',30);
[a2,b2]=kmeans(B2,4,'Replicates',30);
[a3,b3]=kmeans(B3,4,'Replicates',30);
[a4,b4]=kmeans(B4,4,'Replicates',30);
[a5,b5]=kmeans(B5,4,'Replicates',30);
[a6,b6]=kmeans(B6,4,'Replicates',30);
[a7,b7]=kmeans(B7,4,'Replicates',30);
[a8,b8]=kmeans(B8,4,'Replicates',30);
[a9,b9]=kmeans(B9,4,'Replicates',30);
[a10,b10]=kmeans(B10,4,'Replicates',30);






NMI(real1,a1)
NMI(real2,a2)
NMI(real3,a3)
NMI(real4,a4)
NMI(real5,a5)
NMI(real6,a6)
NMI(real7,a7)
NMI(real8,a8)
NMI(real9,a9)
NMI(real10,a10)
save(['/home/godl/IJCAI_CODING/Dataset/Green/se_NMF_f1/group',num2str(zzz),'.mat']);
clear 





disp('group4')
zzz=4;
A=zeros(10000,10000,10);
for o=1:10
    kk=o+30;
k=importdata(['/home/godl/IJCAI_CODING/Dataset/Green/birthdeath/birthdeath.t0',num2str(kk),'.edges']);
G=graph(k(:,1),k(:,2));
A(:,:,o)=adjacency(G);
end
clear G k o;

for i=1:10
    kkk=i+30;
p=dlmread(['/home/godl/IJCAI_CODING/Dataset/Green/birthdeath/birthdeath.t0',num2str(kkk),'.comm']);
[m,n]=size(p);
%p=p+1;
mm=zeros(10000,1);

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


[B1,F1]=senmf1(A1,50);
clear A1 F1 
[B2,F2]=senmf1(A2,50);
clear A2 F2 
[B3,F3]=senmf1(A3,50);
clear A3 F3 
[B4,F4]=senmf1(A4,50);
clear A4 F4 
[B5,F5]=senmf1(A5,50);
clear A5 F5 
[B6,F6]=senmf1(A6,50);
clear A6 F6 
[B7,F7]=senmf1(A7,50);
clear A7 F7 
[B8,F8]=senmf1(A8,50);
clear A8 F8 
[B9,F9]=senmf1(A9,50);
clear A9 F9 
[B10,F10]=senmf1(A10,50);
clear A10 F10

%[a1,b1]=kmeans(B1,4,'Replicates',100);
%[a2,b2]=kmeans(B2,4,'Replicates',100);
%[a3,b3]=kmeans(B3,4,'Replicates',100);
%[a4,b4]=kmeans(B4,4,'Replicates',100);
%[a5,b5]=kmeans(B5,4,'Replicates',100);
%[a6,b6]=kmeans(B6,4,'Replicates',100);
%[a7,b7]=kmeans(B7,4,'Replicates',100);
%[a8,b8]=kmeans(B8,4,'Replicates',100);
%[a9,b9]=kmeans(B9,4,'Replicates',100);
%[a10,b10]=kmeans(B10,4,'Replicates',100);





[a1,b1]=kmeans(B1,4,'Replicates',30);
[a2,b2]=kmeans(B2,4,'Replicates',30);
[a3,b3]=kmeans(B3,4,'Replicates',30);
[a4,b4]=kmeans(B4,4,'Replicates',30);
[a5,b5]=kmeans(B5,4,'Replicates',30);
[a6,b6]=kmeans(B6,4,'Replicates',30);
[a7,b7]=kmeans(B7,4,'Replicates',30);
[a8,b8]=kmeans(B8,4,'Replicates',30);
[a9,b9]=kmeans(B9,4,'Replicates',30);
[a10,b10]=kmeans(B10,4,'Replicates',30);






NMI(real1,a1)
NMI(real2,a2)
NMI(real3,a3)
NMI(real4,a4)
NMI(real5,a5)
NMI(real6,a6)
NMI(real7,a7)
NMI(real8,a8)
NMI(real9,a9)
NMI(real10,a10)
save(['/home/godl/IJCAI_CODING/Dataset/Green/se_NMF_f1/group',num2str(zzz),'.mat']);
clear 





disp('group5')
zzz=5;
A=zeros(10000,10000,10);
for o=1:10
    kk=o+40;
k=importdata(['/home/godl/IJCAI_CODING/Dataset/Green/birthdeath/birthdeath.t0',num2str(kk),'.edges']);
G=graph(k(:,1),k(:,2));
A(:,:,o)=adjacency(G);
end
clear G k o;

for i=1:10
    kkk=i+40;
p=dlmread(['/home/godl/IJCAI_CODING/Dataset/Green/birthdeath/birthdeath.t0',num2str(kkk),'.comm']);
[m,n]=size(p);
%p=p+1;
mm=zeros(10000,1);

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


[B1,F1]=senmf1(A1,50);
clear A1 F1 
[B2,F2]=senmf1(A2,50);
clear A2 F2 
[B3,F3]=senmf1(A3,50);
clear A3 F3 
[B4,F4]=senmf1(A4,50);
clear A4 F4 
[B5,F5]=senmf1(A5,50);
clear A5 F5 
[B6,F6]=senmf1(A6,50);
clear A6 F6 
[B7,F7]=senmf1(A7,50);
clear A7 F7 
[B8,F8]=senmf1(A8,50);
clear A8 F8 
[B9,F9]=senmf1(A9,50);
clear A9 F9 
[B10,F10]=senmf1(A10,50);
clear A10 F10

%[a1,b1]=kmeans(B1,4,'Replicates',100);
%[a2,b2]=kmeans(B2,4,'Replicates',100);
%[a3,b3]=kmeans(B3,4,'Replicates',100);
%[a4,b4]=kmeans(B4,4,'Replicates',100);
%[a5,b5]=kmeans(B5,4,'Replicates',100);
%[a6,b6]=kmeans(B6,4,'Replicates',100);
%[a7,b7]=kmeans(B7,4,'Replicates',100);
%[a8,b8]=kmeans(B8,4,'Replicates',100);
%[a9,b9]=kmeans(B9,4,'Replicates',100);
%[a10,b10]=kmeans(B10,4,'Replicates',100);





[a1,b1]=kmeans(B1,4,'Replicates',30);
[a2,b2]=kmeans(B2,4,'Replicates',30);
[a3,b3]=kmeans(B3,4,'Replicates',30);
[a4,b4]=kmeans(B4,4,'Replicates',30);
[a5,b5]=kmeans(B5,4,'Replicates',30);
[a6,b6]=kmeans(B6,4,'Replicates',30);
[a7,b7]=kmeans(B7,4,'Replicates',30);
[a8,b8]=kmeans(B8,4,'Replicates',30);
[a9,b9]=kmeans(B9,4,'Replicates',30);
[a10,b10]=kmeans(B10,4,'Replicates',30);






NMI(real1,a1)
NMI(real2,a2)
NMI(real3,a3)
NMI(real4,a4)
NMI(real5,a5)
NMI(real6,a6)
NMI(real7,a7)
NMI(real8,a8)
NMI(real9,a9)
NMI(real10,a10)
save(['/home/godl/IJCAI_CODING/Dataset/Green/se_NMF_f1/group',num2str(zzz),'.mat']);
clear 



toc






%{
Z=zeros(10,10000);
Z(1,:)=a1';
Z(2,:)=a2';
Z(3,:)=a3';
Z(4,:)=a4';
Z(5,:)=a5';
Z(6,:)=a6';
Z(7,:)=a7';
Z(8,:)=a8';
Z(9,:)=a9';
Z(10,:)=a10';


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
dlmwrite((['/home/godl/IJCAI_CODING/Algorithm/PisCES-master/syn_NF1/a',num2str(i)]),a,'precision',5,'delimiter',' ');
end
toc
%}
