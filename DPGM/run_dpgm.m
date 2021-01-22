tic
addpath('/home/godl/IJCAI_CODING/Algorithm/AAAI_2018/EPM-master');
N = 1005;%输入节点个数
T = 10;%输入时间

%for opo=1:25
A=zeros(1005,1005,10);
disp('-----------1----------------')
zzz=0;
for o=1:10
%kk=o+2*(zzz);
k=importdata(['/home/godl/IJCAI_CODING/Dataset/email-Eu-core/email.t0',num2str(o),'.edges']);
G=graph(k(:,1),k(:,2));
A(:,:,o)=adjacency(G);
end
clear G k o;
disp('----------2-------------')
B = cell(1, T); %定义元跑数组
for t=1:T
    B{t} = sparse(A(:, :, t));
end
disp('----------3-------------')
clear A
options.burnin = 3;
options.mcmcsamps = 3;
options.display = 1;

Acore = B;
TrainRatio = 1; 
numslices = T;
disp('----------4-------------')
idx_train = cell(numslices, 1); idx_test = cell(numslices, 1); BTrain_Mask = cell(numslices, 1);
%figure(1);
for t = 1:numslices
[idx_train{t}, idx_test{t}] = Create_Mask_big_network(Acore{t}, TrainRatio);
end
disp('----------5-------------')

options.TrainRatio = TrainRatio;
options.idx_train = idx_train;
options.idx_test = idx_test;

numrep = 1;result = cell(numrep, 1);
disp('----------6-------------')
options.Mex = 1;
options.eval=1;
clear Acore BTrain_Mask idx_train idx_test
clear t T TrainRatio numslices 
for rep = 1:numrep
    result{rep} = HGPDR_batch_Gibbs(B, options);
end
disp('----------7-------------')
%最后使用下面的聚类
[a1,b1]=kmeans(result{1,1}.ProbAve{1,1},20);
[a2,b2]=kmeans(result{1,1}.ProbAve{2,1},20);
[a3,b3]=kmeans(result{1,1}.ProbAve{3,1},20);
[a4,b4]=kmeans(result{1,1}.ProbAve{4,1},20);
[a5,b5]=kmeans(result{1,1}.ProbAve{5,1},20);
[a6,b6]=kmeans(result{1,1}.ProbAve{6,1},20);
[a7,b7]=kmeans(result{1,1}.ProbAve{7,1},20);
[a8,b8]=kmeans(result{1,1}.ProbAve{8,1},20);
[a9,b9]=kmeans(result{1,1}.ProbAve{9,1},20);
[a10,b10]=kmeans(result{1,1}.ProbAve{10,1},20);

disp('----------8-------------')
real1=importdata('/home/godl/IJCAI_CODING/Dataset/email-Eu-core/real_lab');

disp('----------9-------------')
Z=zeros(1,1005);
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


NMI(real1,Z(1,:))
NMI(real1,Z(2,:))
NMI(real1,Z(3,:))
NMI(real1,Z(4,:))
NMI(real1,Z(5,:))
NMI(real1,Z(5,:))
NMI(real1,Z(5,:))
NMI(real1,Z(5,:))
NMI(real1,Z(5,:))
NMI(real1,Z(5,:))

%0.873,
disp('----------10-------------')
%end


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
toc