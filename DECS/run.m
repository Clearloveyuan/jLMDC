%% Save dataset
%%GT_Cube 1x10 cell 128x4 128x4 .......每一列对应类标,为哪一个类就为1,其余为0
%%GT_Matrixs 128x10 每一列为一个时刻的类标签
%%W_Cube 1x10 cell 128x128 128x128 为每一个时刻的邻接矩阵
tic
GT_Cube=cell(1,10);
GT_Matrix=zeros(1005,10);
W_Cube=cell(1,10);

%DES_coding
%一共需要三个工作区分别是W_Cube用来存放数据1x10cell 
%GT_Matrix nxT 每一列是一个时刻的聚类结果
%GT_Cube 1x10cell 每一个cell里面是nxT的结果，其中在哪一类，哪一类数据就全部为1，其余的全部为0


%%导入真实数据的邻接矩阵
disp('group1')
zzz=1;
A=cell(1,10);
for o=1:10
k=importdata(['/home/godl/IJCAI_CODING/Dataset/email-Eu-core/email.t0',num2str(o),'.edges']);
G=graph(k(:,1),k(:,2));
A{1,o}=adjacency(G);
end
clear G k o;
W_Cube=A;
disp('W_Cube is ready')
dynMoeaResult=cell(1,10);
%%%获取真实标签的方法:
for i=1:10
p=importdata('/home/godl/IJCAI_CODING/Dataset/email-Eu-core/real_lab');
[m,n]=size(p);
%p=p+1;

mm=zeros(1005,1);

for o=1:m
for k=1:n
if p(o,k)~=0
t=p(o,k);
mm(t,1)=o;
end
end
end
mm(mm==0)=1;
lab=zeros(1005,49);
for kkk=1:1005
loca=mm(kkk,1);
lab(kkk,loca)=1;
end
dynMoeaResult{1,i}=lab;
clear lab
end 
clear i k m mm n o t p;
GT_Cube=dynMoeaResult;
clear loca kkk dynMoeaResult A
disp('GT_Cube is ready')


%%%GT_Matrix
B=zeros(1005,10);
for i=1:10
p=importdata('/home/godl/IJCAI_CODING/Dataset/email-Eu-core/real_lab');
[m,n]=size(p);
%p=p+1;

mm=zeros(1005,1);

for o=1:m
for k=1:n
if p(o,k)~=0
t=p(o,k);
mm(t,1)=o;
end
end
end
mm(mm==0)=1;
B(:,i)=mm;
end 
clear i k m mm n o t p;
GT_Matrix=B;
clear B
disp('GT_Matrix is ready')
save(['/home/godl/IJCAI_CODING/Algorithm/CommunityDetection/CommunityDetection-master/DECS/datasets/email',num2str(zzz),'.mat']);
disp('data is preseving')



flag = 1;
load(['datasets/email',num2str(zzz),'.mat']);
%% Parameter setting
maxgen = 2;         % the maximum number of iterations
pop_size = 100;       % the population size  
num_neighbor = 5;     % the neighbor size for each subproblem in decomposition-based multi-objective optimization
p_mutation = 0.20;    % the mutation rate
p_migration = 0.50;   % the migration rate
p_mu_mi = 0.50;       % the paramater to organize the execution of mutation and migration
PGLP_iter = 5;        % the number of iterations in PGLP
num_repeat = 1;       % the number of repeated run

%% Results at each time step
dynMod = [];          % modularity of detected community structure
dynNmi = [];          % NMI between detected community structure and the ground truth
dynPop = {};          % the population
dynTime = [];         % the running time
DECS_Result = {}; % the detected community structure
disp('action')
for r = 1 : num_repeat
%     global idealp weights neighbors;
    % idealp is reference point (z1, z2) where z1 and z2
    % are the maximum of the 1st and 2nd objective functions
    num_timestep = size(W_Cube, 2);  % W_Cube contains several cells restoring temporal adjacent matrices
    %% DECS only optimizes the modularity at the 1st time step
    timestep_num = 1;
    disp('doing the slow thing')
    [dynMod(1,r), dynPop{1,r}, DECS_Result{1,r}, dynTime(1,r)] = ...
        DECS_1(W_Cube{timestep_num}, maxgen, pop_size, p_mutation, p_migration, p_mu_mi, PGLP_iter); 
    % calculate NMI for synthetic or real-world networks
    disp('1 is doing')
    if flag == 1
        % for synthetic networks
        dynNmi(1,r) = NMI(GT_Matrix(:,1)',DECS_Result{1,r}); 
    else 
        % for real-world networks
        dynNmi(1,r) = NMI(GT_Cube{timestep_num},DECS_Result{1,r});
    end
    disp(['timestep = ', num2str(timestep_num), ', Modularity = ',...
        num2str(dynMod(timestep_num,r)), ', NMI = ', num2str(dynNmi(timestep_num,r))]);
    disp('2 is doing')
    %% DECS optimizes the modularity and NMI in the following time steps
    for timestep_num = 2 : num_timestep  
        [dynMod(timestep_num,r), dynPop{timestep_num,r}, DECS_Result{timestep_num,r}, ... 
            dynTime(timestep_num,r)] = DECS_2(W_Cube{timestep_num}, maxgen, pop_size, ...
            p_mutation, p_migration, p_mu_mi, num_neighbor, DECS_Result{timestep_num-1,r}, PGLP_iter);
        
        if flag == 1
            dynNmi(timestep_num,r) = NMI(DECS_Result{timestep_num,r}, GT_Matrix(:,timestep_num)');
        else
            dynNmi(timestep_num,r) = NMI(DECS_Result{timestep_num,r}, GT_Cube{timestep_num});
        end
        
            disp(['timestep = ', num2str(timestep_num), ', Modularity = ',...
        num2str(dynMod(timestep_num,r)), ', NMI = ', num2str(dynNmi(timestep_num,r))]);
    end
end

avg_dynMod = sum(dynMod,2)/num_repeat;
avg_dynNmi = sum(dynNmi,2)/num_repeat;
avg_dynMod = sum(dynMod,2)/num_repeat;


Z=zeros(10,1005);
Z(1,:)=DECS_Result{1,1};
Z(2,:)=DECS_Result{2,1};
Z(3,:)=DECS_Result{3,1};
Z(4,:)=DECS_Result{4,1};
Z(5,:)=DECS_Result{5,1};
Z(6,:)=DECS_Result{6,1};
Z(7,:)=DECS_Result{7,1};
Z(8,:)=DECS_Result{8,1};
Z(9,:)=DECS_Result{9,1};
Z(10,:)=DECS_Result{10,1};


%{
for i=1:5
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
dlmwrite((['/home/godl/IJCAI_CODING/Algorithm/CommunityDetection/CommunityDetection-master/NF1/a',num2str(i)]),a,'precision',5,'delimiter',' ');
end
%}

save(['/home/godl/IJCAI_CODING/Algorithm/CommunityDetection/CommunityDetection-master/NF1/email',num2str(zzz),'.mat']);

toc





