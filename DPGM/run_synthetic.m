clc;close all
state = 0;
rng('default');
% addpath(genpath('~DRL/'));
if 1
%% synthetic data
N = 65;
T = 6;
A = zeros(N, N, T);
B = cell(1, T);

bd = [0 10 30 65];
C = numel(bd) -1;
t = 1;
for cc = 1:C
    dex = bd(cc)+1:bd(cc+1);
    A(dex, dex, t) = 1;
    A(:, :, t) = A(:, :, t) - diag(diag(A(:, :, t)));
    B{t} = sparse(A(:, :, t));
end
subplot(T, 4, (t-1)* 4 + 1);imagesc(A(:,:,t));

for t = 2:T
    if t == 2
        bd = [0 15 25 35 65];
    elseif t == 3
        bd = [0 18 30 44 65];
    else
        bd = bd + [0 5 5 0];
    end
        
    C = numel(bd) -1;
    for cc = 1:C
        dex = bd(cc)+1:bd(cc+1);
        A(dex, dex, t) = 1;
        A(:, :, t) = A(:, :, t) - diag(diag(A(:, :, t)));
        B{t} = sparse(A(:, :, t));
    end
    subplot(T, 4, (t-1)* 4 + 1);imagesc(A(:,:,t));
    
    if t == 2
        bd = [0 15 35 65];
    elseif t == 3
        bd = [0 20 40 65];
    end
end

options.burnin = 3;
options.mcmcsamps = 3;
options.display = 1;

Acore = B;
TrainRatio = 0; 
numslices = T;
idx_train = cell(numslices, 1); idx_test = cell(numslices, 1); BTrain_Mask = cell(numslices, 1);
%figure(1);
for t = 1:numslices
[idx_train{t}, idx_test{t}] = Create_Mask_big_network(Acore{t}, TrainRatio);

end
options.TrainRatio = TrainRatio;
options.idx_train = idx_train;
options.idx_test = idx_test;

numrep = 1;result = cell(numrep, 1);

options.Mex = 1;
options.eval=1;
for rep = 1:numrep
    result{rep} = HGPDR_batch_Gibbs(B, options);
end
end