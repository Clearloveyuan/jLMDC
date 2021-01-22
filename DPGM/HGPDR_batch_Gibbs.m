function [result] = HGPDR_batch_Gibbs(B, options)

%%
% hierarchical gamma process dynamic relational model
% inference algo: batch Gibbs sampling
%%

burnin = options.burnin;
mcmcsamps = options.mcmcsamps;
maxiter = burnin + mcmcsamps;

N = size(B{1}, 1);
T = numel(B);
K = 200;
Ktrace = zeros(1, maxiter);


if isfield(options,'TrainRatio')
    idx_train = options.idx_train;
    idx_test = options.idx_test;
    BTrain_Mask_triu = cell(T, 1);
    BTrain_Mask = cell(T, 1);
    BTrain = cell(T, 1);
    for t = 1:T
        BTrain_Mask_triu{t} = (zeros(N));
        BTrain_Mask_triu{t}(idx_train{t}) = 1;
        BTrain_Mask{t} = full((BTrain_Mask_triu{t} + BTrain_Mask_triu{t}'));
        BTrain{t} = triu(B{t}, 1);
        BTrain{t}(idx_test{t}) = 0;
    end
end
clear BTrain_Mask_triu BTrain_Mask
r_K = ones(K, 1)./K;
gamma0 = 1; count = 0;
c_0 = 1;
xi = 1;
beta = 1;
lambda_KK = r_K * r_K';
lambda_KK = triu(lambda_KK, 1) + triu(lambda_KK)';
lambda_KK(sparse(1:K, 1:K, true)) = xi * r_K;

phi = zeros(N, K, T);
for t = 1:T
    phi(:, :, t) = gamrnd(1,1, [N, K]);
end
phi_0 = randg(ones(N, K));
% theta_KK = zeros(K);
diagdex = sparse(1:K, 1:K, true);
triu1dex = triu(true(K), 1);
% c_n = ones(N, 1);
c_n_t = ones(N, T);

Probsamps = zeros(N, N, T);

Btriu =cell(1,T);
obs_links = cell(1,T);
dex_lef = cell(1,T);
dex_rig = cell(1,T);
for t = 1:T
    Btriu{t} = triu(BTrain{t}, 1);
    [dex_lef{t}, dex_rig{t}, obs_links{t}] = find(Btriu{t});
    idx{t} = sub2ind([N, N], dex_lef{t}, dex_rig{t});
end


latent_count = cell(1,T);

m_n_k_dot_dot = zeros(K, N, T);
diag_kk = zeros(K, K, T);
Y_n_k_t = zeros(N,K,T);
P_n_k_t = zeros(N,K,T);
P_n_k_0 = zeros(N, K);
e_0 = 1;
h_0 = 1;
g_0= 1;
f_0 = 1;
u_0 = 1;
v_0 = 1;
a_0 = 1;
b_0 = 1;
ProbAve = cell(T,1);
tic;
for iter = 1:maxiter 
    
    %% latent counts
    %% non-mex
    m_n_k_dot_dot = zeros(K, N, T);
    %%
    for t = 1:T
        if iter >options.burnin
        Btriu{t} = triu(BTrain{t}, 1);
        [dex_lef{t}, dex_rig{t}, obs_links{t}] = find(Btriu{t});
        idx{t} = sub2ind([N, N], dex_lef{t}, dex_rig{t});
        end
%     for t = 1:T
        rate = sum(phi(dex_lef{t}, :, t) * lambda_KK .* phi(dex_rig{t}, :, t), 2);
        latent_count{t} = truncated_poisson_rnd(rate);
        %% mex
        if options.Mex ==1
         [m_n_k_dot_dot(:, :, t), diag_kk_temp] = Multrnd_mik1k2j(sparse(dex_lef{t}, dex_rig{t}, latent_count{t}, N, N), phi(:,:,t), lambda_KK);
        else
        %% non-mex
        diag_kk_temp = zeros(K);
         for ij = 1:numel(idx{t})
%             pmf = phi(dex_lef{t}(ij), :, t)' * phi(dex_rig{t}(ij), :, t) .* lambda_KK;
%             pmf = pmf(:);
%             pmf = pmf./sum(pmf);
%             mij_kk = reshape(mnrnd(latent_count{t}(ij), pmf), K, K);
            %%
            pmf = (phi(dex_lef{t}(ij),:, t)'*phi(dex_rig{t}(ij),:, t)).*lambda_KK;
            mij_kk = reshape(multrnd_histc(latent_count{t}(ij),pmf(:)),K,K);
            %%
            m_n_k_dot_dot(:,dex_lef{t}(ij), t) = m_n_k_dot_dot(:,dex_lef{t}(ij), t) + sum(mij_kk,2);
            m_n_k_dot_dot(:,dex_rig{t}(ij), t) = m_n_k_dot_dot(:,dex_rig{t}(ij), t) + sum(mij_kk,1)';
            diag_kk_temp = diag_kk_temp + mij_kk + mij_kk';
        end
        end
        %%
        diag_kk_temp(sparse(1:K, 1:K, true)) = diag_kk_temp(sparse(1:K, 1:K, true))./2;
        diag_kk(:, :, t) = diag_kk_temp;
%         temp = sum(m_n_k_dot_dot(:,:, t), 2);
%         temp = temp + temp';
%         temp(sparse(1:K, 1:K, true)) = temp(sparse(1:K, 1:K, true))./2;
%         diag_kk(:, t) = temp;
    end
    Ktrace(iter) = nnz(sum(sum(m_n_k_dot_dot, 3), 2));
    % Ktrace(iter) = sum(r_K./sum(r_K) > 0.05);
    %% r_K
    theta_KK = zeros(K);
    for t = 1:T
        theta_KK = theta_KK + phi(:, :, t)' * bsxfun(@minus, sum(phi(:, :, t), 1), phi(:, :, t));
    end
    theta_KK(sparse(1:K, 1:K, true)) = theta_KK(sparse(1:K, 1:K, true))./2;
    
    m_dot_k_k_dot_dot = sum(diag_kk, 3);
    L_KK = zeros(K);
    p_kk_tilde_one_minus = zeros(K);
    temp = zeros(K, 1);
    betavec = beta .* ones(1, K);
    for k = randperm(K)
        R_K = r_K';
        R_K(k) = xi;
        p_kk_tilde_one_minus(k, :) = betavec./(betavec + theta_KK(k, :));
        L_KK(k, :) = CRT_sum_mex_matrix(sparse(m_dot_k_k_dot_dot(k,:)), r_K(k) * R_K);
        temp(k) = sum(R_K .* log(max(p_kk_tilde_one_minus(k, :), realmin)));
        r_K(k) = randg(gamma0/K + sum(L_KK(k, :)))./(c_0 - temp(k));
    end
    
    %% xi
    ell = sum(CRT_sum_mex_matrix(sparse(m_dot_k_k_dot_dot(diagdex))', xi * r_K'));
    xi = randg(1e-2 + ell)./(1e-2 - sum(r_K .* log(max(diag(p_kk_tilde_one_minus), realmin))));
    
    %% lambda
    R_KK = r_K * r_K';
    R_KK(diagdex) = xi * r_K;
    lambda_KK = zeros(K);
    lambda_KK(diagdex) = randg(m_dot_k_k_dot_dot(diagdex) + R_KK(diagdex))./(beta + theta_KK(diagdex));
    lambda_KK(triu1dex) = randg(m_dot_k_k_dot_dot(triu1dex) + R_KK(triu1dex))./(beta + theta_KK(triu1dex));
    lambda_KK = lambda_KK + triu(lambda_KK, 1)';
    
    %% beta 
    beta = randg(1 + sum(R_KK(diagdex)) + sum(R_KK(triu1dex)))./(1+ sum(lambda_KK(diagdex)) + sum(lambda_KK(triu1dex)) );
    
    %% gamma0
    L_k_tilde = CRT_sum_mex_matrix(sparse(sum(L_KK, 2)'), gamma0/K);
    p_ttilde_k_one_minus = c_0./(c_0 - temp);
    gamma0 = randg(e_0 + sum(L_k_tilde))./(h_0 - 1/K*sum(log(max(p_ttilde_k_one_minus, realmin))));
    
    %% c_0 
    c_0 = randg(1e-2 + gamma0)./(1e-2 + sum(r_K));
    
    %% phi
    
            %% P_n_k_t :: backward
            temp = sum(phi(:, :, T) * lambda_KK, 1);
            for i = 1:N
                temp_i = temp - phi(i, :, T) * lambda_KK;
                P_n_k_t(i, :, T) = temp_i ./ (c_n_t(i, t) + temp_i);
            end
            for t = (T-1):-1:1
                temp = sum(phi(:, :, t) * lambda_KK, 1);
                for i = 1:N
                    temp_i = temp - phi(i, :, t) * lambda_KK - log(max(1-P_n_k_t(i, :, t+1), realmin));
                    P_n_k_t(i, :, t) = temp_i./(c_n_t(i, t) + temp_i);
                end
            end

            % temp = sum(phi_0 * lambda_KK, 1);
%             for i = 1:N
%                 temp_i = temp - phi_0(i, :) * lambda_KK;
%                 P_n_k_0(i, :) = temp_i./(f_0 + temp_i);
%             end
            % P_n_k_0 = -log() 
            
            %% Y_n_k_t :: backward
            Y_n_k_t(:,:,T) = CRT_matrix(m_n_k_dot_dot(:, :, T)', phi(:, :, T-1));
            for t = (T-1):-1:2
                Y_n_k_t(:,:,t) = CRT_matrix(m_n_k_dot_dot(:, :, t)' + Y_n_k_t(:, :, t+1), phi(:, :, t-1));
            end
            Y_n_k_t(:, :, 1) = CRT_matrix(m_n_k_dot_dot(:, :, 1)' + Y_n_k_t(:, :, 2),phi_0);
            Y_n_k_0 = CRT_matrix(Y_n_k_t(:, :, 1), g_0(ones(N, K)));
            
            %% phi :: forward
            phi_0 = randg(g_0 + Y_n_k_t(:, :, 1))./(f_0 - log(max(1 - P_n_k_t(:, :, 1), realmin)));
            
            %%
            phi_lambda = phi(:, :, 1) * lambda_KK;
            temp = sum(phi_lambda, 1);
            for i = randperm(N)
                temp = temp - phi_lambda(i, :);
                phi(i, :, 1) = randg(phi_0(i, :) + Y_n_k_t(i, :, 2) + m_n_k_dot_dot(:, i, 1)')./(c_n_t(i, 1) + temp - log(max(1 - P_n_k_t(i, :, 2), realmin)));
                temp = temp + phi(i, :, 1) * lambda_KK;
            end
            %%
            
            for t = 2:(T-1)
                %%
                phi_lambda = phi(:, :, t) * lambda_KK;
                temp = sum(phi_lambda, 1);
                for i = randperm(N)
                    temp = temp - phi_lambda(i, :);
                    phi(i, :, t) = randg(phi(i, :, t-1) + Y_n_k_t(i, :, t+1) + m_n_k_dot_dot(:, i, t)')./(c_n_t(i, t) + temp - log(max(1 - P_n_k_t(i, :, t+1), realmin)));
                    
%                     if isnan(sum(phi(i, :, t))) || isinf(sum(phi(i, :, t)))
%                         xx = 1;
%                     end
                    
                    temp = temp + phi(i, :, t) * lambda_KK;
                    
%                     if isnan(temp)
%                         xx = 1;
%                     end
                end
                %%
            end
            
            %%
            phi_lambda = phi(:, :, T) * lambda_KK;
            temp = sum(phi_lambda, 1);
            for i = randperm(N)
                  temp = temp - phi_lambda(i, :);
                  phi(i, :, T) = randg(phi(i, :, T-1) + m_n_k_dot_dot(:, i, T)')./(c_n_t(i, T) + temp);
%                   if sum(phi(i, :, T))>1e3
%                       xx = 1;
%                   end
%                   if isnan(sum(phi(i, :, T))) || isinf(sum(phi(i, :, T)))
%                       xx = 1;
%                   end
                  temp = temp + phi(i, :, T)  * lambda_KK;
%                   if isnan(temp)
%                         xx = 1;
%                     end
            end
            
            %% 
            g_0 = randg(u_0 + sum(sum(Y_n_k_0)))./(v_0 - sum(sum(log(max(1 - P_n_k_0, realmin)))));
            
            %%
            f_0 = randg(g_0 * N * K + a_0)./(b_0 + sum(sum(phi_0)));
            
            %%
%             for i = 1:N
%                 c_n(i) = randg(1 + sum(phi_0(i, :)) + sum(sum(phi(i, 1:(T-1), :))))./(1 + sum(sum(phi(i, :, :))));
% %                 if c_n(i) < eps
% %                     xx = 1;
% %                 end
%             end

%             c_n_t(:, 1) = randg(1e-2 + sum(phi_0, 2))./(1e-2 + sum(phi(:, :, 1), 2));
%             for t = 2:T
%                 c_n_t(:, t) = randg(1e-2 + sum(phi(:, :, t-1), 2))./(1e-2 + sum(phi(:, :, t), 2));
%             end

%             c_n_t(:, 1) = gamrnd(1e-2 + sum(phi_0, 2), 1./(1e-2 + sum(phi(:, :, 1), 2)));
%             for t = 2:T
%                 c_n_t(:, t) = gamrnd(1e-2 + sum(phi(:, :, t-1), 2), 1./(1e-2 + sum(phi(:, :, t), 2)));
%             end

            
            for t = 1:T
                Prob = phi(:,:,t) * lambda_KK * phi(:,:,t)' + eps;
                Prob = 1 - exp(- Prob);
                if iter > burnin
                    Probsamps(:, :, t) = Probsamps(:, :, t) + Prob;
                    ProbAve{t} = Probsamps(:, :, t)/(iter - burnin);
                else
                    ProbAve{t} = Prob;
                end
                
                if iter > options.burnin && options.TrainRatio ~=1
                    BTrain{t}(idx_test{t}) = rand() < ProbAve{t}(idx_test{t});
                end
                
            end
            fprintf('iter: %03d, K:%02d,gamma0:%03f.\n', [iter Ktrace(iter) gamma0]);

            
            if mod(iter, 100) == 0 
                fprintf('iter: %03d, K:%02d,gamma0:%03f.\n', [iter Ktrace(iter) gamma0]);
                if options.display == 1
                    Tplot = 6;
                for t = 1:Tplot
                    subplot(Tplot, 4, (t-1)* 4 + 1);
                    imagesc(BTrain{t}+BTrain{t}');colormap('jet');
                    subplot(Tplot, 4, (t-1)* 4 + 2);
                    imagesc(ProbAve{t});colormap('jet');
                    subplot(Tplot, 4, (t-1)* 4 + 3);
                    imagesc(log(phi(:,:,t) * lambda_KK + 1e-2));colormap('jet');
                    subplot(Tplot, 4, (t-1)* 4 + 4);
                    bar(r_K./sum(r_K));colormap('jet');
                end
                drawnow;
                end
            end
end

result.timecost = toc
result.r_K = r_K;
result.K = sum(r_K./sum(r_K) > 0.05);
result.phi = phi;
result.lambda_KK = lambda_KK;
result.ProbAve = ProbAve;

%% quantify results
if isfield(options,'TrainRatio') && options.eval == 1
    if options.TrainRatio == 1
        idx_test = idx_train;        
    end
        rate_test = cell(T,1);
        links = cell(T, 1);
        coll_rate = [];
        coll_links = [];
        for t = 1:T
            B{t} = full(B{t});
            rate_test{t} = ProbAve{t}(idx_test{t});
            links{t} = B{t}(idx_test{t});
            coll_rate = [coll_rate;rate_test{t}];
            coll_links = [coll_links;links{t}];
        end
        % coll_rate = double(coll_rate > 1e-3);
        %figure(234);
        %subplot(1,2,1);
        %[XX, YY, TT, AUCroc] = perfcurve(coll_links, coll_rate, 1);
        %plot(XX,YY);
        %axis([0 1 0 1]), grid on, xlabel('FPR'), ylabel('TPR'), hold on;
        %x = [0:0.1:1];plot(x,x,'b--'), hold off; title(['AUCroc = ', num2str(AUCroc)])
        
        %subplot(1,2,2);
        %[prec, tpr, fpr, thresh] = prec_rec(coll_rate, coll_links,  'numThresh',3000);
        %plot([0; tpr], [1 ; prec]); % add pseudo point to complete curve
        %xlabel('recall');
        %ylabel('precision');
        %title('precision-recall graph');
        %AUCpr = trapz([0;tpr],[1;prec]);
        %F1= max(2*tpr.*prec./(tpr+prec));
        %title(['AUCpr = ', num2str(AUCpr), '; F1 = ', num2str(F1)])
        %fprintf('aucroc: %f, aucprec: %f.\n', [AUCroc, AUCpr]);
        %result.aucroc = AUCroc;
        %result.aucprec = AUCpr;
end
%%

end