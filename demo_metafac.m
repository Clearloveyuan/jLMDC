% DEMO_METAFAC - demonstrate the metafac algorithm
% See also METAFAC

% Author: Yu-Ru Lin <yu-ru.lin@asu.edu>, August 2008

fprintf('\nDemo MetaFac algorithm...');
clear RG; clear Vdims;
K = 40; % no. of factors (or clusters)
prior={}; % test without prior
w=1; % test with equal weights for all relations

fprintf('\n => generate 3-order data tensor for R1');
D = sptenrand([1000 1000 1000],1e+4);
RG{1} = {D,1:3,w};
Vdims = [size(D)];
Vdims,RG
fprintf('\n => run MetaFac for R1');
[S,U,iters,cvtime]=metafac(RG, Vdims, K, prior);

fprintf('\n => generate 2-order data tensor for R2 (which shares 1st entity with R1)');
D = sptenrand([1000 1000 ],1e+3);
RG{2} = {D,[1 length(Vdims)+1],w};
Vdims = [Vdims size(D,2)];
Vdims,RG
fprintf('\n => run MetaFac for R12');
[S,U,iters,cvtime]=metafac(RG, Vdims, K, prior);

fprintf('\n => generate 2-order data tensor for R3 (which shares 2nd entity with R1)');
D = sptenrand([1000 1000 ],1e+3);
RG{3} = {D,[2 length(Vdims)+1],w};
Vdims = [Vdims size(D,2)];
Vdims,RG
fprintf('\n => run MetaFac for R123');
[S,U,iters,cvtime]=metafac(RG, Vdims, K, prior);

