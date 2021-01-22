function [S,U,iters,cvtime] = metafac(RG,Vdims,rank,prior)
% METAFAC - relatoinal non-negative tensor decomposition
%   This function is a wrapper for the MetaFac implementation,
%   metafac_imp.p
%
% Output:
%   S - 1xK vector; the core tensor (PARAFAC-like model, super diagonal)
%   U - cell array of factors {Ui}; each Ui is a NixK matrix and each
%   column of Ui sum to 1
%
% Input:
%   RG - relational graph (metagraph)
%     RG = {{Xj,[1,2,4],wj}, {...}, ...} is a cell array where each cell
%     has 3 components for a relation:
%       Xj - the data tensor for relation Rj
%       [1,2,...,m]: the indices of incident facets (corresponding to the
%       indices in Vdims)
%       wj - weight of Rj
%   Vdims - the dimensions of all facets
%     Vdims = [N1,N2,...,Ni] is a vector containing the dimensions of
%     the first, second, ..., the i-th facet
%   rank - number of factors to be decomposed
%   prior - prior distribution of S and U; 
%     prior = {pS,pU,pa} is a cell arrary containing three components
%       pS - the prior core tensor (same format as the output S)
%       pU - prior factors in a cell array {pUi} (same format as the output U)
%       pa - the weight for prior
%
% Configuration:
%   The algorithm can be configured by editing metafac_config.txt
%
% Example:
%   See demo_metafac.m
%
% See also DEMO_METAFAC

% Author: Yu-Ru Lin <yu-ru.lin@asu.edu>, August 2008
% Reference:
%   Yu-Ru Lin et. al. , "MetaFac: Commmunity Discovery via Relational
%   Hypergraph Factorization", SIGKDD 2009

[S,U,iters,cvtime] = metafac_imp(RG,Vdims,rank,prior);
