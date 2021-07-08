function [tinput, clstmodel] = smsc(input, support, supportopt, varargin) % output, support, varargin) % no options/ output??
% SMSC support function dynamics based labeling method using stable-equilibrium
%
% Description:
%  This method first decomposes a given data set into several disjoint groups,
%  where each group is represented by a stable equilibrium point to which all 
%  its members converge. 
%  To differentiate between groups that belong to different clusters, the 
%  complete graph (CG) labeling strategy, restricted to the set of the
%  stable equilibrium points {si} can be emploted using an adjacecy matrix
%  Aij between pairs of si and sj.Aij = 1, if R(y) < = R for all the points 
%  y on the line segment connecting si and sj and Aij = 0, otherwise. 
%  A pair of groups A(si) and A(sj)is then assigned to the same connected 
%  component if si and sj belong the the same connected components of the 
%  graph induced by A.
%
% * The source code is available under the GNU LESSER GENERAL PUBLIC
% LICENSE, version 2.1. 

%% Step 1 : Training SVDD
[tinput, supportmodel] = modelsupport(input, support, supportopt); 

%% Step 2 : Labeling Data for Clustering

clstmodel.support_model = supportmodel; 

tic;
% find stable equilibrium points
[rep_locals,locals,local_val,match_local]=findSEPs(tinput',supportmodel);
clstmodel.local=locals';    % dim x N_local

% calculates the adjacent matrix
[adjacent] = findAdjMatrix(clstmodel.local,supportmodel);

% Finds the cluster assignment of each data point
local_clusters_assignments = findConnectedComponents(adjacent);
clstmodel.local_clusters_assignments = local_clusters_assignments;
clstmodel.cluster_labels = local_clusters_assignments(match_local)';

labeling_time_S_MSC = toc