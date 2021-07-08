function [tinput, clstmodel] = fmsc(input, support, supportopt, options, varargin) 
% FMSC Fast Support-Based Cluster Labeling Method
%
% Description:
%  To expedite the clustering process, in the initial stage, this method
%  decompose the whole data sample into separate small balls. In the main
%  stage, we iteratively merge small balls at each iteration.
%  There are two types of intermediate merging procedures. When there exist
%  a SEV that is closer than pre-specified merging parameter RHO, then
%  that centers are regarded as already converged to that SEV. Otherwise,
%  check if there exist other centers that are within the range of merging
%  parameter and then merge these centers into one if this is the case.
%  After obtaining all SEVs from centers, we can assign cluster label to
%  each SEVs by investigating the adjacency matrix.          
%
% * The source code is available under the GNU LESSER GENERAL PUBLIC
% LICENSE, version 2.1. 

%% Step 1 : Training SVDD
[tinput, supportmodel] = modelsupport(input, support, supportopt);

%% Step 2 : Labeling Data for Clustering
clstmodel.options = options;
clstmodel.support_model = supportmodel; 


% Partitioning Data into Small Ball
R1 = options.R1;
R2 = options.R2;
[centers, clst] = partition(tinput',R1);
% NC = size(Centers,1);

% calculates the adjacent matrix
[all_locals, unique_locals, match_locals] = labelFSVC(centers,supportmodel,R2); 

clstmodel.local = unique_locals';
[adjacent] = findAdjMatrix(unique_locals',supportmodel);
local_clusters_assignments = findConnectedComponents(adjacent);
cluster_labels = local_clusters_assignments(match_locals);

% Finds the cluster assignment of each data point
clstmodel.local_clusters_assignments = local_clusters_assignments;
clstmodel.cluster_labels = matchBallIndex(clst, cluster_labels, centers, unique_locals);
clstmodel.centers = centers;
clstmodel.match_locals = match_locals;

end






