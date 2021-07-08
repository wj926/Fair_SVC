function [tinput, clstmodel] = cgsc(input, supporttype, supportopt, varargin)
% CGSVC Support Vector Clusteing using Complete-Graph Based Labeling Method
%
% Description:
%  To determine whether a pair of xi and xj is in the same contour, 
%  it can be used a complete-graph(CG) strategy that relies on the fact 
%  that any path connecting two data points in different contours must 
%  exit the contours in data space, which is equivalent that the image 
%  of the path in feature space exits the minimum enclosing sphere.
%  CG strategy makes adjacency matrix Aij between pairs of points 
%  xi and xj as follows :
%   A(ij) = 1 if for all y on the line segment connecting xi and xj
%           0 otherwise
%
% * The source code is available under the GNU LESSER GENERAL PUBLIC
% LICENSE, version 2.1. 


%% Step 1 : Training SVDD
[tinput, supportmodel] = modelsupport(input, supporttype, supportopt); % output, support);
clstmodel.support_model = supportmodel; 

%% Step 2 : Labeling Data for Clustering
tic;
adjacent = findAdjMatrix(tinput,supportmodel);
% Finds the cluster assignment of each data point
clusters = findConnectedComponents(adjacent);
clstmodel.cluster_labels=double(clusters);

labeling_time_CG_SC = toc