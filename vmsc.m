function [tinput, clstmodel] = vmsc(input, support, supportopt, voronoiopt, varargin)
% VMSC Voronoi Cell-based Kernel Support Clustering
%
% Description:
%  To achieve fast cluster labeling and sparse estimation of the support,
%  this method employs two main approximations for fast processing. First,
%  the support estimate for a large number of data points is approximated
%  by the support estimate from its small sampled data set. Second, the
%  partitioned cluster boundaries that use the basin cells are approximated
%  by the Voronoi cells. The Voronoi cells allow the assignment of a
%  cluster label of each data to that of the nearest REP without an
%  explicit call for the kernel support function. 
%  (Slow Phase) To label the remaining data more accurately, if the
%  ambiguity ratio is greater than of equal to the fixed constant, use
%  dynamic system and find the REP from the point. 
% 
% * The source code is available under the GNU LESSER GENERAL PUBLIC
% LICENSE, version 2.1. 

%% SAMPLING

sampleRate = voronoiopt.samplerate;
[dataNum,dim] = size(input);
sampleNum = dataNum * sampleRate; % sample number for each class

[sampleData,sampleIndex] = datasample(input,sampleNum,'Replace',false);
testIndex = setdiff(1:dataNum,sampleIndex);
testData = input(testIndex,:);

%% Construct a support function using Samples and Find REPs for Samples

% Selecting Labeling Method
method = voronoiopt.voronoiMethod; 

switch method
    case 'T-MSC'
        [tsam, clst_model] = tmsc(sampleData, support, supportopt, voronoiopt);
    case 'F-MSC'
        [tsam, clst_model] = fmsc(sampleData, support, supportopt, voronoiopt);
end

sepList=clst_model.local';
sepLabel=clst_model.local_clusters_assignments;

% transfrom input
switch support
    case 'SVDD'
        testData = my_normalize(testData');
    case 'GP'
        testData = my_normalize2(testData);
end

C = knnclassify(testData',sepList,sepLabel,1);


%% Slow Phase

eta = voronoiopt.eta; % 0.5 <= eta <= 1
[idx,dist] = knnsearch(sepList, testData','k',2);
PP = sepLabel(idx);
clstind = (PP(:,1)~=PP(:,2));
distind = ((dist(:,1)./dist(:,2)) > eta);
wrongind = (clstind & distind);

tempLocal = findSEPs(testData(:,wrongind)',clst_model.support_model);
CNew = knnclassify(tempLocal, sepList, sepLabel,1);
C(wrongind) = CNew;

labels = zeros(dataNum,1);
labels(sampleIndex,1) = clst_model.cluster_labels;
labels(testIndex,1) = C;

tinput = zeros(dim,dataNum);
tinput(:,sampleIndex) = tsam;
tinput(:,testIndex) = testData;

clstmodel.options = voronoiopt;
clstmodel.support_model = clst_model.support_model;
clstmodel.local = clst_model.local;
clstmodel.local_clusters_assignments = clst_model.local_clusters_assignments;
clstmodel.cluster_labels = labels;



