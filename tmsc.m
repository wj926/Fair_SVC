

function [tinput, clstmodel] = tmsc(input, support, supportopt, options, fairness,varargin ) 
% TMSC Transition Point Based Adjacency Decision Labeling Method
%
% Description:
%  In this method, an SEV, s_a, is said to be adjacent to another SEV, s_b,
%  if there exists an index-one saddle equilibrium vector d in the
%  intersection of boundaries of s_a's and s_b's basin. Such an index-one
%  saddle equilibrium vector, d, is called a transition equilibrium vector
%  (TEV) between s_a and s_b.
%  From a practical viewpoint, the notions of adjacent SEVs and TEVs enable
%  us to build a weighted graph Gr = (Vr,Er) describing the connections
%  between the SEVs with the following elements:  
%   1. The vertices Vr of Gr consist of SEVs, si, in V with f(s_i) < r.
%   2. The edge Er of Gr is defined as follows : (s_i,s_j) in Er with the
%   edge weight distance, dist(s_i, s_j) = f(d), if there is a TEV, d,
%   between s_i and s_j with dist(s_i, s_j)  < r.
%
%   Moreover, when we use the TEVs to build a weighted graph, we can
%   utilize dynamic dissimilarity measure to control the number of clusters
%   efficiently(hierarchicalLabelTSVC).
%
% * The source code is available under the GNU LESSER GENERAL PUBLIC
% LICENSE, version 2.1. 


%% Step 1 : Training SVDD
sensitive = input(:,3);
input = input(:,1:2);
[tinput, supportmodel] = modelsupport(input, support, supportopt);
%% Step 2 : Labeling Data for Clustering
clstmodel.options = options;
clstmodel.support_model = supportmodel; 

X = tinput; 
fHierarchical = options.hierarchical; 

[N, attr] = size(X);

% Find SEPs
[rep_locals,locals,local_val,match_local]=findSEPs(tinput',supportmodel);
nOfLocals=size(locals,1);
clstmodel.local=locals';

% Find transition points and label the SEPs
epsilon=options.epsilon;
[ts]=findTPs(locals,supportmodel,epsilon);   
nOfTS=length(ts.f);


%% Cluster assignment of each data point

% --- Automatic determination of cluster number based on the cluster boundary
if ~fHierarchical

    disp('Automatic determination of cluster numbers based on the SVDD boundearies defined by R^2');
    adjacent = zeros(nOfLocals,nOfLocals);
    tmp=find(ts.f<supportmodel.r+10^(-7));
    
    if ~isempty(tmp)    % only check the connectivity of TSs inside the sphere
        for j=1:length(tmp)
            adjacent(ts.neighbor(tmp(j),1),ts.neighbor(tmp(j),2))=1;
            adjacent(ts.neighbor(tmp(j),2),ts.neighbor(tmp(j),1))=1;
        end
        %% To connect nodes which can be connected via directly connected edges.
        for i=1:nOfLocals
            for j=1:i
                if (adjacent(i,j) == 1)
                    adjacent(i,:) = (adjacent(i,:) | adjacent(j,:));
                end
            end
            adjacent(i,i) = 1;
        end
        local_clusters_assignments = findConnectedComponents(adjacent);
    end  

    % model update
    clstmodel.ts.x=ts.x(tmp',:);
    clstmodel.ts.f=ts.f(tmp',:);
    clstmodel.ts.purturb=ts.purturb(tmp',:);
    clstmodel.ts.neighbor=ts.neighbor(tmp',:);
    clstmodel.ts.cuttingLevel=supportmodel.r;
    clstmodel.ts.wholeneighbor=ts.wholeneighbor(tmp',:);
    
    % cluster assignment into entire data points
    clstmodel.local_clusters_assignments = local_clusters_assignments;
    clstmodel.cluster_labels = local_clusters_assignments(match_local)';

    
    
else    
    K = options.K;

    ts.nOfLocals=nOfLocals;
    ts.sensitive=sensitive;
    switch fairness
        case 'worstupdate'
        [h_ts,h_cluster_labels,h_local_ass] = hierarchicalLabelTSVC_worstupdate(locals,match_local,ts,K,sensitive);
        case 'balanceincrease'
         [h_ts,h_cluster_labels,h_local_ass] = hierarchicalLabelTSVC_balanceincrease(locals,match_local,ts,K,sensitive);
        case 'nofair'
         [h_ts,h_cluster_labels,h_local_ass] = hierarchicalLabelTSVC_nofair(locals,match_local,ts,K,sensitive);
    end
    clstmodel.ts = h_ts;
    clstmodel.local_clusters_assignments = h_local_ass;
    clstmodel.cluster_labels = h_cluster_labels;

end



