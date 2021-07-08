function [tinput, clstmodel] = tmsc(input, support, supportopt, options, varargin) 
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
    clstmodel.ts.x=ts.x(tmp,:);
    clstmodel.ts.f=ts.f(tmp,:);
    clstmodel.ts.purturb=ts.purturb(tmp,:);
    clstmodel.ts.neighbor=ts.neighbor(tmp,:);
    clstmodel.ts.cuttingLevel=supportmodel.r;
    
    % cluster assignment into entire data points
    clstmodel.local_clusters_assignments = local_clusters_assignments;
    clstmodel.cluster_labels = local_clusters_assignments(match_local)';

    
    
else    
    K = options.K;

    %Calculate Balance
    %기존 밸런스 
    Noffirst = zeros(nOfLocals,1);
    Noftotal = zeros(nOfLocals,1);
    balance = zeros(nOfLocals,1);
    for i=1:nOfLocals
        Noffirst(i) = sum(sensitive(match_local==i));
        Noftotal(i) = length(sensitive(match_local==i));
        balance(i) = min([Noffirst(i)/Noftotal(i),1-Noffirst(i)/Noftotal(i)]);
    end
    
    %balance의 증가율
    ts.balance = zeros(nOfTS,1);
    newbalance = zeros(nOfTS,1);
    ts.neighbor = ts.neighbor;
    for i=1:nOfTS
        idx1 =  ts.neighbor(i,1);
        idx2 =  ts.neighbor(i,2);       
        newbalance(i) = min([(Noffirst(idx1)+Noffirst(idx2))/(Noftotal(idx1)+Noftotal(idx2)),1-(Noffirst(idx1)+Noffirst(idx2))/(Noftotal(idx1)+Noftotal(idx2))]);
        %제일 개선되는 방향
        %ts.balance(i) = max([newbalance(i)-balance(idx1),newbalance(i)-balance(idx2)]);
        
        %%%제일 안좋은애로 진행
        ts.balance(i) = max([0.5-balance(idx1),0.5-balance(idx2)]);
    end
    ts.balance = -1*ts.balance;
    bbbb= ts.balance
    [h_ts,h_cluster_labels,h_local_ass] = hierarchicalLabelTSVC(locals,match_local,ts,K,sensitive);
    
    clstmodel.ts = h_ts;
    clstmodel.local_clusters_assignments = h_local_ass;
    clstmodel.cluster_labels = h_cluster_labels;

end
