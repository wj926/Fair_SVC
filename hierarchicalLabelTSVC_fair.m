function [out_ts,cluster_labels,local_ass] = hierarchicalLabelTSVC_fair(locals,match_local, ts, K ,sensitive)
% HIERARCHICALLABELTSVC 
%
% Description: 
%  In TMSC, we can tilize dynamic dissimilarity measure to control the
%  number of clusters efficiently 
%
% * The source code is available under the GNU LESSER GENERAL PUBLIC
% LICENSE, version 2.1. 

    nOfLocals=size(locals,1);
    nOfTS=length(ts.f);
   
    local_clusters_assignments=[];
    %일단 이렇게 바꿈 
    ts.f = ts.f;
    %ts.f = ts.f;
    f_sort=sort(ts.f);
    
    adjacent = zeros(nOfLocals,nOfLocals,nOfTS);
    a=[];
    flag=0;
    for m=1:nOfTS
        cur_f=f_sort(end+1-m);    % cutting level:large --> small  (small number of clusters --> large number of clusters)
        %cur_f=f_sort(i);         % cutting level: small --> large (large number of clusters --> small number of clusters)
        
        tmp=find(ts.f<cur_f);
        if ~isempty(tmp) % TSs inside the sphere
            for j=1:length(tmp)
                adjacent(ts.neighbor(tmp(j),1),ts.neighbor(tmp(j),2),m)=1;
                adjacent(ts.neighbor(tmp(j),2),ts.neighbor(tmp(j),1),m)=1;
            end
            %% To connect nodes which can be connected via directly connected edges.
            for i=1:nOfLocals
                for j=1:i
                    if (adjacent(i,j,m) == 1)
                        adjacent(i,:,m) = (adjacent(i,:,m) | adjacent(j,:,m));
                    end
                end
                adjacent(i,i) = 1;                
            end
         
        end % end of current TS

        %%아래다가 넣어야함 원래
        a=[a;cur_f];
        my_ts.x=ts.x(tmp,:);
        my_ts.f=ts.f(tmp,:);
        my_ts.purturb=ts.purturb(tmp,:);
        my_ts.neighbor=ts.neighbor(tmp,:);
        my_ts.cuttingLevel=cur_f;
        ind=find(ts.f==cur_f);
        my_ts.levelx=ts.x(ind(1),:);
        tmp_ts{m}=my_ts;

        assignment = findConnectedComponents(adjacent(:,:,m));

        
        if max(assignment)==K
            disp('We can find the number of K clusters');         
            % clstmodel update
            out_ts=tmp_ts{m};    
            % cluster assignment into entire data points
            local_ass = assignment;
            cluster_labels = local_ass(match_local)';       
            flag=1;
            break;
        end

        local_clusters_assignments = [local_clusters_assignments assignment];                        
    end % end of K-control
    
    % cannot find k clusters     
    if flag==0
        disp('Cannot find cluster assignments with K number of clusters, instead that we find cluster assignments the with the nearest number of clusters to K !');
        [dummy,ind]=min(dist2(max(local_clusters_assignments)',K));
        
        %ts=[];
        out_ts=tmp_ts{ind(1)};
        local_clusters_assignments=local_clusters_assignments(:,ind(1));
        local_ass = local_clusters_assignments;
        cluster_labels = local_ass(match_local);
    end
end

