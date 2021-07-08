function [out_ts,cluster_labels,local_ass] = hierarchicalLabelTSVC_balanceincrease(locals,match_local, ts, K ,sensitive)
% HIERARCHICALLABELTSVC 
%
% Description: 
%  In TMSC, we can tilize dynamic dissimilarity measure to control the
%  number of clusters efficiently 
%
% * The source code is available under the GNU LESSER GENERAL PUBLIC
% LICENSE, version 2.1. 
tmp = [];
    %SEP별 정확한 개수랑 balance계산하기 
    nOfLocals=size(locals,1);
    originalNoffirst = zeros(nOfLocals,1);
    originalNoftotal = zeros(nOfLocals,1);
    balance = zeros(nOfLocals,1);
    for i=1:nOfLocals
        originalNoffirst(i) = sum(sensitive(match_local==i));
        originalNoftotal(i) = length(sensitive(match_local==i));
        originalbalance(i) = min([originalNoffirst(i)/originalNoftotal(i),1-originalNoffirst(i)/originalNoftotal(i)]);
    end
    
    
    %초기 balance의 증가율
    nOfTS=length(ts.f);
    ts.balance = zeros(nOfTS,1);
    changednewbalance = zeros(nOfTS,1);
    ts.neighbor = ts.neighbor;
    
    for i=1:nOfTS
        idx1 =  ts.neighbor(i,1);
        idx2 =  ts.neighbor(i,2);       
        changednewbalance(i) = min([(originalNoffirst(idx1)+originalNoffirst(idx2))/(originalNoftotal(idx1)+originalNoftotal(idx2)),1-(originalNoffirst(idx1)+originalNoffirst(idx2))/(originalNoftotal(idx1)+originalNoftotal(idx2))]);
        %제일 개선되는 방향
        ts.balance(i) = max([changednewbalance(i)-originalbalance(idx1),changednewbalance(i)-originalbalance(idx2)]);

        %%%제일 안좋은애로 진행
        %ts.balance(i) = max([0.5-originalbalance(idx1),0.5-originalbalance(idx2)]);
    end
    
    nOfLocals=size(locals,1);
    nOfTS=length(ts.f);
   
    local_clusters_assignments=[];
    
    %%%%%%%%%%%%%%%%%%초기기준은 여기서 바꿈!!
    %%이것이 기본 svc(지금 잘 안됨)
    %%ts.sorting = ts.f;
    
    %이거는 우리의 방법 적용
    ts.sorting = -1000*ts.balance + ts.f;
    
    %ts.f = ts.f;
    f_sort=sort(ts.sorting);
        
    adjacent = zeros(nOfLocals,nOfLocals,nOfTS);
    a=[];
    flag=0;
    
    %초기버전의 iter마다 계속 바뀌는 first,total,balance 할당한다. 
    updatefirst = originalNoffirst;
    updatetotal = originalNoftotal;
    updatebalance = originalbalance;
    
    for m=1:nOfTS
        %돌때마다 balance 증가율 계산한다. 
        nOfTS=length(ts.f);
        ts.balance = zeros(nOfTS,1);
        changednewbalance = zeros(nOfTS,1);
        ts.neighbor = ts.neighbor;
        
        for i=1:nOfTS
            idx1 =  ts.neighbor(i,1);
            idx2 =  ts.neighbor(i,2);
            %%만약 합치려는애가 이미 합쳐진놈이면 penalty엄청 줘서 안합쳐지게 진행
            if (updatefirst(idx1)==updatefirst(idx2) && updatetotal(idx1)==updatetotal(idx2))
                ts.balance(i)=-100000;
                
            else
                %제일 개선되는 방향
                %ts.balance(i) = max([changednewbalance(i)-updatebalance(idx1),changednewbalance(i)-updatebalance(idx2)]);
                changednewbalance(i) = min([(updatefirst(idx1)+updatefirst(idx2))/(updatetotal(idx1)+updatetotal(idx2)),1-(updatefirst(idx1)+updatefirst(idx2))/(updatetotal(idx1)+updatetotal(idx2))]);
                %ts.balance(i) = max([changednewbalance(i)-updatebalance(idx1),changednewbalance(i)-updatebalance(idx2)]);
                ts.balance(i) = (changednewbalance(i)-updatebalance(idx1))+(changednewbalance(i)-updatebalance(idx2));

                %%%제일 안좋은애로 진행
                %ts.balance(i) = max([0.5-updatebalance(idx1),0.5-updatebalance(idx2)]);
            
            end
        end
        
        %%%%%%%여기도 기준 바꿔야함 
        %기존 svc
        %ts.sorting  = ts.f;
        
        %우리의 방법으로 sorting
        ts.sorting  = -1000*ts.balance + ts.f;
        %ts.f = ts.f;
        [f_sort,idx]=sort(ts.sorting);
    
        %제일 작은 sorting의 neighbor의 index를 구해야함. 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %neighborone = ts.neighbor(idx(m),1)
        %neighbortwo = ts.neighbor(idx(m),2)

        neighborone = ts.neighbor(idx(1),1);
        neighbortwo = ts.neighbor(idx(1),2);
        
        %둘이 연결 시킨다.  
        adjacent(neighborone,neighbortwo,:)=1;
        adjacent(neighbortwo,neighborone,:)=1;
        %% To connect nodes which can be connected via directly connected edges.
        for i=1:nOfLocals
            for j=1:i
                if (adjacent(i,j,:) == 1)
                    adjacent(i,:,:) = (adjacent(i,:,:) | adjacent(j,:,:));
                end
            end
            adjacent(i,i) = 1;                
        end

        tmp = [tmp idx(1)];
        my_ts.x=ts.x(tmp',:);
        my_ts.f=ts.f(tmp',:);
        my_ts.purturb=ts.purturb(tmp',:);
        my_ts.neighbor=ts.neighbor(tmp',:);
        my_ts.wholeneighbor=ts.neighbor;
        tmp_ts{m}=my_ts;
        assignment = findConnectedComponents(adjacent(:,:,m));

        %새롭게 fist, second 업데이트 시킨다. 
        
        for i=1:max(assignment)
            assigncluster = find(assignment==i);
            tmpfirst = sum(originalNoffirst(assigncluster));
            tmptotal = sum(originalNoftotal(assigncluster));
            for j=assigncluster
                updatefirst(j)=tmpfirst;
                updatetotal(j)=tmptotal;
            end
        end
        
        
        
        %새롭게 balance 업데이트
        for i=1:nOfLocals
            updatebalance(i) = min([updatefirst(i)/updatetotal(i),1-updatefirst(i)/updatetotal(i)]);
        end
        
        if max(assignment)==K
            %disp('We can find the number of K clusters');         
            % clstmodel update
            out_ts=tmp_ts{m};
            %out_ts=[0];
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
        %out_ts=tmp_ts{ind(1)};
        out_ts=[0];
        
        local_clusters_assignments=local_clusters_assignments(:,ind(1));
        local_ass = local_clusters_assignments;
        cluster_labels = local_ass(match_local);
    end
end

