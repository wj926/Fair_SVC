function [balance_increase, worstbalance] = calculate_balance(assignmentofclusters, sensitive, options, varargin) 
    %���� �뷱�� 
    nOfClusters = max(nOfClusters);
    Noffirst = zeros(nOfClusters,1);
    Noftotal = zeros(nOfClusters,1);
    balance = zeros(nOfClusters,1);
    for i=1:nOfClusters
        Noffirst(i) = sum(sensitive(assignmentofclusters==i));
        Noftotal(i) = length(sensitive(assignmentofclusters==i));
        balance(i) = min([Noffirst(i)/Noftotal(i),1-Noffirst(i)/Noftotal(i)]);
    end
    
    %balance�� ������
    ts.balance = zeros(nOfTS,1);
    newbalance = zeros(nOfTS,1);
    ts.neighbor = ts.neighbor;
    for i=1:nOfTS
        idx1 =  ts.neighbor(i,1);
        idx2 =  ts.neighbor(i,2);       
        newbalance(i) = min([(Noffirst(idx1)+Noffirst(idx2))/(Noftotal(idx1)+Noftotal(idx2)),1-(Noffirst(idx1)+Noffirst(idx2))/(Noftotal(idx1)+Noftotal(idx2))]);
        %���� �����Ǵ� ����
        %ts.balance(i) = max([newbalance(i)-balance(idx1),newbalance(i)-balance(idx2)]);
        
        %%%���� �������ַ� ����
        ts.balance(i) = max([0.5-balance(idx1),0.5-balance(idx2)]);
    end
