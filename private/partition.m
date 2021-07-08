function [C, clst] = partition(X,r)
% PARTITION Partitioning data points into small balls
% [C, clst] = partition(X,r)
% partitions data space X into small balls with radius r
% X [N x dim] : whole input data
% r : radius of balls
%
% C : center points of each balls
% clst [N x 1] : ball index of each data point
%
% * The source code is available under the GNU LESSER GENERAL PUBLIC
% LICENSE, version 2.1. 



[N d]=size(X); % N : number of data points    d : dimension of input
clst=zeros(N,1);

x_ind=[1:1:N];

X2=X;   % temp X
i=1;    % initial # of ball
while 1
    randind=randperm(size(X2,1));
    C(i,:)=X2(randind(1),:);    % random odering
    dst=sqrt(dist2(C(i,:),X2));
    ind=find(dst<r);
    if isempty(ind)
        error(message('R1 is too small! R1 must be positive real number'));
    end
    clst(x_ind(ind),1)=i;
    
    X2(ind,:)=[];
    x_ind(ind)=[];
    
    i=i+1;
    if isempty(X2)
        return;
    end
end

end