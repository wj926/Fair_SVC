function final_label=matchBallIndex(clst,label,C,EVs)

[n1 d1]=size(C);
[n2 d2]=size(EVs);
final_label=zeros(size(clst,1),1);

for i=1:n1
    ind=find(clst==i);
    final_label(ind,:)=label(i,1);
end
    