%% MERGE2
%==========================================================================
%
%   Intermediate Merging of Samples Approaching to SEPs
%
%   
%   Implemented by Kyu-Hwan Jung
%   June 14, 2008
%
% * The source code is available under the GNU LESSER GENERAL PUBLIC
% LICENSE, version 2.1. 
%==========================================================================

function [State, arrival]=mergeBall(all_locals, arrival, R2)

N=size(all_locals,1);
Done=zeros(N,1);
State=all_locals;

for i=1:N
    if Done(i,1) == 0
        Done(i,1)=1;
        cur_sample=all_locals(i,:);
        inds=find(Done==0);
        dst=sqrt(dist2(cur_sample,all_locals(inds,:)));
        inds2=find(dst<R2);
        
        if length(inds2)>0 & arrival(i,1)==1
            inds3=find(arrival(inds(inds2))==0);
            if length(inds3)>0
                State((inds(inds2(inds3))),:)=repmat(cur_sample,length(inds3),1);
                arrival((inds(inds2(inds3))),1)=ones(length(inds3),1);
                Done((inds(inds2(inds3))),1)=ones(length(inds3),1);
            end
            
        elseif length(inds2)>0 & arrival(i,1)==0
            
            inds4=find(arrival(inds(inds2))==1);
            if length(inds4)==0
                State(inds(inds2),:)=repmat(cur_sample,length(inds2),1);
                Done(inds(inds2),1)=ones(length(inds2),1);
                
            else
                State(i,:)=all_locals(inds(inds2(inds4(1))),:);
                arrival(i,:)=1;
            end
        end
    end
    
end

end


