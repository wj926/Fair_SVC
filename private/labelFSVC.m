%% FSVC2
%==========================================================================
%
%   Fast Labeling using Support Function
%
%   Return Values:
%      N_locals: local min corresponding to each sample
%      local: unique local mins
%
%
%   April 26, 2010
%   Implemented by Kyu-Hwan Jung
%
% * The source code is available under the GNU LESSER GENERAL PUBLIC
% LICENSE, version 2.1. 
%==========================================================================

function [all_locals,unique_locals,match_locals] = labelFSVC(centers,model,R2)

[N dim]=size(centers);

N_locals=[];
local_val=[];
arrival=zeros(N,1);

all_locals=centers;
converge=0;
iter=1;
while(converge==0)
    ind=find(arrival==0);
    run_Samples=all_locals(ind,:);
    
    switch model.support_type
        
        case 'SVDD'
            
            for i=1:length(ind)
                x0=run_Samples(i,:);
                options = optimset('Display','off','LargeScale','on','GradObj','on','Hessian','on','MaxIter',5);
                if length(x0)<=2
                    [temp val]=fminsearch(@my_R,x0,[],model);
                else
                    [temp val]=fminunc(@my_R,x0,options,model);
                end
                all_locals(ind(i),:)=temp;
                local_val(ind(i),:)=val;
                
                if sum((temp-x0).^2)<10^(-3)
                    arrival(ind(i),1)=1;                % convergence check
                end
            end
            
        case 'GP'
            
            for i=1:length(ind)
                x0=run_Samples(i,:);
                options = optimset('Display','off','LargeScale','on','GradObj','on','Hessian','on','MaxIter',5);
                if length(x0)<=2
                    [temp val]=fminsearch(@my_R_GP,x0,[],model);
                else
                    [temp val]=fminunc(@my_R_GP,x0,options,model);
                end
                all_locals(ind(i),:)=temp;
                local_val(ind(i),:)=val;
                
                if sum((temp-x0).^2)<10^(-3)
                    arrival(ind(i),1)=1;                % convergence check
                end
            end       
            
    end
    
    
    if sum(arrival) > 0
        [all_locals, arrival]=mergeBall(all_locals, arrival, R2);
    end
    
    if sum(arrival)==N
        converge=1;    % Convergence Check
    end
    
    iter=iter+1;
end


%all_locals=round((10*all_locals)); % round off errors
[unique_locals,I,match_locals]=unique(round((10*all_locals)),'rows');
unique_locals=all_locals(I,:);

%unique_locals=all_locals(I,:);
local_val=local_val(I,:);

end