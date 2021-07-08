%==========================================================================
%
%   Function to find local minimum of trained kernel radius function
%
%   Return Values:
%       N_locals: local min corresponding to each sample
%      local: unique local mins
%
%
%==========================================================================
% January 13, 2009
% Implemented by Daewon Lee
%
% * The source code is available under the GNU LESSER GENERAL PUBLIC
% LICENSE, version 2.1. 
%==========================================================================

function [N_locals,local,local_val,match_local]=findSEPs(X,model)

[N dim]=size(X);

N_locals=[];
local_val=[];

switch model.support_type
    
    case 'GP' % for GP support function
        
        for i=1:N
            x0=X(i,:);
            options = optimset('Display','off','GradObj','on');
            [temp val]=fminsearch(@my_R_GP,x0,[],model);
            N_locals=[N_locals; temp];
            local_val=[local_val;val];
        end
        
    case 'SVDD'
        
        for i=1:N
            x0=X(i,:);
            options = optimset('Display','off','LargeScale','on','GradObj','on');
            if length(x0)<=2
                [temp val]=fminsearch(@my_R,x0,[],model);   % Nelder-Mead
                %             [temp val]=newfminsearch(@my_R,x0);   % Nelder-Mead
            else
                [temp val]=fminunc(@my_R,x0,options,model); % trust region method
            end
            N_locals=[N_locals; temp];
            local_val=[local_val;val];
        end
end

[local,I,match_local]=unique(round((10*N_locals)),'rows');
local=N_locals(I,:);

% for i=1:length(I)
%     tmp=find(match_local==i);
%     if length(tmp)<=1
%         one_ind=I(match_local(tmp(1)));
%         [dummy,ind]=sort(dist2(local,N_locals(one_ind,:)));
%         N_locals(one_ind,:)=repmat(N_locals(ind(2),:),length(tmp),1);
%     end
%
% end
% [local,I,match_local]=unique(round((10*N_locals)),'rows');
% local=N_locals(I,:);
%local=N_locals(I,:);

local_val=local_val(I,:);


