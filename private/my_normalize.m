%==========================================================================
% Normalizing input for SVDD
%
% Implemented by Kyu-Hwan Jung, April 26, 2010
%
% * The source code is available under the GNU LESSER GENERAL PUBLIC
% LICENSE, version 2.1. 
%==========================================================================

function X_normal = my_normalize(Xin)
% input : Xin [dim x num_data]
% output : X_normal [dim x num_data]

[dim, n]=size(Xin);

means=repmat(mean(Xin,2),1,n);
stds = repmat(std(Xin,0,2),1,n);

X_normal=(Xin-means)./stds;


