function mscpath(toolboxroot)
if nargin < 1
   toolboxroot = pwd ; % get toolbox directory
end

% disp('Adding path for the Statistical Pattern Recognition Toolbox...');
% path for UNIX
p = [
    '$kernels:',...
    '$optimization:',...
    '$visual:'
%      '$kernels:',...
%      '$kernels/extraction:',...
%      '$kernels/preimage:'
    ];

p=translate(p,toolboxroot);

% adds path at the start
addpath(p);


function p = translate(p,toolboxroot)
%TRANSLATE Translate unix path to platform specific path
%   TRANSLATE fixes up the path so that it's valid on non-UNIX platforms
%
% This function was derived from MathWork M-file "pathdef.m"

cname = computer;
% Look for VMS, this covers VAX_VMSxx as well as AXP_VMSxx.
%if (length (cname) >= 7) & strcmp(cname(4:7),'_VMS')
%  p = strrep(p,'/','.');
%  p = strrep(p,':','],');
%  p = strrep(p,'$toolbox.','toolbox:[');
%  p = strrep(p,'$','matlab:[');
%  p = [p ']']; % Append a final ']'

% Look for PC
if strncmp(cname,'PC',2)
  p = strrep(p,'/','\');
  p = strrep(p,':',';');
  p = strrep(p,'$',[toolboxroot '\']);

% Look for MAC
%elseif strncmp(cname,'MAC',3)
%  p = strrep(p,':',':;');
%  p = strrep(p,'/',':');
%  m = toolboxroot;
%  if m(end) ~= ':'
%    p = strrep(p,'$',[toolboxroot ':']);
%  else
%    p = strrep(p,'$',toolboxroot);
%  end
else
  p = strrep(p,'$',[toolboxroot '/']);
end
