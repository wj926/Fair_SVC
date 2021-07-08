%==========================================================================
%
%   Plotting 2D Clustering Results
%
%
%   Implemented by Daewon Lee
%   and Modified by Kyu-Hwan Jung for GP, April 26, 2010
%   and Modified by Sujee Lee
%
%
% * The source code is available under the GNU LESSER GENERAL PUBLIC
% LICENSE, version 2.1. 
%==========================================================================

function varargout=plotmsc(input,clstmodel)
supportmodel = clstmodel.support_model;
minx=min(input(1,:));
maxx=max(input(1,:));
miny=min(input(2,:));
maxy=max(input(2,:));

xrange=maxx-minx;
yrange=maxy-miny;

[Ax,Ay] = meshgrid(linspace(minx-0.1*xrange,maxx+0.1*xrange,100),linspace(miny-0.1*yrange,maxy+0.1*yrange,100));

dist = kradius([Ax(:)';Ay(:)'],supportmodel);

figure;
ppatterns(input);
title('Original Data')

if isequal(supportmodel.support_type,'SVDD')
    
    if isequal(clstmodel.labeling_method,'S-MSC') || isequal(clstmodel.labeling_method,'T-MSC') || isequal(clstmodel.labeling_method,'F-MSC')
        figure; hold on;
        test.X=input;
        test.y=clstmodel.cluster_labels;
        ppatterns2(input(:,supportmodel.sv_ind),'ro',10);   % SV
        plot(clstmodel.local(1,:),clstmodel.local(2,:),'rs','MarkerSize',10); % SEP

        %%%%%%일단 잠시 아래 지움
        if isfield(clstmodel,{'ts'})
            plot(clstmodel.ts.x(:,1),clstmodel.ts.x(:,2),'r+','MarkerSize',10);
            legend('SVs','SEPs','TSs','Location','SouthEast');
        else
            legend('SVs','SEPs','Location','SouthEast');
        end
        ppatterns2(test);

%         if isequal(clstmodel.labeling_method,'T-SVC') %0902
%             contour( Ax, Ay, reshape(dist,100,100),[clstmodel.ts.cuttingLevel clstmodel.ts.cuttingLevel],'k'); 
%         else
%             contour( Ax, Ay, reshape(dist,100,100),[supportmodel.r supportmodel.r],'k'); 
%         end
%        contour( Ax, Ay, reshape(dist,100,100),[supportmodel.r supportmodel.r],'k'); 
        title(strcat(['Support Type: ', supportmodel.support_type,', Labelling method: ', clstmodel.labeling_method]));
        
        for i=1:length(clstmodel.ts.wholeneighbor)
            firstpoint = clstmodel.local(:,clstmodel.ts.wholeneighbor(i,1));
            secondpoint = clstmodel.local(:,clstmodel.ts.wholeneighbor(i,2));
            line([firstpoint(1) secondpoint(1)],[firstpoint(2) secondpoint(2)],'Color','blue');
        end       
        
        for i=1:length(clstmodel.ts.neighbor)
            firstpoint = clstmodel.local(:,clstmodel.ts.neighbor(i,1));
            secondpoint = clstmodel.local(:,clstmodel.ts.neighbor(i,2));
            line([firstpoint(1) secondpoint(1)],[firstpoint(2) secondpoint(2)],'Color','red');
        end    

        
        hold off
    elseif isequal(clstmodel.labeling_method,'V-MSC')
                figure; hold on;
        test.X=input;
        test.y=clstmodel.cluster_labels;
        ppatterns2(test);
        title(strcat(['Support Type: ', supportmodel.support_type,', Labelling method: ', clstmodel.labeling_method]));
        hold off
    else
        figure; hold on;
        test.X=input;
        test.y=clstmodel.cluster_labels;
        ppatterns2(input(:,supportmodel.sv_ind),'ro',10);   % SV
        ppatterns2(input(:,supportmodel.bsv_ind),'k.',5);  % BSV
        ppatterns2(test);
        contour( Ax, Ay, reshape(dist,100,100),[supportmodel.r supportmodel.r],'k');
        legend('SVs','BSV (outliers)','Location','SouthEast');
        title(strcat(['Support Type: ', supportmodel.support_type,', Labelling method: ', clstmodel.labeling_method]));
        hold off
    end
    
else
    if isequal(clstmodel.labeling_method,'S-MSC') || isequal(clstmodel.labeling_method,'T-MSC')  || isequal(clstmodel.labeling_method,'F-MSC')
        figure; hold on;
        test.X=input;
        test.y=clstmodel.cluster_labels;
        plot(clstmodel.local(1,:),clstmodel.local(2,:),'rs','MarkerSize',10); % SEP
        if isfield(clstmodel,{'ts'})
            plot(clstmodel.ts.x(:,1),clstmodel.ts.x(:,2),'r+','MarkerSize',10);
            legend('SEPs','TSs','Location','SouthEast');
        else
            legend('SEPs','Location','SouthEast');
        end
        ppatterns2(test);
%         if isequal(clstmodel.labeling_method,'T-SVC') %0902
%             contour( Ax, Ay, reshape(dist,100,100),[clstmodel.ts.cuttingLevel clstmodel.ts.cuttingLevel],'k'); 
%         else
%             contour( Ax, Ay, reshape(dist,100,100),[supportmodel.r supportmodel.r],'k'); 
%         end
        contour( Ax, Ay, reshape(dist,100,100),[supportmodel.r supportmodel.r],'k'); 
        title(strcat(['Support Type: ', supportmodel.support_type,', Labelling method: ', clstmodel.labeling_method]));
        hold off
    else
        figure; hold on;
        test.X=input;
        test.y=clstmodel.cluster_labels;
        ppatterns2(test);
%         contour( Ax, Ay, reshape(dist,100,100),[supportmodel.r supportmodel.r],'k');
        title(strcat(['Support Type: ', supportmodel.support_type,', Labelling method: ', clstmodel.labeling_method]));
        hold off
    end
    
end
