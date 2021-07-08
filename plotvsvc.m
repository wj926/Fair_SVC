
function varargout = plotvsvc(input,clstmodel)
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
    
 
        figure; hold on;
        test.X=input;
        test.y=clstmodel.cluster_labels;
%         ppatterns2(clstmodel.local','ro',10);   % Local
%         ppatterns2(input(:,supportmodel.bsv_ind),'k.',5);  % BSV
        ppatterns2(test);
%         contour( Ax, Ay, reshape(dist,100,100),[supportmodel.r supportmodel.r],'k');
%         legend('SVs','BSV (outliers)','Location','SouthEast');
        title(strcat(['Support Type: ', supportmodel.support_type,', Labelling method: ', clstmodel.labeling_method]));
        hold off

    
else
    
        figure; hold on;
        test.X=input;
        test.y=clstmodel.cluster_labels;
        ppatterns2(test);
        contour( Ax, Ay, reshape(dist,100,100),[supportmodel.r supportmodel.r],'k');
        title(strcat(['Support Type: ', supportmodel.support_type,', Labelling method: ', clstmodel.labeling_method]));
        hold off

    
end
