function plotConnections(connectMatrix,siteLocations,threshold,thickness)
% PLOTCONNECTIONS Draw arrows representing strength of connectivity between 
% sites on a map (already plotted, using e.g. plotMeshPDensSimple.m)
%
% Inputs:   connectMatrix - the matrix of connectivities between sites
%           siteLocations - a 2d matrix of site locations (x,y in rows)
%           threshold - minimum value of connectivity to plot
%           thickness - multiplier to make lines thicker
%
% Outputs:  NONE
%
    
    % A normalised matrix
    connectMatrix2=connectMatrix/max(max(connectMatrix));
    
    
    [row col]=find(connectMatrix>threshold);

    for r=1:size(row,1)
        val=connectMatrix2(row(r),col(r));
        
        minwidth=2;
        linewidth=max(minwidth,thickness*log(val+1));
        % add something here to change color of lines plotted
        
        p=plot([siteLocations(row(r),1) siteLocations(col(r),1)],[siteLocations(row(r),2) siteLocations(col(r),2)],'Color',[tanh(20*val),0,0],'LineWidth',linewidth);
        %p=plot([siteLocations(row(r),1) siteLocations(col(r),1)],[siteLocations(row(r),2) siteLocations(col(r),2)],'-k','LineWidth',linewidth);
        
        x1=double(siteLocations(row(r),1));
        x2=double(siteLocations(col(r),1));
        y1=double(siteLocations(row(r),2));
        y2=double(siteLocations(col(r),2));

        alpha       = 0.15;   % head length
        beta        = 0.07;   % head width
        den         = x2 - x1 + eps;                                % make sure no devision by zero occurs
        teta        = atan( (y2-y1)/den ) + pi*(x2<x1) - pi/2;      % angle of arrow
        cs          = cos(teta);                                    % rotation matrix
        ss          = sin(teta);
        R           = [cs -ss;ss cs];
        line_length = sqrt( (y2-y1)^2 + (x2-x1)^2 );                % sizes
        head_length=0.01;
        head_width=0.01;
        x0          = x2*cs + y2*ss;                                % build head coordinates
        y0          = -x2*ss + y2*cs;
        coords      = R*[x0 x0+head_width/2 x0-head_width/2; y0 y0-head_length y0-head_length];
        patch( coords(1,:),coords(2,:),[0 0 0] );

        
    end
end