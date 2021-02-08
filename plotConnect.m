function plotConnect(connectMatrix,siteLocations,varargin)
% PLOTCONNECTIONS Draw arrows representing strength of connectivity between 
% sites on a map (already plotted, using e.g. plotMeshPDensSimple.m)
%
% Inputs:   connectMatrix - the matrix of connectivities between sites
%           siteLocations - a 2d matrix of site locations (x,y in rows)
%           
%           varargin    threshold - minimum value of connectivity to plot
%           thickness - multiplier to make lines thicker
%
% Outputs:  NONE
%
    
    narginchk(1,14);
    zScale = [];
    minThreshold = 0;
    maxThreshold = 1;
    
    cutoffMin = -10;
        
    lineWidth=2;
    os=1;
    colorBar=0; % Should a colorbar be added?
    logScale=0;
    alphaV=1;
    
    % Not presently used; always added to an empty pDens plot
%     xl = [0,0];
%     yl = [0,0];
%     add=0;

    for i = 1:2:length(varargin) % only bother with odd arguments, i.e. the labels
        switch varargin{i}
            case 'logScale'
                logScale =  varargin{i+1};
                if (logScale==1)
                    minThreshold = -5;
                    maxThreshold = 0;
                end
            case 'minThreshold'
                minThreshold = varargin{i+1};
            case 'maxThreshold'
                maxThreshold = varargin{i+1};
            case 'os'
                os = varargin{i+1};
            case 'lineWidth'
                lineWidth = varargin{i+1};
            case 'colorBar'
                colorBar =  varargin{i+1};
            
            case 'zScale'
                zScale = varargin{i+1};
            case 'alphaV'
                alphaV = varargin{i+1};
%             case 'xl'
%                 xl = varargin{i+1};
%             case 'yl'
%                 yl = varargin{i+1};
%             case 'add'
%                 add = varargin{i+1};
            case 'cutoffMin'
                cutoffMin = varargin{i+1};
        end
    end
    warning('***If specifying logScale and thresholds, ensure thresholds are specified after logscale***');
    if ~isempty(zScale)
        minThreshold = zScale(1);
        maxThreshold = zScale(2);
    end
    
    % Normalise the matrix, if the values are larger than 1, and print a
    % warning
    connectMatrix2 = connectMatrix;
    maxValue=max(max(connectMatrix2));
    if (maxValue>1 && logScale==0)
        warning('Maximum connectivity value > 1; scaling matrix to [0:1]');
        connectMatrix2=connectMatrix/max(max(connectMatrix));
        %maxValue=1;
    end
    
    % Get only the non-zero connectivities
    row=[];
    col=[];
    if (logScale==0)
        [row,col]=find(connectMatrix>minThreshold);
    else
        [row,col]=find(~isinf(log10(connectMatrix)));
        %[row,col]=find(log10(connectMatrix)>minThreshold);
    end
    
    for r=1:size(row,1)
        %val=connectMatrix2(row(r),col(r));

        val=connectMatrix2(row(r),col(r));
        %disp(['raw value = ' num2str(val)]);


        if (logScale==1)
            val=log10(val);
        end
        %disp(['transformed value 1 = ' num2str(val)]);

        % Fix the maximum value to be plotted
        if (val>maxThreshold)
            val=maxThreshold;
        end
        %if (val<minThreshold)
            %val=minThreshold;
        %end

        x1=double(siteLocations(row(r),1));
        x2=double(siteLocations(col(r),1));
        y1=double(siteLocations(row(r),2));
        y2=double(siteLocations(col(r),2));

        %alpha       = 0.15;   % head length
        %beta        = 0.07;   % head width
        den         = x2 - x1 + eps;                                % make sure no division by zero occurs
        teta        = atan( (y2-y1)/den ) + pi*(x2<x1) - pi/2;      % angle of arrow
        cs          = cos(teta);                                    % rotation matrix
        ss          = sin(teta);
        R           = [cs -ss;ss cs];
        %line_length = sqrt( (y2-y1)^2 + (x2-x1)^2 );                % sizes
        if (os==0 || os==2)
            head_length=0.02;
            head_width=0.02;
        else
            head_length=1000;
            head_width=1000;
        end
        x0          = x2*cs + y2*ss;                                % build head coordinates
        y0          = -x2*ss + y2*cs;
        coords      = R*[x0 x0+head_width/2 x0-head_width/2; y0 y0-head_length y0-head_length];


        if (val>minThreshold)
            % Scale the value onto [0:1]. This will mean that: 
            % - The minimum threshold is displayed black
            % - The maximum threshold is displayed red
            val = (val-minThreshold)/(maxThreshold-minThreshold);
            %disp(['transformed value 2 = ' num2str(val)]);

            p=plot([siteLocations(row(r),1) siteLocations(col(r),1)],[siteLocations(row(r),2) siteLocations(col(r),2)],...
                'Color',[val,0,0],'LineWidth',lineWidth);
            p.Color(4)=alphaV;
            %p=plot([siteLocations(row(r),1) siteLocations(col(r),1)],[siteLocations(row(r),2) siteLocations(col(r),2)],'-k','LineWidth',linewidth);
            patch(coords(1,:),coords(2,:),[0 0 0],'FaceColor',[val 0 0],'EdgeColor',[val 0 0],'FaceAlpha',alphaV,'EdgeAlpha',0);
        else
            % Value below the minimum threshold - plot a dashed black line
            if (val>cutoffMin)
                p=plot([siteLocations(row(r),1) siteLocations(col(r),1)],[siteLocations(row(r),2) siteLocations(col(r),2)],'--k');
                p.Color(4)=alphaV;
                patch(coords(1,:),coords(2,:),[0 0 0],'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'FaceAlpha',0,'EdgeAlpha',alphaV);
            end
        end
    end
    
    if (colorBar==1)
        % Sort out the colorbar based on provided thresholds
        colVals=(0:0.25:1)';
        zer=zeros(size(colVals));
        cmap=[colVals,zer,zer];
        %cb1=colorbar('Position',[.87 .13 .05 .77]);
        if (logScale==0)
            caxis([minThreshold maxThreshold]);
            colormap(cmap)
            cb1 = colorbar;
            cb1.Label.String = 'P(connection)';
        else
            caxis([minThreshold maxThreshold]);
            colormap(cmap)
            cb1 = colorbar;
            cb1.Label.String = 'log_{10}(P(connection))';
        end
    end
end