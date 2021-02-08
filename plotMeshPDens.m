function [out,z2,cb1] = plotMeshPDens(mesh,varargin)
% PLOTMESHPDENS  plot an FVCOM model mesh. The plot can include some 
% kind of density information (specified in "meshDensity"), site location points
% (specified in "startLocations"). Various options are available, depending
% on the plot that is desired, see Arguments.
%
% Arguments:    mesh
%               varargin:
%                   meshDensity: density information; must have same
%                       number of values as there are model elements
%                   startLocations: a list of points coordinates to add to 
%                       the plot (dimensions N x 2)
%                   os: logical, is the plot in ordnance survey (or other
%                       projected coordinate system in metres) (=1) or to
%                       be in degrees (long/lat) (=0)
%                   logScale: logical, is density information to be
%                       plotted on a logarithmic scale?
%                   xl: x limits
%                   yl: y limits
%                   add: logical, add to existing plot?
%                   areaScale: logical, is density information to be
%                       scaled by element area?
%                   colorBar: logical, display a colorbar?
%                   colorBarLabel: logical, display a colorbar label?
%                   zScale: limits for z range
%                   firstElementIndex: Allow switching between Java/C 
%                       arrays (first index 0) and Matlab/R arrays (first 
%                       index 1)
%
% Author: Thomas Adams, Scottish Association for Marine Science, 2019.
%

    narginchk(1,19);

    meshdens = [];
    startlocs = [];
    xl = [0,0];
    yl = [0,0];
    os = 2;
    logScale = 0;
    filename = 'noprint';
    add=0;
    areaScale=1;
    colorBar=0;
    colorBarLabel='density (m^{-3})';
    zScale = [];
    firstElementAdjust=0;
    plotEdges = 0;
    removeZeros = 1;
    faceAlpha = 0.8;
    bgColour=[0.6 0.6 0.6];
    minZ=[];
    cMap=[];

    for i = 1:2:length(varargin) % only bother with odd arguments, i.e. the labels
        switch varargin{i}
            case 'meshDensity'
                meshdens = varargin{i+1};
            case 'meshDensitySparse'
                meshDensSparse = varargin{i+1};
                meshdens = zeros(size(mesh.uvnode,1),1);
                meshdens(meshDensSparse(:,1)+1) = meshDensSparse(:,2);
                %find(meshdens)
            case 'startLocations'
                startlocs = varargin{i+1};
            case 'os'
                os = varargin{i+1};
            case 'logScale'
                logScale = varargin{i+1};
            case 'xl'
                xl = varargin{i+1};
            case 'yl'
                yl = varargin{i+1};
            case 'add'
                add = varargin{i+1};
            case 'areaScale'
                areaScale = varargin{i+1};
            case 'colorBar'
                colorBar = varargin{i+1};
            case 'colorBarLabel'
                colorBarLabel = varargin{i+1};
            case 'cMap'
                cMap = varargin{i+1};
            case 'zScale'
                zScale = varargin{i+1};
            case 'firstElementIndex'
                firstElementAdjust = varargin{i+1};
            case 'plotEdges'
                plotEdges = varargin{i+1};
            case 'removeZeros'
                removeZeros = varargin{i+1};
            case 'FaceAlpha'
                faceAlpha = varargin{i+1};
            case 'BackgroundColour'
                bgColour = varargin{i+1}; 
            case 'minZ'
                minZ = varargin{i+1};
            
        end
    end
    
    if ~isempty(cMap)
        colormap(cMap);
    end
        
        
    
    % Work out element areas for scaling
    if (os==1)
        nodexy=mesh.nodexy_os;
        uvnode=mesh.uvnode_os;
    elseif (os==2)
        nodexy=mesh.nodexy;
        uvnode=mesh.uvnode;
    else
        nodexy=mesh.nodexy_deg;
        uvnode=mesh.uvnode_deg;
    end
    
    col=0.6;
    
    if (add==0)
        f1=patch('Vertices',nodexy,'Faces',mesh.trinodes,'FaceVertexCdata',0,'FaceColor',bgColour);

        
        if plotEdges==1
            set(f1,'EdgeColor',[0.7 0.7 0.7],'FaceAlpha',0.3)
        else
            set(f1,'EdgeColor','none','FaceAlpha',0.3)
        end
    end
        
    if ~isempty(meshdens)
        
        if firstElementAdjust==1
            meshdens=[0,meshdens(1:(end-1))];
        end
        
        if (areaScale==1)
            % Distances above are in m. So areas are in m^2.
            % NOTE: areas are always calculated using the OS coordinates, as these
            % are in metres
            x = [mesh.nodexy_os(mesh.trinodes(:,1),1) mesh.nodexy_os(mesh.trinodes(:,2),1) mesh.nodexy_os(mesh.trinodes(:,3),1)];
            y = [mesh.nodexy_os(mesh.trinodes(:,1),2) mesh.nodexy_os(mesh.trinodes(:,2),2) mesh.nodexy_os(mesh.trinodes(:,3),2)];
            area = 0.5*abs((x(:,1)-x(:,3)).*(y(:,2)-y(:,1))-(x(:,1)-x(:,2)).*(y(:,3)-y(:,1)));
            % Debugging
            %fprintf('element areas');
            %min(area)
            %max(area)
            z=meshdens./area;            
        else
            z=meshdens;
            area=ones(length(meshdens),1);
        end
        
        if (removeZeros==1)
            z(z==0)=nan;
        end

        out=[meshdens(~isnan(z)), area(~isnan(z)), z(~isnan(z))];
        
        
        %z
        
        % Debugging
        %fprintf('density values');
        %min(z)
        %max(z)
        %size(find(~isnan(z)))
             
        % plot the mesh data, only for those cells with non-zero values

        if ~isempty(minZ)
            z(z<minZ)=nan;
        end
        
        if logScale==1
            f2=patch('Vertices',nodexy,'Faces',mesh.trinodes,'FaceVertexCdata',log10(z),'FaceColor','flat');
        else
            f2=patch('Vertices',nodexy,'Faces',mesh.trinodes,'FaceVertexCdata',z,'FaceColor','flat');
        end
        
        %set(f2,'EdgeColor','none','FaceAlpha',0)
        set(f2,'EdgeColor','none','FaceAlpha',faceAlpha);
        
        z2=z(uvnode(:,1)>xl(1) & uvnode(:,1)<xl(2) & uvnode(:,2)>yl(1) & uvnode(:,2)>xl(2));
        % Debugging
        %size(z2)
        %max(z2)
    else
        warning('meshdens not present, or mis-specified Name parameter in arguments');
    end

    hold on
    
    % add scatter points
    if ~isempty(startlocs)
        scatter(startlocs(:,1),startlocs(:,2),20,'filled')
    end

    if (colorBar>0)
        % Sort out the colorbar based on provided thresholds
        if (logScale==0)
            if isempty(zScale)
                caxis([min(z) max(z)]);
            else
                caxis(zScale);
            end
            %colormap(cmap)
            if colorBar==1
                cb1 = colorbar('Position',[.87 .13 .05 .77]);
                cb1.Label.String = colorBarLabel;
            elseif colorBar==2
                cb1 = colorbar('Position',[.90 .13 .03 .77]);
                cb1.Label.String = colorBarLabel;
            end
        else
            if isempty(zScale)
                caxis([-8 1.5]);
            else
                caxis(zScale);
            end
            %colormap(cmap)
            cb1 = colorbar('Position',[.87 .13 .05 .77]);
            cb1.Label.String = 'log_{10}(density) (m^{-3})';
        end
    end
    
    if (os==1)
        xlabel('Easting (m)');
        ylabel('Northing (m)');
    else
        xlabel('Longitude');
        ylabel('Latitude');
    end

    
    if (xl~=0)
        set(gca,'xlim',xl,'ylim',yl)
    end
    
    % print to file
    if(strcmp(filename,'noprint')~=true)
        print('-painters','-dpng','-r600',filename)
    end

end