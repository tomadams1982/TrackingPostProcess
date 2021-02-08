function plotConnectivityMatrices(connectivity,siteLocations,varargin)
% PLOTCONNECTIVITYMATRICES - Plot connectivity matrices, either linear or
% log or both

    linear=1;
    logarithmic=1;
    linearOnly=0;
    logOnly=0;
    startInd=1;
    endInd=size(connectivity,3);
    linearScale=[0 0.5];
    logScale=[-6 0];
    
    for i = 1:2:length(varargin) % only bother with odd arguments, i.e. the labels
        switch varargin{i}
            case 'startInd'
                startInd = varargin{i+1};
            case 'endInd'
                endInd = varargin{i+1};
            case 'logOnly'
                logOnly = varargin{i+1};
            case 'linearOnly'
                linearOnly = varargin{i+1};
            case 'logScale'
                logScale = varargin{i+1};
            case 'linearScale'
                linearScale = varargin{i+1};
        end
    end
    
    logOnly
    
    nColours=length(logScale(1):logScale(2));
    
    if linearOnly==1
        logarithmic=0;
    end
    if logOnly==1
        linear=0;
    end

    %[B,I]=sort(allLocations.Latitude(1:6));
    I=1:size(connectivity,1);
    mx=1;
    
    %names=startlocs.SEPA_Site(I);
    if strcmp(class(siteLocations),'cell')
        names=siteLocations;
    elseif any(strcmp('SEPA_Site', siteLocations.Properties.VariableNames))
        names=siteLocations.SEPA_Site(I);
    elseif any(strcmp('Var1', siteLocations.Properties.VariableNames))
        names=siteLocations.Var1;
    else
        names={I};
    end

    vals1=mean(connectivity(I,I,startInd:endInd),3)/mx;
    vals2=log10(mean(connectivity(I,I,startInd:endInd),3));
    vals2(isinf(vals2))=NaN;

    % Absolute relative values
    if linear==1 && logarithmic==1
        subplot(1,2,1)
    end

    if linear==1
        cmap=flipud(gray(nColours));
        cmap=cmap(2:end,:);
        colormap(cmap)
        imagesc(vals1,'AlphaData',~isnan(vals2),linearScale);

    %     colormap(flipud(gray))
    %     imagesc(mean(connectivity(I,I,1:nSimWeeks),3)/mx)

        set(gca,'Ydir','Normal')
        cb1=colorbar;
        cb1.Label.String = 'P(connection)';
        xlabel('destination site')
        xticks(1:length(names))
        xticklabels(names)
        ylabel('source site')
        yticks(1:length(names))
        yticklabels(names)
        xtickangle(90)
        % Add title if making two plots
        if logarithmic==1
            title('(a)')
        end
    end
        
    if linear==1 && logarithmic==1
        subplot(1,2,2)
    end
    
    if logarithmic==1
        mx=1;

        cmap=flipud(gray(nColours));
        cmap=cmap(2:end,:);
        colormap(cmap)
        %imagesc(vals1,'AlphaData',~isnan(vals2))
        imagesc(vals2,'AlphaData',~isnan(vals2),logScale);
        %mx=max(max(mean(connectivity(I,I,:),3)));
        %imagesc(log10(mean(connectivity(I,I,1:nSimWeeks),3)/mx))

        set(gca,'Ydir','Normal')
        cb1=colorbar;
        cb1.Label.String = 'log_{10}(P(connection))';
        xlabel('destination site')
        xticks(1:length(names))
        xticklabels(names)
        ylabel('source site')
        yticks(1:length(names))
        yticklabels(names)
        xtickangle(90)
        % Add title if making two plots
        if linear==1
            title('(b)')
        end
    end

    %print('-painters','-dpng','-r600','figures\connectivityMean_textNames.png')

end