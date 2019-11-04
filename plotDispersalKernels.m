function plotDispersalKernels(Mesh,particleDensity,siteLocations,nrow,ncol,varargin)
    
    maxDist=15000;
    histInterval=1000;
    
    if mod(ncol,2)==0
        xlabelSubplot=floor(ncol/2.0)+(nrow-1)*ncol;
    else
        xlabelSubplot=ceil(ncol/2.0)+(nrow-1)*ncol;
    end
    if mod(nrow,2)==0
        ylabelSubplot=(floor(nrow/2.0)-1)*ncol + 1;
    else
        ylabelSubplot=(ceil(nrow/2.0)-1)*ncol + 1;
    end
    
    xlabelSubplot
    ylabelSubplot

    for i = 1:2:length(varargin) % only bother with odd arguments, i.e. the labels
        switch varargin{i}
            case 'maxDist'
                maxDist = varargin{i+1};
            case 'histInterval'
                histInterval = varargin{i+1};
        end
    end
    
    clf
    for site=1:size(siteLocations,1)

        dens=mean(particleDensity(:,site,:),3);
        [dHist,dists,percentileDist]=dispersalHistFromDensity(Mesh,dens,[siteLocations.Easting(site),siteLocations.Northing(site)],...
            'histInterval',histInterval,'maxDist',maxDist);

        dHist(:,1)=dHist(:,1)/1000;

        subplot(nrow,ncol,site)
        h=histogram('BinEdges',[dHist(:,1);dHist(end,1)+histInterval/1000]','BinCounts',dHist(:,2)'/sum(dHist(:,2)));
        title(siteLocations.SEPA_Site(site));
        if (site==ylabelSubplot)
            ylabel('Proportion of particle time')
        end
        if (site==xlabelSubplot)
            xlabel('Distance (km)')
        end

        yMax=0.8;
        
        siteLocations.SEPA_Site(site)
        percentileDist/1000
        
        xlim([0 maxDist/1000])
        ylim([0 yMax])
        
        line([percentileDist/1000 percentileDist/1000], [0 yMax],'Color','red');

    end

end