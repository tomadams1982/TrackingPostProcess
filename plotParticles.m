function plotParticles(particleTable,varargin)
    
    plotTracks=1;
    for i = 1:2:length(varargin) % only bother with odd arguments, i.e. the labels
        switch varargin{i}
            case 'plotTracks'
                plotTracks = varargin{i+1};
        end
    end

    if (plotTracks==0)
        scatter(particleTable.x,particleTable.y,'.');
    else
        ps = unique(particleTable.ID);
        for p = 1:length(ps)
            Tp = particleTable(particleTable.ID==ps(p),:);
            plot(Tp.x,Tp.y);
        end
    end
end