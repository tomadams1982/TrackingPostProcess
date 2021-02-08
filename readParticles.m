function [particleTable] = readParticles(datadir,varargin)

    maxDays=1;
    for i = 1:2:length(varargin) % only bother with odd arguments, i.e. the labels
        switch varargin{i}
            case 'maxDays'
                maxDays = varargin{i+1};
        end
    end
            
    %datadir = 'C:\Users\sa01ta\OneDrive - SAMS\Documents\OFF-AQUA\liceTracking\liceTrack210113_maxMerged\';
    fileList=dir([datadir '/locations_*'])
    
    if (maxDays>1)
        md=maxDays;
    else
        md=length(fileList);
    end
    if (md>length(fileList))
        md=length(fileList);
        disp(['Number of days available = ',num2str(md)]);
    end
    
    particleTable=readtable([datadir,fileList(1).name]);
    if (md>1)
        for f=1:md
            fprintf(['-- ' num2str(f) ' --\n']);
            newDay=readtable([datadir,fileList(f).name]);
            particleTable=[particleTable;newDay];
        end
    end
end