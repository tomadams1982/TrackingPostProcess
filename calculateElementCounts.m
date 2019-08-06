function [countsJuvenile,countsCopepodid] = calculateElementCounts(locationsDir,sourceLocations,Mesh,weekStartDates,varargin)
%CALCULATEELEMENTCOUNTS Sum element densities from particle tracking locations files.
%   This function reads model output files, and calculates a weekly sum of
%   the particle densities in each hydrodynamic model mesh element, split
%   by source site. It does this one source site at a time in order to make
%   an array which has source site as the second index.
%
% Inputs:   locationsDir    - the directory containing the locations files
%           sourceLocations - a cell array of strings of the source site
%                               IDs (listed in a column in locations files)
%           Mesh            - hydrodynamic model mesh
%           weekStartDates  - datetime array of start dates (i.e. list of
%                               weeks)
%           varargins:
%               multiplier  - a scaling to apply during the calculation of
%                               particle density (default 1)
%               extension   - extension of the locations files (defaults to
%                               .dat)
%
% Outputs:  countsJuvenile  - counts of particles in the juvenile
%                               (non-infective) stage 1, aggregated by 
%                               element, source site, and week (Nmodelelements X
%                               NsourceSites X Nweeks)
%           countsCopepodid - counts of particles in the copepodid
%                               (infective) stage 2, aggregated by element
%                               source site, and week (Nmodelelements X
%                               NsourceSites X Nweeks)

    narginchk(4,8)
    
    multiplier=1;
    extension='.dat';
    
    for i = 1:2:length(varargin) % only bother with odd arguments, i.e. the labels
        switch varargin{i}
            case 'extension'
                extension = varargin{i+1};
            case 'multiplier'
                multiplier = varargin{i+1};
        end
    end
    
    validateattributes(locationsDir, {'char'}, ...
        {'scalartext'}, ...
        mfilename, 'the locationsDir', 1)
    validateattributes(sourceLocations, {'cell'}, ...
        {'vector'},...
        mfilename, 'the sourceLocations', 2)
    validateattributes(Mesh, {'struct'}, ...
        {'scalar'},...
        mfilename, 'the Mesh', 3)
    validateattributes(weekStartDates, {'datetime'}, ...
        {'vector'}, ...
        mfilename, 'the startDate', 4)
    validateattributes(multiplier, {'double'}, ...
        {'scalar'}, ...
        mfilename, 'the multiplier', 5)
    validateattributes(extension, {'char'}, ...
        {'scalartext'}, ...
        mfilename, 'the extension', 5)

    countsJuvenile=zeros(size(Mesh.uvnode,1),size(sourceLocations,1),length(weekStartDates));
    countsCopepodid=zeros(size(Mesh.uvnode,1),size(sourceLocations,1),length(weekStartDates));

    for t=1:length(weekStartDates)

        endDate=weekStartDates(t)+6;
        fprintf('Week %d - %s to %s (inclusive)\n',t,datestr(weekStartDates(t),'yyyymmdd'),datestr(endDate,'yyyymmdd'));

        parfor site=1:length(sourceLocations)

            sourceLocSingle = sourceLocations(site);
            disp(char(sourceLocSingle))
            
            % Count the juveniles
            [c1,~]=elementCounts(Mesh,locationsDir,weekStartDates(t),endDate,'extension',extension,'status',1,'sourceLocations',sourceLocSingle,'scaling',multiplier);
            countsJuvenile(:,site,t)=c1(:,3); 
            % Count the copepodids
            [c2,~]=elementCounts(Mesh,locationsDir,weekStartDates(t),endDate,'extension',extension,'status',2,'sourceLocations',sourceLocSingle,'scaling',multiplier);
            countsCopepodid(:,site,t)=c2(:,3); 

        end     
    end


end

