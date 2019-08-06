function [distHist1,distHist2,distances] = dispersalDistances(mesh,startLocations,startDate,endDate,startTime,endTime,status,locationsToCount,histInterval,maxDist)
% DISPERSALDISTANCES - Calculate dispersal distances from a set of source sites (rows of
% startLocations and indices in locationsToCount correspond to numeric values in data.startLocation).
%
% Note that this expects the particle location data to be formatted as rows with source sites in sequential fashion
% (that is, row 1 particle started at site 1, row 2 at site 2, etc). Numeric values of startlocations are stored but not
% read here.
%
%   Inputs: mesh    - model mesh; this would only be needed if calculating
%                       seaway distance (much more computationally
%                       intensive)
%           startLocations - Table containing information on start locations
%                       of particles
%           startDate - start date for elementCounts in format datetime
%           endDate - end date for elementCounts in format datetime
%           startTime - start time in hours
%           endTime - start time in hours
%           status  - numeric status of particles to count (main categories:
%                       1 = free, 2 = viable (able to settle))
%           locationsToCount - list of site IDs to count particles from
%           histInterval - interval for histogram (m) (omit to prevent creation
%                       of histogram)
%           maxDist     - maximum distance (m) bin for histogram
%
%   Outputs:    distHist1    - Histogram of distances, mortality adjusted
%               distHist2    - Histogram of distances, not mortality adjusted
%               distances   - particle distances from each source site listed in locationsToCount
%                 

    narginchk(4, 10)
    if nargin < 10
        if nargin < 8
            if nargin < 7
                if nargin < 6
                    % Issue a custom warning, including the ID.
                    warning('dispersalDistances:startTime', 'Fewer than 6 arguments; setting startTime to 0')
                    startTime = 0;
                    warning('dispersalDistances:endTime', 'Fewer than 6 arguments; setting endTime to 24')
                    endTime = 0;
                end
                warning('dispersalDistances:status', 'Fewer than 7 arguments; setting status to all')
                status = 0;
            end
            warning('dispersalDistances:default startlocs', 'Fewer than 8 arguments; setting locations counted to all')
            locationsToCount = [];
        end
        warning('dispersalDistances: no histogram', 'Fewer than 10 arguments, no histogram produced (must be done internally to include density)')
    end
    
    % Check the inputs. (see
    % C:\Users\SA01TA\Documents\code\matlab\class_SAMS\class_SAMS\work\Days01_02\fitTrendModel.m)
    % Presently not checking "mesh" on input. Is there a smart way to do
    % this rather than allowing a fail later?
    validateattributes(startDate, {'datetime'}, ...
        {'vector'}, ...
        mfilename, 'the startDate', 2)
    if endDate<startDate
        error('elementCounts.m: endDate supplied (arg 3) was not later startDate (agr 2)');
    end
%     validateattributes(endDate, {'datetime'}, ...
%         {'>=',startDate,'vector'}, ...
%         mfilename, 'the endDate', 3)
    validateattributes(startTime, {'double'}, ...
        {'>=',0,'<=',24}, ...
        mfilename, 'the startTime', 4)
    validateattributes(endTime, {'double'}, ...
        {'>=',0,'<=',24}, ...
        mfilename, 'the endTime', 5)
    validateattributes(status, {'double'}, ...
        {'scalar'}, ...
        mfilename, 'the status', 6)
    %validateattributes(locationsToCount, {'cell'}, ...
    %    mfilename, 'the locationsToCount', 7)
%     validateattributes(scaling, {'double'}, ...
%         {'size',[size(locationsToCount) 1]}, ...
%         mfilename, 'the scaling', 8)
    
    % load files for the days required
    % Check that output for the required days exists
    locfiles=dir('locations*.out');
    % pLocations will contain the files to be used in the element count
    pLocations = {};
    daysFound=0;
    for date=startDate:endDate
        for i=1:size(locfiles)
            if (strcmp(['locations_',datestr(date,'yyyymmdd'),'.out'],locfiles(i).name)==1)
                %fprintf([locfiles(i).name '\n']);
                daysFound=daysFound+1;
                pLocations{daysFound}=fullfile(pwd,locfiles(i).name);
            end
        end
    end
    
    % Here, pLocations is a cell array of filenames
    dsLocations = datastore(pLocations,'Type','tabulartext','FileExtensions','.out','ReadVariableNames',1);
    dsLocations.Delimiter = ' ';
    dsLocations.NumHeaderLines = 0;
    dsLocations.TextscanFormats = {'%d','%d','%{yyyyMMdd}D','%f','%s','%f','%f','%d','%d','%f'};
    dsLocations.ReadSize= size(startLocations,1)*1000;
    dsLocations.SelectedVariableNames = {'hour','startLocation','x','y','status','density'};
    
    % The loop is very slightly faster than reading everything in one go,
    % for two days, and eliminates the huge array
    % https://blogs.mathworks.com/loren/2014/12/03/reading-big-data-into-matlab/
    %tic
	
	if isempty(locationsToCount)
		locationsToCount=1:size(startLocations,1)
	end
    counter = 0;
    rowsCounted = 0;
    distances=[];
	
    if nargin==10
        distHist1=zeros(maxDist/histInterval,1);
        distHist2=zeros(maxDist/histInterval,1);
    end
    
    reset(dsLocations)
        
    % The code below is what actually makes the calculations
    while hasdata(dsLocations)
        % Read in Chunk
        locationDataChunk = read(dsLocations);
        counter = counter + 1;
        rowsCounted = rowsCounted + dsLocations.ReadSize;
        
        % Filter by selected sites 
        %rows=find(ismember(locationDataChunk.startLocation,startLocations(locationsToCount,:).SEPA_Site));
        %locationDataChunk=locationDataChunk(rows,:);
        
        for i=1:size(locationsToCount,1)
            rows=find(ismember(locationDataChunk.startLocation,startLocations(locationsToCount(i),:).SEPA_Site));
            
            % Get the distance from source location
            d=pdist2([startLocations.Easting(locationsToCount(i)) startLocations.Northing(locationsToCount(i))],...
                [locationDataChunk.x(rows) locationDataChunk.y(rows)]);
            
            % Get the corresponding density at this point in time
            dens=locationDataChunk.density(rows);
            
            if (nargin<10)
                % If no producing a histogram, just output distances directly
                distances=[distances,d];
            else
                % Place a value in the relevant histogram bin, scaled by
                % density
                bin=floor(d/histInterval+1);
                for j=1:length(bin)
                    if bin(j)<length(distHist1)
                        distHist1(bin(j))=distHist1(bin(j))+dens(j);
                        distHist2(bin(j))=distHist2(bin(j))+1;
                    else
                        distHist1(end)=distHist1(end)+dens(j);
                        distHist2(end)=distHist2(end)+1;
                    end
                end
            end            
            
        end    
    end
end