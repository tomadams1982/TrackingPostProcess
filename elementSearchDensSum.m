function [dSum,dMean,dCount] = elementSearchDensSum(locationsDir,startDate,endDate,searchElements,varargin)
% ELEMENTCOUNTS  count the hours that a particular set of elements is occupied by a
% particles during a set time window (optionally filtering by source site)
%
% Making a connectivity calcluation by this method can underestimate in
% comparison to the arrival count produced by the particle tracking model
% for two reasons:
% - Not all particles going within specified distance are picked up, as it is
% based on element ID. Elements are only included in the search if their
% centroid is within the distance (note that this can also lead to
% particles being counted if they enter an element which has a centroid
% within the threshold, even if they do not pass within the threshold
% themselves)
% - This is only done based on the hourly text files that are written out,
% but the particle tracking model counts particles that go within the
% radius at any moment (on any of the sub-hourly timesteps).

%   Inputs: locationsDir -
%           startDate - start date for elementCounts in format datetime
%           endDate - end date for elementCounts in format datetime
%           varargin:
%               startTime - start time in hours
%               endTime - start time in hours
%               status - numeric status of particles to count (main categories:
%                   1 = free, 2 = viable (able to settle))
%               sourceLocations - list of site IDs to count particles from
%               extension - the file extension to look within
%               sourceLocations - a list of source site IDs from which to
%                   count particles
%               firstAppearanceOnly - Should the first appearance only be
%                   counted (the default), or should subsequent values be
%                   counted?
%               printIDs - Should the particle IDs be printed to screen
%                   (debugging)
%               firstElementIndex: Allow switching between Java/C 
%                   arrays (first index 0) and Matlab/R arrays (first 
%                   index 1)
%
%   Outputs:    counts - a matrix containing [elementIndices,elementCounts, elementDensities] 
%

    % Has the user given the right number of inputs?
    narginchk(4, 14)
    
    startTime=0;
    endTime=24;
    status=0;
    extension='.dat';
    sourceLocations=[];
    firstAppearanceOnly=1;
    printIDs=0;
    firstElementIndex=0;
    
    for i = 1:2:length(varargin) % only bother with odd arguments, i.e. the labels
        switch varargin{i}
            case 'startTime'
                startTime = varargin{i+1};
            case 'endTime'
                endTime = varargin{i+1};
            case 'status'
                status = varargin{i+1};
            case 'sourceLocations'
                sourceLocations = varargin{i+1};
            case 'extension'
                extension = varargin{i+1};
            case 'first'
                firstAppearanceOnly = varargin{i+1};
            case 'printIDs'
                printIDs = varargin{i+1};
            case 'firstElementIndex'
                firstElementIndex = varargin{i+1};
        end
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
    
    startOffset=1;
    if firstElementIndex==1
        startOffset=0;
    end

    locationsDir;
    
    % load files for the days required
    % Check that output for the required days exists
    %locfiles=dir('locations*.out');
    locfiles=dir([locationsDir '/locations*' extension]);
    % pLocations will contain the files to be used in the element count
    pLocations = {};
    daysFound=0;
    for date=startDate:endDate
        for i=1:size(locfiles)
            if (strcmp(['locations_',datestr(date,'yyyymmdd'),extension],locfiles(i).name)==1)
                %fprintf([locfiles(i).name '\n']);
                daysFound=daysFound+1;
                pLocations{daysFound}=fullfile(locationsDir,locfiles(i).name);
            end
        end
    end
    
    pLocations';
    
    % Here, pLocations is a cell array of filenames
    dsLocations = datastore(pLocations,'Type','tabulartext','FileExtensions',extension,'ReadVariableNames',1);
    dsLocations.Delimiter = ' ';
    dsLocations.NumHeaderLines = 0;
    dsLocations.TextscanFormats = {'%d','%d','%{yyyyMMdd}D','%f','%s','%f','%f','%d','%d','%f'};
    dsLocations.ReadSize= 500000;
    dsLocations.SelectedVariableNames = {'ID','startLocation','elem','status','density'};
    
    % The loop is very slightly faster than reading everything in one go,
    % for two days, and eliminates the huge array
    % https://blogs.mathworks.com/loren/2014/12/03/reading-big-data-into-matlab/
    %tic
    counter = 0;
    reset(dsLocations)
    
    if size(sourceLocations,1) > 1
        error('elementSearch.m: either do not specify source, or select a single source location');
    end
    
    while hasdata(dsLocations)
        % Read in Chunk
        locationDataChunk = read(dsLocations);
        counter = counter + 1;
        
        % Count things only from the desired startLocations
        if ~isempty(sourceLocations)
            locationDataChunk = locationDataChunk(ismember(locationDataChunk.startLocation,sourceLocations),:);
        end 
        
        if ~isempty(searchElements)
            locationDataChunk = locationDataChunk(ismember(locationDataChunk.elem+startOffset,searchElements),:);
        end 
        % Count things only with the desired status
        if status~=0
            locationDataChunk = locationDataChunk(ismember(locationDataChunk.status,status),:);
        end
        
        % print out unique sources
%         if (~isempty(locationDataChunk))
%             unique(locationDataChunk.startLocation);
%         end
        
        % print out first (highest density) occurence of each unique particle ID within the elements, 
        % calculate sum, mean and count of density values
        if (~isempty(locationDataChunk))
            d=[];
            if firstAppearanceOnly==1
                ids=unique(locationDataChunk.ID);
                if (printIDs)
                    ids
                end
                for i=1:length(ids)
                    d(i)=locationDataChunk.density(find(locationDataChunk.ID==ids(i),1,'first'));
                end
            else
                d=locationDataChunk.density;
            end
            dSum=sum(d);
            dMean=mean(d);
            dCount=length(d);
        else
            dSum=0;
            dMean=0;
            dCount=0;
        end        
    end
    
    %counter
    %toc
end









