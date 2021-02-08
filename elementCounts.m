function [counts,sourceCounts] = elementCounts(mesh,locationsDir,startDate,endDate,varargin)
% ELEMENTCOUNTS  count the hours that each model element is occupied by a
% particle during a set time window (optionally scaled by recorded density)
%
%   Inputs: mesh - model mesh
%           locationsDir - directory containing the locations files
%           startDate - start date for elementCounts in format datetime
%           endDate - end date for elementCounts in format datetime
%           varargin:
%               startTime - start time in hours
%               endTime - start time in hours
%               status - numeric status of particles to count (main categories:
%                   1 = free, 2 = viable (able to settle))
%               sourceLocations (or locationsToCount) - list of site IDs to count particles from
%               scaling -
%               firstElementIndex: Allow switching between Java/C 
%                   arrays (first index 0) and Matlab/R arrays (first 
%                   index 1)
%
%   Outputs:    counts - a matrix containing [elementIndices,elementCounts, elementDensities]
% 				sourceCounts -	 
%

    % Has the user given the right number of inputs?
    narginchk(4, 12)
    
    startTime=0;
    endTime=24;
    status=0;
    extension='.dat';
    sourceLocations=[];
    scaling=1;
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
            case 'locationsToCount'
                sourceLocations = varargin{i+1};
            case 'extension'
                extension = varargin{i+1};
            case 'scaling'
                scaling = varargin{i+1};
            case 'firstElementIndex'
                firstElementIndex = varargin{i+1};
        end
    end
    
    sourceLocations;
    
    % Catch the case where mesh has been read in as mesh=load('mesh.mat');
    if(isfield(mesh,'mesh'))
        mesh=mesh.mesh;
    end
    
    % Check the inputs. (see
    % C:\Users\SA01TA\Documents\code\matlab\class_SAMS\class_SAMS\work\Days01_02\fitTrendModel.m)
    % Presently not checking "mesh" on input. Is there a smart way to do
    % this rather than allowing a fail later?
    validateattributes(startDate, {'datetime'}, ...
        {'vector'}, ...
        mfilename, 'the startDate', 2)
    if endDate<startDate
        error('elementCounts.m: endDate supplied (arg 4) was not later startDate (arg 3)');
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
        
    validateattributes(firstElementIndex, {'double'}, ...
        {'>=',0,'<=',1}, ...
        mfilename, 'the firstElementIndex')
    
    startOffset=1;
    if firstElementIndex==1
        startOffset=0;
    end

    % load files for the days required
    % Check that output for the required days exists
    locationsDir;
    %disp([locationsDir '/locations*' extension])
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
    
    counts=zeros(size(mesh.uvnode,1),3);
    counts(:,1)=1:size(mesh.uvnode,1);
    sourceCounts=zeros(size(mesh.uvnode,1),2);
    sourceCounts(:,1)=1:size(mesh.uvnode,1);
    sources=cell(size(mesh.uvnode,1),1);
    
    if ~isempty(pLocations)
        
        % Here, pLocations is a cell array of filenames
        dsLocations = datastore(pLocations,'Type','tabulartext','FileExtensions',extension,'ReadVariableNames',1);
        dsLocations.Delimiter = ' ';
        dsLocations.NumHeaderLines = 0;
        %dsLocations.TextscanFormats = {'%d','%d','%{yyyyMMdd}D','%f','%s','%f','%f','%d','%d','%f'};
        dsLocations.TextscanFormats = {'%d','%d','%{yyyyMMdd}D','%f','%s','%f','%f','%d','%d','%f','%d','%f','%f'};
        dsLocations.ReadSize= 500000;
        dsLocations.SelectedVariableNames = {'startLocation','elem','status','density'};

        % The loop is very slightly faster than reading everything in one go,
        % for two days, and eliminates the huge array
        % https://blogs.mathworks.com/loren/2014/12/03/reading-big-data-into-matlab/
        %tic
        counter = 0;
        reset(dsLocations)

        while hasdata(dsLocations)
            % Read in Chunk
            locationDataChunk = read(dsLocations);
            counter = counter + 1;

            % Count things only from the desired startLocations
            if isempty(sourceLocations)==false
                locationDataChunk = locationDataChunk(ismember(locationDataChunk.startLocation,sourceLocations),:);
            end 
            % Count things only with the desired status
            if status~=0
                locationDataChunk = locationDataChunk(ismember(locationDataChunk.status,status),:);
            end

            % multiply values by scaling prior to summation
            if ~isempty(scaling)
                for i=1:size(sourceLocations,1)
                    a=ismember(locationDataChunk.startLocation,sourceLocations{i});
                    locationDataChunk.density(a)=locationDataChunk.density(a)*scaling(i);
                    % check calculation
                    %fprintf('%d %s %f %f\n',i,locationsToCount{i},scaling(i),sum(locationDataChunk.density(a)));
                end
            end

            % Find the occupied elements, and count the number of sources that
            % supplied them
            occupiedElems=unique(locationDataChunk.elem)+startOffset;
            for i=1:size(occupiedElems,1)
                %sources{i}=unique(locationDataChunk.startLocation(locationDataChunk.elem==i));
                sources{occupiedElems(i)}=[sources{occupiedElems(i)};unique(locationDataChunk.startLocation(locationDataChunk.elem+startOffset==occupiedElems(i)))];
                %sources{occupiedElems(i)}
            end

            % Print out the location data that has been identified - to check
            % data split correctly
            %locationDataChunk

            % Calculate summed element counts and densities
            counts(:,2) = counts(:,2)+accumarray(locationDataChunk.elem+startOffset,1,[size(mesh.uvnode,1) 1],@sum); 
            counts(:,3) = counts(:,3)+accumarray(locationDataChunk.elem+startOffset,locationDataChunk.density,[size(mesh.uvnode,1) 1],@sum);

        end

%         for i=1:size(sourceCounts,1)
%             sourceCounts(i,2) = size(sources{i},1);
%         end
    else
        warning('No files found within requested date range');
    end
    
    %counter
    %toc
end









