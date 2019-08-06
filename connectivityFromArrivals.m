function [connectivity,partDens,partIDs] = connectivityFromArrivals(directory,siteIDs,startDate,endDate,nparts,varargin)
% CONNECTIVITYFROMARRIVALS  calculate the connectivity between pairs of
% sites from lists of arrivals.
%
%   Inputs: directory - the directory containing the arrivals files
%           extension - extension of the arrivals files e.g. '.out', '.dat'
%           siteIDs - a list of the habitat sites (cell array of strings;
%                      can be empty if required)
%           startDate - start date for connectivity calculation ARRIVALS in format YYYYMMDD
%           endDate - end date for connectivity calculation ARRIVALS in format YYYYMMDD
%           nparts - number of particles per site per hour
%           varargin:
%               scaling - optional weighting of connectivity by origin sites
%               firstAppearanceOnly - logical; should the first
%                           appearance be counted alone (1), or all appearances(0)?
%               extension - Extension for the arrivals files; '.out' or '.dat' etc.
%               allOutput - Should individual particle densities and IDs be
%                           produced?
%
%   Outputs:    connectivity - a matrix of dimensions length(siteIDs) X length(siteIDs) 
%                           (SOURCE X DESTINATION), containing summed connection density 
%                           values, divided by nparts
%               partDens - a cell array of the same size as connectivity.
%                           Each entry is a list of the density of each
%                           individual particle that made the connection.
%                           Only produced if allOutput == 1.
%               partIDs - a cell array as the same size as the two above.
%                           Contains particle ID numbers corresponding to 
%                           the density values in partDens.
%                           Only produced if allOutput == 1.

    narginchk(5, 11)
    
    scaling=[];
    firstAppearanceOnly=1;
    extension='.dat';
    allOutput=0;
    
    for i = 1:2:length(varargin) % only bother with odd arguments, i.e. the labels
        switch varargin{i}
            case 'scaling'
                scaling = varargin{i+1};
            case 'first'
                firstAppearanceOnly = varargin{i+1};
            case 'extension'
                extension = varargin{i+1};
            case 'allOutput'
                allOutput = varargin{i+1};
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
        error('connectivityFromArrivals.m: endDate supplied (arg 3) was not later startDate (agr 2)');
    end
%     validateattributes(endDate, {'double'}, ...
%         {'>=',startDate,'vector', 'real', 'nonnan', 'finite'}, ...
%         mfilename, 'the endDate', 3)
%     validateattributes(scaling, {'double'}, ...
%         {'size',[size(siteIDs,1) 1]}, ...
%         mfilename, 'the scaling', 8)
    
    % load files for the days required
    % Check that output for the required days exists
    arrivalfiles=dir([directory '\arrivals*' extension]);
    % pLocations will contain the files to be used in the element count
    pArrivals = {};
    daysFound=0;
    for date=startDate:endDate
        for i=1:size(arrivalfiles)
            if (strcmp(['arrivals_',datestr(date,'yyyymmdd'),extension],arrivalfiles(i).name)==1)
                %fprintf([arrivalfiles(i).name '\n']);
                daysFound=daysFound+1;
                pArrivals{daysFound}=fullfile(directory,arrivalfiles(i).name);
            end
        end
    end

    dsArrivals = datastore(pArrivals,'Type','tabulartext','FileExtensions',extension,'ReadVariableNames',1);
    dsArrivals.Delimiter = ' ';
    dsArrivals.NumHeaderLines = 0;
    dsArrivals.TextscanFormats = {'%d','%{yyyyMMdd}D','%f','%s','%{yyyyMMdd}D','%f','%s','%f','%f'};
    %dsArrivals.TextscanFormats = {'%d','%s','%s','%f'};
    dsArrivals.ReadSize= 500000;
    dsArrivals.SelectedVariableNames = {'ID','startLocation','endLocation','density'};

    counter = 0;
    connectivity=zeros(size(siteIDs,1),size(siteIDs,1));
    
    
    % Create a list of particles that have made each pairwise connection
    connectParts=cell(size(siteIDs,1),size(siteIDs,1));
    partDens=cell(size(siteIDs,1),size(siteIDs,1));
    partIDs=cell(size(siteIDs,1),size(siteIDs,1));

    reset(dsArrivals)
    while hasdata(dsArrivals)
        % Read in Chunk
        %fprintf('in data read\n');
        locationDataChunk = read(dsArrivals);
        counter = counter + 1;
        
        %locationDataChunk(ismember('KPOOL',locationDataChunk.endLocation),:)
        
        % option 1 - loop over rows in data
        if (1==1)
            for i=1:size(locationDataChunk,1)
                fromID=find(ismember(siteIDs,locationDataChunk.startLocation{i}));
                %locationDataChunk.startLocation{i};
                toID=find(ismember(siteIDs,locationDataChunk.endLocation{i}));
                %locationDataChunk.endLocation{i};
                
                % Get the particle ID
                pID = locationDataChunk.ID(i);
                
                ismember(pID,[]);
                
                % Check against the list in connectParts{} in order to
                % prevent the same particle arriving multiple times at the
                % same site
                if (~isempty(fromID) && ~isempty(toID))
                    if ~ismember(pID,connectParts{fromID,toID})
                        val=locationDataChunk.density(i)/nparts;
                        %fprintf('from %s (%d) to %s (%d) val %f\n',locationDataChunk.startLocation{i},from,locationDataChunk.endLocation{i},to,val); 
                        if ~isempty(scaling)
                            val=val*scaling(fromID);
                        end
                        connectivity(fromID,toID)=connectivity(fromID,toID)+val;
                        
                        % Add it to the list of particles that have arrived at this site
                        if firstAppearanceOnly==1
                            connectParts{fromID,toID}=[vertcat(connectParts{fromID,toID}),pID];
                        end
                        if allOutput==1
                            partDens{fromID,toID}=[vertcat(partDens{fromID,toID}),locationDataChunk.density(i)];
                            partIDs{fromID,toID}=[vertcat(partIDs{fromID,toID}),pID];
                        end
                    end
                end
                
                
            end
        end
        
        % option 2 - find relevant rows for all present sources and
        % destinations, and add densities. Possibly using accumarray
        
    end
end