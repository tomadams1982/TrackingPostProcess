function [connectivity] = connectivityFromElementDensity(Mesh,siteLocations,weekStartDates,directories,varargin)
%CONNECTIVITYFROMELEMENTDENSITY - Calculate connectivity between sites from
%   the element density files, based on a list of site location IDs, which
%   will be searched for in the files.
%
    nparts=5*ones(size(siteLocations,1),1); % The number of particles released per site per hour
    nhours=168; % The number of hours to be used for aggregation of values - NOTE calculation is presently fixed to aggregate weekly
    % List elements within a certain radius
    radius=500;
    destinationFile=[];
    first=1; % Record the first appearance only, default for connectivity calculation. If want to record overall pressure, =0 might be better
    status=2; % For sea lice, use copepodid densities only
    
    extensions=cell(size(siteLocations,1),1);
    for i = 1:length(extensions)
        extensions{i} = '.dat';
    end
    
    sources=1:size(siteLocations,1); % Numerical list of IDs of source sites
    destinations=1:size(siteLocations,1); % Numerical list of IDs of destination sites

    for i = 1:2:length(varargin) % only bother with odd arguments, i.e. the labels
        switch varargin{i}
            case 'sources'
                sources = varargin{i+1};
            case 'destinations'
                destinations = varargin{i+1};
            case 'extensions'
                extensions = varargin{i+1};
            case 'nparts'
                nparts = varargin{i+1};
%             case 'nhours' % DON'T ALLOW THIS TO BE CHANGED PRESENTLY as
%             aggregation is fixed at one week
%                 nhours = varargin{i+1};
            case 'radius'
                radius = varargin{i+1};
            case 'firstAppearance'
                first = varargin{i+1};
            case 'nParts'
                nparts = varargin{i+1};            
        end
    end
    
    connectivity=zeros(size(siteLocations,1),size(siteLocations,1),length(weekStartDates));

%     locationsDir(1:16)={'locations_filtered/'};
%     locationsDir(17:20)={'locations_filtered/dfClydeSites_160215/'};
%     locationsDir(21)={'W:/sa01ta/fishfarm_long/dfClydeExtra/'};
%     extension(1:16)={'.out'};
%     extension(17:21)={'.dat'};

    % Check for sources of particles going in those elements
    % Get the sum of the densities, divide by nparts=5 and by nhours=168
    for i=1:length(sources)
        tic
        disp(['SOURCE ' num2str(sources(i))]);
        %startlocsSSF(sourceSite,:);
        
        for j=1:length(destinations)
            disp(['DESTINATION ' num2str(destinations(j))]);
            %startlocsSSF(destSite,:);
            destLoc=[siteLocations.Easting(destinations(j)),siteLocations.Northing(destinations(j))];
            near=find(pdist2(destLoc,Mesh.uvnode_os)<radius);

            % Check for sources of particles going in those elements
            % Get the sum of the densities, divide by nparts=5 and by nhours=168
            dSum=zeros(length(weekStartDates),1);
            parfor t=1:length(weekStartDates)
                %weekStartDate=startDate+(t-1)*7;
                weekStartDates(t);
                %endDate=startDate;
                endDate=weekStartDates(t)+6;
                %fprintf('Week %d - %s to %s (inclusive)\n',t,datestr(weekStartDate,'yyyymmdd'),datestr(endDate,'yyyymmdd'));
                fprintf('%d ',t);

                %allLocations.SEPA_Site(sources(i))
                %sources(i)
                %directories{sources(i)}
                
                [dSum(t),dMean,dCount]=elementSearchDensSum(directories{sources(i)},weekStartDates(t),endDate,near,...
                    'sourceLocations',siteLocations.SEPA_Site(sources(i)),'extension',extensions{sources(i)},'status',2,'first',first);
                %disp([num2str(dSum/(nparts*nhours)) ' ' num2str(dMean) ' ' num2str(dCount)]);

            end
            connectivity(sources(i),destinations(j),:) = dSum/(nparts(i)*nhours);

            fprintf('\n');
        end
        % Save work in progress
        try
            save('connectivityTmp.mat','connectivity')
        catch
            warning('Could not save interim file');
        end
        toc
    end

end