function [siteValues] = extractMssSiteValues(siteStockingDataFile,siteNames,weekStartDates,varargin)
%EXTRACTMSSSITEVALUES Extract biomass for a list of sites over a specified
%time interval.
%
% Inputs:   siteStockingDataFile - the full path to the MSS data file
%                                   (normally called "se_monthly_reports.csv") 
%           siteNames       - cell array of site IDs (strings)
%           weekStartDates  - date time array of dates upon which to
%                               get the biomass (or other value)
%           varargin:
%               variable    - name of the variable to extract (if not
%                               biomass)
%
% Outputs:  siteValues      - array of values from file (Nsites X Nweeks)
%
    
    narginchk(3,5);
    
    variable='ActualBiomassOnSite_tonnes_';
    
    for i = 1:2:length(varargin) % only bother with odd arguments, i.e. the labels
        switch varargin{i}
            case 'variable'
                variable = varargin{i+1};
        end
    end
    
    validateattributes(siteStockingDataFile, {'char'}, ...
        {'scalartext'}, ...
        mfilename, 'the siteStockingDataFile', 1)
    validateattributes(siteNames, {'cell'}, ...
        {'vector'}, ...
        mfilename, 'the siteNames', 2)
     validateattributes(weekStartDates, {'datetime'}, ...
        {'vector'}, ...
        mfilename, 'the weekStartDates', 3) 
    validateattributes(variable, {'char'}, ...
        {'scalartext'}, ...
        mfilename, 'the variable', 2)

    siteStockingData=readtable(siteStockingDataFile);
    % Filter by site IDs
    %siteStockingData = siteStockingData(ismember(siteStockingData.SEPASite,siteNames),:);

    % Extract the biomasses
    siteValues = zeros(length(siteNames),length(weekStartDates));
    for site=1:size(siteNames)
        for week=1:length(weekStartDates)
            % Filter by site
            stockingDataChunk = siteStockingData(ismember(siteStockingData.SEPASite,siteNames(site)),:);
            % Filter by date
            ds=str2num(datestr(weekStartDates(week),'yyyymm'));
            stockingDs=str2num(datestr(stockingDataChunk.Year,'yyyymm'));
            ind = find(ismember(stockingDs,ds));
            if (length(ind)>1)
                warning('More than 1 row returned for this site and date');
                stockingDataChunk(ind,:)
            end
            % Get the biomass (or whatever the specified variable is)
            siteValues(site,week) = str2num(eval(['stockingDataChunk.' variable '{ind}']));
        end
    end
end

