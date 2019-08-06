function scaling = fhmraAverageCountToSites(sites,FHMRAdata,date1,date2)
% FHMRAAVERAGECOUNTTOSITES - Create scaling matrix from a table of startlocations, list of dates, and
% table of lice count values for FHMRAs
% 
% Inputs:   sites - a table containing sites, with ID, x,y location and
%                   management area assignation
%           FHMRAdata - matrix of lice counts
%           date1 - date on which to get data. If empty will get mean 
%           date2 - if supplied, mean value will be calculated over range
%           [date1,date2]
% Outputs:  scaling - lice count values pertaining to each site
    
    narginchk(2,4)
    if (nargin < 4)
        warning('fhmraCountScaling:dateWarning2','Less than 4 arguments, getting value for specific date');
        if (nargin < 3)
            warning('fhmraCountScaling:dateWarning1','Less than 3 arguments, getting mean value for each area');
            date1=[];
        end
        date2=[];
    end
          
    % Set the rows to extract counts from the table based on supplied dates
    if ~isempty(date2)
        % Set to first of month since this is what at "dayless" date like
        % FHMRAdata.Date defaults to
        date1.Day=1;
        date2.Day=1;         
        row = find(FHMRAdata.Date>=date1 & FHMRAdata.Date<=date2);
    elseif ~isempty(date1)
        date1.Day=1;
        row = find(FHMRAdata.Date==date1);
    else
        row = (1:size(FHMRAdata.Date))';
    end
    row;
    
    % Select the column of the FHMRA counts that applies to each site in
    % the list, and average count values over the given date period
    scaling=zeros(size(sites,1),1);
    for i=1:size(sites,1)
        col=find(ismember(FHMRAdata.Properties.VariableNames,['x' sites.FHMRA{i}]));
        
        if ~isempty(FHMRAdata(row,col))
            %FHMRAdata(row,col)
            count=mean(table2array(FHMRAdata(row,col)));
            scaling(i)=count;
        else
            count=0;
            scaling(i)=0;
        end
        
        %fprintf('i %d: col %d startlocs.FHMRA %s count %d\n',i,col,sites.FHMRA{i},count);
    end
end