function [weekStartDates] = createWeekStartDates(startDateString,nWeeks,varargin)
%WEEKSTARTDATES Make an array of start dates to use for aggregating arrays
%in a later calculation.
%   
    narginchk(2,4)

    nDaysPerWeek=7;

    for i = 1:2:length(varargin) % only bother with odd arguments, i.e. the labels
        switch varargin{i}
            case 'nDays'
                nDaysPerWeek = varargin{i+1};
        end
    end
    
    validateattributes(startDateString, {'char'}, ...
        {'scalartext'}, ...
        mfilename, 'the startDateString', 1)
    validateattributes(nWeeks, {'double'}, ...
        {'scalar'}, ...
        mfilename, 'nWeeks', 2)
    validateattributes(nDaysPerWeek, {'double'}, ...
        {'scalar'}, ...
        mfilename, 'nDaysPerWeek', 3)

    startDate=datetime(startDateString,'InputFormat','yyyyMMdd');

    weekStartDates = NaT(1,nWeeks);
    for t=1:nWeeks
        weekStartDates(t)=startDate+(t-1)*nDaysPerWeek;
    end

end

