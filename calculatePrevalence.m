function [prevalence] = calculatePrevalence(density,varargin)

    threshold = 0;
    nWeeks = 52;
    for i = 1:2:length(varargin) % only bother with odd arguments, i.e. the labels
        switch varargin{i}
            case 'threshold'
                threshold = varargin{i+1};
            case 'nWeeks'
                threshold = varargin{i+1};
        end
    end
    threshold
    nWeeks
    prevalence=sum(density>threshold,2)/nWeeks;
end