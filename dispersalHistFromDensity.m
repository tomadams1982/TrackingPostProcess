function [distHist,dists,percentileDist] = dispersalHistFromDensity(mesh,density,sourceLocation,varargin)
% DISPERSALHISTFROMDENSITY - Create a histogram of particle transport
% distances (as a proportion of total density) from the source location and
% the particle density distribution.
%
% Calculation is based upon the centroid distance of elements of the
% density array from the source location
%
% Inputs:   mesh        - Model mesh defining the spatial pattern of
%                           density elements
%           density     - Matrix of densities
%           sourceLocation - [x,y] coordinate (OSGB1936 metres) of the source location
%           varargin:
%               histInterval - Distance interval to use for histogram (metres)
%               maxDist      - The maximum distance to record, otherwise add to
%                           last bin
%               percentile   - The percentile of the density at which to
%                           provide a distance (calculated by cumulative
%                           sum)
%
% Outputs:  distHist
%           dists
%           percentileDist
%
% Author: Thomas Adams, SAMS, 2019

    histInterval = 5000;
    maxDist = 80000;
    percentile = 95;

    for i = 1:2:length(varargin) % only bother with odd arguments, i.e. the labels
        switch varargin{i}
            case 'histInterval'
                histInterval = varargin{i+1};
            case 'maxDist'
                maxDist = varargin{i+1};
            case 'percentile'
                percentile = varargin{i+1};
        end
    end

    % Setup output matrix
    distHist=zeros(maxDist/histInterval,2);
    distHist(:,1)=(0:(size(distHist,1)-1))*histInterval;

    % Get the indices of the non-zero elements of the density array, and
    % put the values in an array
    ind = find(density);
    values = density(ind);
    % Get the locations of those element centroids in the mesh
    xy = [mesh.uvnode_os(ind,1),mesh.uvnode_os(ind,2)];

    % Get the distances from source location
    dists=pdist2([sourceLocation(1) sourceLocation(2)],...
        [xy(:,1) xy(:,2)]);

    % Place a value in the relevant histogram bin, scaled by
    % density
    bin=floor(dists/histInterval+1);
    for j=1:length(bin)
        if bin(j)<size(distHist,1)
            distHist(bin(j),2)=distHist(bin(j),2)+values(j);
        else
            distHist(end,2)=distHist(end,2)+values(j);
        end
    end
    
    % Sort the values to get the percentile distance
    [d2,ind2]=sort(dists);
    val2=values(ind2);
    csum=cumsum(val2);
    percentileInd=find(csum<(percentile/100)*sum(values),1,'last');
    percentileDist=d2(percentileInd);
 
    size(dists)
    
    %disp(['histSum = ' num2str(sum(distHist(:,2))) ' ' num2str(sum(values))]);
    disp(['Percentile dist = ' num2str(percentileDist)]);
    
end