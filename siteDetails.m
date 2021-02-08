function [locFMA,locDMA,locFHMRA] = siteDetails(x,y,basedir,write)
% SITEDETAILS - get the relevant site information for the set of sites to
% be used in a particle tracking simulation
% 
% Inputs    - locations
%
% Outputs   - siteLocations [numericID,x,y,SEPA_ID,MSS_ID,FMA,DMA,FHMRA]
%           - site biomasses [numericID,biomass(dates)]
%           - site lice counts [numericID,lice/fish(dates)]
    
    % read management area shapefiles
    shapedir = [basedir '\Sealice_NorthMinch\Management_areas\'];
    FMA = shaperead([shapedir '160407_units\FMA_OS.shp']);
    DMA = shaperead([shapedir '160407_units\DMA_OS.shp']);
    FHMRA = shaperead([shapedir '160407_units\FHMRA_OS.shp']);
    
    if isempty(x)
        loc=readtable('C:\Users\sa01ta\OneDrive - SAMS\Documents\Sealice_NorthMinch\SiteLice\siteCounts2018.csv');
        x=loc.Easting;
        y=loc.Northing;
    end
    
    % read site locations and filter by date
%     sitedir = 'C:\Users\SA01TA\Documents\Sealice_NorthMinch\Management_areas\Mapping\';
%     sitefile = [sitedir 'se_monthly_reports.csv'];
%     
%     dates={'01-Sep-17','01-Aug-17'};
%     monthlyData = read_site_MSS_AquaData(sitefile,dates);
    
    locFMA=cell(size(x,1),1);
    locDMA=cell(size(x,1),1);
    locFHMRA=cell(size(x,1),1);

    % check which areas sites belong to
    % FMA
    for i=1:size(FMA,1)
        [in,~]=inpolygon(x,y,FMA(i).X,FMA(i).Y);
        i;
        find(in);
        if ~isempty(find(in))
            locFMA(in)={FMA(i).NAME};
        end
    end
    locFMA(cellfun('isempty',locFMA))={'0'};
    % DMA
    for i=1:size(DMA,1)
        [in,~]=inpolygon(x,y,DMA(i).X,DMA(i).Y);
        i;
        find(in);
        if ~isempty(find(in))
            locDMA(in)={DMA(i).areacode};
        end
    end
    locDMA(cellfun('isempty',locDMA))={'0'};
    % FHMRA
    for i=1:size(FHMRA,1)
        [in,~]=inpolygon(x,y,FHMRA(i).X,FHMRA(i).Y);
        i;
        find(in);
        if ~isempty(find(in))
            locFHMRA(in)={FHMRA(i).AreaID};
        end
    end
    locFHMRA(cellfun('isempty',locFHMRA))={'0'};

    if write==1
        %dlmwrite('siteAreaIDs.csv',[cell2mat(locFMA) cell2mat(locDMA) cell2mat(locFHMRA)]);
        %writecell([locFMA locDMA locFHMRA],'siteAreaIDs.csv');
        writetable(cell2table([locFMA locDMA locFHMRA]),'siteAreaIDs.csv');
        %dat=[locFMA locDMA locFHMRA];
    end

end