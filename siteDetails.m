function siteDetails = siteDetails()
% SITEDETAILS - get the relevant site information for the set of sites to
% be used in a particle tracking simulation
% 
% Inputs    - Month
%           - Year
%
% Outputs   - siteLocations [numericID,x,y,SEPA_ID,MSS_ID,FMA,DMA,FHMRA]
%           - site biomasses [numericID,biomass(dates)]
%           - site lice counts [numericID,lice/fish(dates)]
    
    % read management area shapefiles
    
    shapedir = 'C:\Users\SA01TA\Documents\Sealice_NorthMinch\Management_areas\';
    FMA = shaperead([shapedir '160407_units\FMA_OS.shp']);
    DMA = shaperead([shapedir '160407_units\DMA_OS.shp']);
    FHMRA = shaperead([shapedir '160407_units\FHMRA_OS.shp']);
    
    % read site locations and filter by date
    sitedir = 'C:\Users\SA01TA\Documents\Sealice_NorthMinch\Management_areas\Mapping\';
    sitefile = [sitedir 'se_monthly_reports.csv'];
    
    dates={'01-Sep-17','01-Aug-17'};
    monthlyData = read_site_MSS_AquaData(sitefile,dates);
    

    
    
    % check which sites are in the domain
    
    
    % check which areas sites belong to
    
    

end