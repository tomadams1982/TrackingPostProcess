% Load the model mesh
mesh=load('C:\Users\sa01ta\Documents\Sealice_NorthMinch\hydro_mesh_run\mesh.mat');

% Load the file containing site locations
dsLocations = datastore('C:\Users\SA01TA\Documents\Sealice_NorthMinch\Kames2018\all_sites_header.dat',...
    'Type','tabulartext','ReadVariableNames',1);
dsLocations.Delimiter = '\t';
dsLocations.NumHeaderLines = 0;
dsLocations.TextscanFormats = {'%s','%d','%d','%f','%f','%f','%f','%s','%s','%s'};
startlocs = readall(dsLocations);

%sitefile='C:\Users\SA01TA\Documents\Sealice_NorthMinch\site_locations\170522_OS_fishfarm3\startlocations_ff3.dat';
%startlocs=load(sitefile);

% Add the path required for plotting
addpath('C:\Users\SA01TA\Documents\code\matlab\trackingProcess');