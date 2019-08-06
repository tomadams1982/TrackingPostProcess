%% Load the model mesh
load('mesh.mat')
load('modelOutputs.mat')

% Plot model mesh in OS coordinates
clf
plotMeshPDens(mesh,'os',1);

% Plot model mesh in lon/lat
clf
plotMeshPDens(mesh,'os',1);

% Plot model mesh in lon/lat, specifying axis limits
clf
plotMeshPDens(mesh,'os',0,'xl',[-5.5,-5],'yl',[57.8,58.15]);

% Plot a predicted density on the model mesh, with a few options, then add
% some site locations
clf
plotMeshPDens(mesh,'meshDensity',dens,'xl',[-5.5,-5],'yl',[57.8,58.15],'os',0,'logScale',1,'colorBar',1,'zScale',[-5 1]);
hold on
scatter(startlocs.Longitude,startlocs.Latitude,20,'ob','filled')

% Plot connectivity matrix over the model mesh
clf
plotMeshPDens(mesh.mesh,'xl',[-5.8,-5.1],'yl',[56.3,56.9],'os',os);
hold on
plotConnect(mean(connectivity,3),table2array(startlocs(:,2:3)),...
    'logScale',1,'os',os,...
    'colorBar',1);