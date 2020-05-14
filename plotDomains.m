
%basedir='C:\Users\SA01TA\Documents';
basedir='C:\Users\SA01TA\OneDrive - SAMS\Documents';
addpath([basedir '\code\matlab\trackingProcess']);

% Load and plot FVCOM mesh area and convex hull
mesh=load([basedir '\Sealice_NorthMinch\hydro_mesh_run\mesh.mat']);
mesh=mesh.mesh;

% Need to comment "figure" command in the below script
%hold on
plotMeshPDens(mesh,'os',0)
%plotMeshPDensSimple(mesh,[],[],0,0,'noprint');
%%
hold on

%romsdir='C:\Users\SA01TA\';
romsdir='C:\Users\sa01ta\OneDrive - SAMS\hydroOut\NEA_ROMS\2019\'
% Load and plot NEA ROMS subdomain and convex hull
% http://milas.marine.ie/thredds/catalog/IMI_ROMS_HYDRO/NEATLANTIC_NATIVE_2KM_40L_1H/ANALYSIS/catalog.html
%file='F:\hydroOut\NEA_ROMS\2017_OLD\NEATL_SAMS_2017123000.nc';
file=[romsdir '\NEATL_2019010223.nc'];
%ncdisp(file)
x=ncread(file,'lon_u');
y=ncread(file,'lat_u');
x2=ncread(file,'lon_v');
y2=ncread(file,'lat_v');

u=ncread(file,'ubar');
dom=u;
dom(~isnan(dom))=1;
load([basedir '\particle_track\convexHull_ROMS.dat'])
hold on
p=pcolor(x,y,dom);
cmap=flipud(summer(9));
colormap(cmap(8:9,:))
p.EdgeColor = 'none';
set(p,'facealpha',0.3)
hold on
plot([convexHull_ROMS(:,1);convexHull_ROMS(1,1)],[convexHull_ROMS(:,2);convexHull_ROMS(1,2)])

hold on
plotMeshPDens(mesh,'os',0)
%plotMeshPDensSimple(mesh,[],[],0,0,'noprint');

% Plot FVCOM convex hull
load([basedir '\particle_track\convexHull_FVCOM.dat'])

plot([convexHull_FVCOM(:,1);convexHull_FVCOM(1,1)],[convexHull_FVCOM(:,2);convexHull_FVCOM(1,2)])
plot([convexHull_FVCOM(:,1)],[convexHull_FVCOM(:,2)])

%%
xlim([-8.5,-4])
ylim([53,60])

%% Save figure
print('-painters','-dpng','-r600','domains_FVCOM_ROMS_close.png')

%%

bn=ncread('C:\Users\SA01TA\Documents\particle_track\WestCOMS_mesh.nc','boundaryNodesAll');
nodexy=ncread('C:\Users\SA01TA\Documents\particle_track\WestCOMS_mesh.nc','nodexy');