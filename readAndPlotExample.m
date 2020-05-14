basedir = 'C:\Users\sa01ta\OneDrive - SAMS\Documents\'

% Folder where I have readMeshNC
addpath([basedir 'COMPASS\scripts'])
% Folder where I have plotMeshPDens
addpath(genpath([basedir 'code\matlab\trackingProcess']));

% Standard reading
mesh = readMeshNC([basedir 'particle_track\WestCOMS2_Mesh.nc'])

mesh = readMeshNC([basedir 'particle_track\WestCOMS_mesh.nc'],'readOS',1)


% Also plot the mesh
mesh = readMeshNC([basedir 'particle_track\WestCOMS2_Mesh.nc'],'plotMesh',1)

% Also read Ordnance Survey coordinates
mesh = readMeshNC([basedir 'particle_track\WestCOMS2_Mesh.nc'],'readOS',1)

% If you have a mesh with only latitude/longitude (degrees) OR Ordnanace Survey
% coordinates (metres), you can create the extra fields manually using the SEPA OS
% toolbox. It automatically detects what format the inputs are in.
% Download from: https://github.com/OceanMetSEPA/os_toolbox
addpath(genpath([basedir '\code\matlab\os_toolbox-master']));
% If nodexy and uvnode are in OS
[mesh.nodexy_deg(:,1),mesh.nodexy_deg(:,2)]=OS.convertAndTransform(mesh.nodexy(:,1),mesh.nodexy(:,2));
[mesh.uvnode_deg(:,1),mesh.uvnode_deg(:,2)]=OS.convertAndTransform(mesh.uvnode(:,1),mesh.uvnode(:,2));
% If nodexy and uvnode are in long/lat
[mesh.nodexy_os(:,1),mesh.nodexy_os(:,2)]=OS.convertAndTransform(mesh.nodexy(:,1),mesh.nodexy(:,2));
[mesh.uvnode_os(:,1),mesh.uvnode_os(:,2)]=OS.convertAndTransform(mesh.uvnode(:,1),mesh.uvnode(:,2));
         
%%

% Plot the mesh using the field nodexy, uvnode
plotMeshPDens(mesh,'os',2)
% Plot the mesh using the field nodexy_os, uvnode_os (may not exist)
plotMeshPDens(mesh,'os',1)
% Plot the mesh using the field nodexy_deg, uvnode_deg (may not exist)
plotMeshPDens(mesh,'os',0)



% Plot the mesh using the field nodexy_os, uvnode_os, and apply spatial
% limits
plotMeshPDens(mesh,'os',1,'xl',[100000,200000],'yl',[650000,720000])


%%
hydrofile = 'C:\Users\sa01ta\OneDrive - SAMS\hydroOut\minch_FVCOM\netcdf_2019\minch2_20190101_0007.nc';

u = ncread(hydrofile,'u');
s = ncread(hydrofile,'salinity');
t = ncread(hydrofile,'temp');

%%
% This can be used for u,v with appropriate adjustments to arguments
plotMeshPDens(mesh,'os',2,'meshDensity',abs(u(:,1,1)),'areaScale',0,'FaceAlpha',1,'colorBar',1);

% To plot T, S, need to do something else
doc patch

%%

% specify a point
xy = [-5.598 56.43];
% find nearest location
%addpath('C:\Users\sa01ta\OneDrive - SAMS\Documents\code\matlab\particletracking')
distAll=pdist2(xy,mesh.uvnode);
[dist,nearestID]=min(distAll);

% verify point is correct
hold on
scatter(mesh.uvnode(nearestID,1),mesh.uvnode(nearestID,2))

% plot a timeseries
figure
plot(squeeze(u(nearestID,1,:)))


%% Plot new psteps output
outputDir='C:\Users\sa01ta\OneDrive - SAMS\Documents\OFF-AQUA\liceTracking\pstepsTest\20200430_1\';
pstepFile = [outputDir 'pstepsMature_20190106_143.dat'];
outputDir='C:\Users\sa01ta\OneDrive - SAMS\Documents\OFF-AQUA\liceTracking\pstepsTest\20200430_2\';
pstepFile = [outputDir 'pstepsMature_20190106_130.dat'];


psteps = load(pstepFile);

meshDens = zeros(size(mesh.uvnode,1),1);
meshDens(psteps(:,1)+1) = psteps(:,2);
%meshDens(psteps(:,1)) = sum(psteps(:,2:3),2);
plotMeshPDens(mesh,'meshDensity',meshDens);



plotMeshPDens(mesh,'meshDensitySparse',psteps);
plotMeshPDens(mesh,'meshDensitySparse',[psteps(:,1) sum(psteps(:,2:3),2)]);






