%% testOutputScenario

% set to cover right hand monitor
xl=[0,0];yl=[0,0];
% Loch Pooltiel limits 1
%xl=[60000,180000];yl=[800000,900000];
% Loch Pooltiel limits 2 (closer - connectivity)
%xl=[90000,130000];yl=[810000,890000];
% Test1 limits
%xl=[100000,220000];yl=[660000,760000];

% Eddrachillis
xl=[200000,230000];yl=[925000,960000];

%% Initial setup
% - add paths
% - Go to directory for computation
% - Set start date
% - set number of weeks
clear all
addpath('C:\Users\SA01TA\Documents\code\matlab\trackingProcess');
addpath('C:\Users\SA01TA\Documents\code\matlab\cmocean');

%cd W:\sa01ta\fishfarm_long\20160301_135_160301_134
cd W:\sa01ta\fishfarm_long\20160301_135_160301_134_maxage336
%cd W:\sa01ta\fishfarm_long\20160614_135_160614_134_maxage336
%cd W:\sa01ta\fishfarm_long\20160614_135_160614_134_maxage336_dev70h

startDateString='20160301'; % start date for W:\sa01ta\fishfarm_long\20160301_135_160301_134
startDate=datetime(startDateString,'InputFormat','yyyyMMdd');
nweeks=18;
nSites=160;

% Global options
getFHMRAcounts=1; % Should FHMRA monthly counts be used, or treatment limit?
oldnumeric=0; % Does the site locations file contain SiteIDs in the old numeric (0:nsites-1) or SEPA ID?

% Calculate the particle count scaling
%   = treatment limit * Nfishpertonne / (Nmodellicereleasedperhour * Mdepthdistributed)
% In count scaling, the number is also divided by the number of hours
% covered by the requested time window, and multiplied by biomass. 
% Fish per tonne is based upon the mean over the whole production cycle from 
% Kames stocking plan (04/2018), though varies from >8000 early in cycle to <1000
% later on

% Set values for scaling the calculated densities and connectivities
nhours=168; % Number of hours in a week, to divide the summed counts by.
nFishPerTonne=2300;
liceDepthThickness=5;
nPartsPerSite=5;
fecundityDaily=28.2; % L447
hoursPerDay=24; % L447
treatLimit=0.5;

mult = ones(nSites,1)*treatLimit*nFishPerTonne/(liceDepthThickness*nPartsPerSite);
mult_newSite= treatLimit*nFishPerTonne/(liceDepthThickness*nPartsPerSite);

% Site positions in table
existingSiteRows=1:156;
newSiteRow=[157];
%newSiteRow=[157:160]; % position of new site in table
newSiteMeanBio=1006.5;
newSiteMaxBio=2575.2;

%% Coord ref system
os=1;


%% load the site location data
%loadMeshAndSites
% Load the model mesh
mesh=load('C:\Users\sa01ta\Documents\Sealice_NorthMinch\hydro_mesh_run\mesh.mat');

if oldnumeric==0
    % Load the file containing site locations
    dsLocations = datastore('C:\Users\SA01TA\Documents\Sealice_NorthMinch\Kames2018\all_sites_header_WGS84.dat',...
        'Type','tabulartext','ReadVariableNames',1);
    dsLocations.Delimiter = '\t';
    dsLocations.NumHeaderLines = 0;
    dsLocations.TextscanFormats = {'%s','%d','%d','%f','%f','%f','%f','%s','%s','%s'};
    startlocs = readall(dsLocations);
else
    % OR old numeric
    startlocs=load('W:\sa01ta\fishfarm_long\20130719_120_130719_120\startlocs.dat');
end

% %%%%%%% Setup some new site stuff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

startlocs.meanBio(newSiteRow)=newSiteMeanBio;
startlocs.maxBio(newSiteRow)=newSiteMaxBio;

%% Read in FHMRA data if required
if getFHMRAcounts==1
    fhmraCountFile='C:\Users\SA01TA\Documents\Sealice_NorthMinch\SSPO_reports\180423_countsTranspose.csv';
    dsCounts = datastore(fhmraCountFile,'Type','tabulartext','ReadVariableNames',1);
    dsCounts.Delimiter = ',';
    dsCounts.NumHeaderLines = 0;
    dsCounts.TextscanFormats = {'%{MMM-yy}D','%f','%f','%f','%f','%f','%f','%f','%f','%f','%f',...
        '%f','%f','%f','%f','%f','%f','%f','%f','%f','%f',...
        '%f','%f','%f','%f','%f','%f','%f','%f','%f','%f'};
    FHMRAdata = readall(dsCounts);
end

%% Load in site lice data
dsSiteLice = datastore('C:\Users\SA01TA\Documents\Sealice_NorthMinch\SiteLice\Counts_wCodes.csv',...
        'Type','tabulartext','ReadVariableNames',1);
dsSiteLice.Delimiter = ',';
dsSiteLice.NumHeaderLines = 0;
dsSiteLice.TextscanFormats = {'%s','%s','%s','%s','%f','%s'};
siteLice = readall(dsSiteLice);

% Remove the NA entries
siteLice=siteLice(~strcmp(siteLice.SEPA_Site,'NA'),:);

% Check which sites are duplicated
%[~,idxu,idxc] = unique(siteLice(:,2));
% count unique values (use histc in <=R2014b)
%[count, ~, idxcount] = histcounts(idxc,numel(idxu));
% Where is greater than one occurence
%idxkeep = count(idxcount)>1;
% Extract from C
%siteLice(idxkeep,:)
% Remove the duplicate entries - probably ARDT1 (58) should be removed
% anyway as ID was guessed....
siteLice([54,160,177],:)=[];

% Add a column containing site counts to the startlocs table
[idxA, idxB] = ismember(startlocs.SEPA_Site, siteLice.SEPA_Site);
startlocs.siteCount=NaN(size(startlocs,1),1);
startlocs(idxA,'siteCount') = siteLice(idxB(idxA),'Count');

% Make an FHMRA aggregated mean column based on site counts
% Doesn't work as non-numeric
%accumarray(startlocs.FHMRA,startlocs.siteCount,[size(startlocs,1) 1],@sum)
% NOTE "grpstats" treats NaNs as non-existent. Is this what you really
% want? I think it is...
fhmraSumCounts = grpstats(startlocs,'FHMRA','mean','DataVars','siteCount');
[idxA, idxB] = ismember(startlocs.FHMRA, fhmraSumCounts.FHMRA);
startlocs.fhmraSumCount=NaN(size(startlocs,1),1);
startlocs(idxA,'fhmraSumCount') = fhmraSumCounts(idxB(idxA),'mean_siteCount');

% Set the new site counts to zero (end up being non-zero due to matching an
% NA in row 6).
startlocs.fhmraSumCount(157:160)=NaN;



%% set site ID lists for spliting up particle locations and making element
% counts
if oldnumeric==0
    allsiteIDs=startlocs.SEPA_Site;
    existingsiteIDs=startlocs.SEPA_Site(existingSiteRows);
    newsiteID=startlocs.SEPA_Site(newSiteRow);
else
    % OR old numeric
    %allsiteIDs=cellfun(@num2str, num2cell(startlocs(:,1)), 'UniformOutput', false);
    allsiteIDs=cellfun(@num2str, num2cell([0:size(startlocs,1)-1]'), 'UniformOutput', false);
end

%% plot
%locXY=double(startlocs{:,2:3});
%locXY=startlocs(:,2:3);
plotMeshPDensSimple(mesh.mesh,0,[startlocs.Easting startlocs.Northing],os,'noprint');
hold on
for i=1:size(startlocs,1)
    text(double(startlocs.Easting(i)), double(startlocs.Northing(i)), allsiteIDs{i});
end

%% %%%%%%% Calculate element counts and element source counts %%%%%%%%%%%%%%%%%%%%%%%%

load elementCounts160507_3wks
% Just keep  countsExisting countsExisting_juv
clear countsSingle countsSingle_juv countsAll countsAll_juv sourceCountsAll sourceCountsAll_juv 

fprintf('***** Element counts *****\n');

%countsSingle=zeros(size(mesh.mesh.uvnode,1),nweeks);
%countsSingle_juv=zeros(size(mesh.mesh.uvnode,1),nweeks);
countsSingle_test1=zeros(size(mesh.mesh.uvnode,1),nweeks);
countsSingle_juv_test1=zeros(size(mesh.mesh.uvnode,1),nweeks);
countsSingle_test2=zeros(size(mesh.mesh.uvnode,1),nweeks);
countsSingle_juv_test2=zeros(size(mesh.mesh.uvnode,1),nweeks);

%countsAll=zeros(size(mesh.mesh.uvnode,1),nweeks);
%countsAll_juv=zeros(size(mesh.mesh.uvnode,1),nweeks);

sourceCountsAll=zeros(size(mesh.mesh.uvnode,1),nweeks);
sourceCountsAll_juv=zeros(size(mesh.mesh.uvnode,1),nweeks);

%%
countsExisting_siteScale=zeros(size(mesh.mesh.uvnode,1),nweeks);
countsExisting_fhmraScale=zeros(size(mesh.mesh.uvnode,1),nweeks);

%%

% Start here if fails partway
tic
for t=7:7
    weekStartDate=startDate+(t-1)*7;
    weekStartDate;
    %endDate=startDate;
    endDate=weekStartDate+6;
    fprintf('Week %d - %s to %s (inclusive)\n',t,datestr(weekStartDate,'yyyymmdd'),datestr(endDate,'yyyymmdd'));

    % Multiply by FHMRA counts in that month
    if getFHMRAcounts==1
        % Apply FHMRA counts to the sites
        siteCounts = fhmraAverageCountToSites(startlocs,FHMRAdata,weekStartDate,[]);  
        mult = siteCounts*nFishPerTonne/(liceDepthThickness*nPartsPerSite);
        mult(newSiteRow)=mult_newSite;
    end
    
    % ---------------- specific existing site count scenarios ------------
%     status=2;
%     locationsToCount=existingsiteIDs;
%     
%     mult = startlocs.siteCount(existingSiteRows)*nFishPerTonne/(liceDepthThickness*nPartsPerSite);
%     scaling=startlocs.maxBio(existingSiteRows).*mult/nhours;
%     [c3,sc3]=elementCounts(mesh,weekStartDate,endDate,0,24,status,locationsToCount,scaling);
%     countsExisting_siteScale(:,t)=c3(:,3);
%     
%     mult = startlocs.fhmraSumCount(existingSiteRows)*nFishPerTonne/(liceDepthThickness*nPartsPerSite);
%     scaling=startlocs.maxBio(existingSiteRows).*mult/nhours;
%     [c3,sc3]=elementCounts(mesh,weekStartDate,endDate,0,24,status,locationsToCount,scaling);
%     countsExisting_fhmraScale(:,t)=c3(:,3);
    
    
    
    
    
    % -------------------- single site -------------------------
%     startlocs([newSiteRow],:);
%     locationsToCount=allsiteIDs([newSiteRow]);
%     %locationsToCount=allsiteIDs([1,2]);
%     scaling=table2array(startlocs([newSiteRow],6)).*mult(newSiteRow)/nhours;

%     status=1;
%     [c1,sc1]=elementCounts(mesh,weekStartDate,endDate,0,24,status,locationsToCount,scaling);
%     countsSingle_juv(:,t)=c1(:,3); 
%     status=2;
%     [c2,sc2]=elementCounts(mesh,weekStartDate,endDate,0,24,status,locationsToCount,scaling);
%     countsSingle(:,t)=c2(:,3); 
%     
%     clear c1 c2;
    
    scal=scaling(3);
    locToCnt=locationsToCount(3);
    status=1;
    [c1,sc1]=elementCounts(mesh,weekStartDate,endDate,0,24,status,locToCnt,scal);
    countsSingle_juv_test1(:,t)=c1(:,3); 
    status=2;
    [c2,sc2]=elementCounts(mesh,weekStartDate,endDate,0,24,status,locToCnt,scal);
    countsSingle_test1(:,t)=c2(:,3); 
    
    clear c1 c2;

    scal=scaling(4);
    locToCnt=locationsToCount(4);
    status=1;
    [c1,sc1]=elementCounts(mesh,weekStartDate,endDate,0,24,status,locToCnt,scal);
    countsSingle_juv_test2(:,t)=c1(:,3); 
    status=2;
    [c2,sc2]=elementCounts(mesh,weekStartDate,endDate,0,24,status,locToCnt,scal);
    countsSingle_test2(:,t)=c2(:,3); 

    clear c1 c2;
    
    % -------------------- all sites ---------------------------
%     startlocs;
%     locationsToCount=allsiteIDs;
%     % remove the sites to be ignored
%     locationsToCount(setdiff(1:size(startlocs,1),[existingSiteRows newSiteRow]))=[];
%     scaling=table2array(startlocs([existingSiteRows newSiteRow],6)).*mult([existingSiteRows newSiteRow])/nhours;
%     %scaling=startlocs(:,4)*mult/nhours;
% 
%     status=1;
%     [c3,sc3]=elementCounts(mesh,weekStartDate,endDate,0,24,status,locationsToCount,scaling);
%     countsAll_juv(:,t)=c3(:,3);
%     sourceCountsAll_juv=sc3(:,2);
%     status=2;
%     [c4,sc4]=elementCounts(mesh,weekStartDate,endDate,0,24,status,locationsToCount,scaling);
%     countsAll(:,t)=c4(:,3);
%     sourceCountsAll=sc4(:,2);  
%     clear c3 c4;
    
    % Advance one week
    weekStartDate=weekStartDate+7;
end

% -------------------- existing sites ---------------------------
% countsExisting=countsAll-countsSingle;
% countsExisting_juv=countsAll_juv-countsSingle_juv;
% 
% %
% save elementCounts160507_3wks countsSingle countsSingle_juv countsAll countsAll_juv countsExisting countsExisting_juv sourceCountsAll sourceCountsAll_juv 

% countsAll_test1=countsExisting+countsSingle_test1;
% countsAll_juv_test1=countsExisting_juv+countsSingle_juv_test1;
% 
% countsAll_test2=countsExisting+countsSingle_test2;
% countsAll_juv_test2=countsExisting_juv+countsSingle_juv_test2;

toc

save countsExisting_varyScale countsExisting_siteScale countsExisting_fhmraScale

%save work180530


%%
% ---------------------Connectivity calculation -------------------------------
%directory='W:\sa01ta\fishfarm_long\20160301_135_160301_134_maxage336\'
%startDateString='20160301'; % start date for W:\sa01ta\fishfarm_long\20160301_135_160301_134
%startDate=datetime(startDateString,'InputFormat','yyyyMMdd');
%nweeks=18;

directory='W:\sa01ta\fishfarm_long\20160614_135_160614_134_maxage336\'
startDateString='20160614'; % start date for W:\sa01ta\fishfarm_long\20160301_135_160301_134
startDate=datetime(startDateString,'InputFormat','yyyyMMdd');
%nweeks=15;

fprintf('***** Connectivity *****\n');
% Creating a datastore fails when there is no data in the columns....
%startDateString='20130802';

nhours=168;
nparts=5; % number of particles per site per hour

connectivity=zeros(size(allsiteIDs,1),size(allsiteIDs,1),nweeks);

%weekStartDate=startDate;
weekStartDate=startDate+1*7;
tic
for t=2:nweeks
    weekStartDate;
    %endDate=startDate;
    endDate=weekStartDate+6;
    fprintf('Week %d - %s to %s (inclusive)\n',t,datestr(weekStartDate,'yyyymmdd'),datestr(endDate,'yyyymmdd'));
    
    %siteIDs=allsiteIDs;
    
    % Scaling is not currently used in connectivity calculation - can be
    % applied afterwards in any case
    %scaling=startlocs.maxBio*mult/nhours;
    %scaling=startlocs(:,4)*mult/nhours;
    
    c1 = connectivityFromArrivals(directory,allsiteIDs,weekStartDate,endDate,nparts,[]);
    
    connectivity(:,:,t)=c1;
    
    clear c1;
    
    % Advance one week
    weekStartDate=weekStartDate+7;
end
toc

% Need to divide by the number of hours represented by each file, as each
% file represents an opportunity for NpartsPerHour X Nhours to successfully
% go from each site A to site B
connectivity = connectivity/nhours;

%
save connectivity_redone connectivity


%%

load('C:\Users\sa01ta\Documents\Sealice_NorthMinch\180516_siteScoping\connectivity_160301_160704.mat')
connectivity1 = connectivity(:,:,2:end);
load('C:\Users\sa01ta\Documents\Sealice_NorthMinch\180516_siteScoping\connectivity_160621_161017.mat')
connectivity2 = connectivity(:,:,4:end);
connectivity = cat(3,connectivity1,connectivity2);




%% ---------------- plot the mean density -----------------------------------
%locXY=double(startlocs{:,2:3});

% Read data from two different directories and combine

%cd W:\sa01ta\fishfarm_long\20160301_135_160301_134
cd W:\sa01ta\fishfarm_long\20160301_135_160301_134_maxage336
load elementCounts160507_3wks.mat
countsSingle_juv1=countsSingle_juv(:,3:end);
countsSingle_cop1=countsSingle(:,3:end);
countsExisting_juv1=countsExisting_juv(:,3:end);
countsExisting_cop1=countsExisting(:,3:end);
countsAll_juv1=countsAll_juv(:,3:end);
countsAll_cop1=countsAll(:,3:end);

cd W:\sa01ta\fishfarm_long\20160614_135_160614_134_maxage336
%cd W:\sa01ta\fishfarm_long\20160614_135_160614_134_maxage336_dev70h
load elementCounts160507_3wks.mat
countsSingle_juv2=countsSingle_juv(:,6:end); % files in previous directory go to 13/07/2016
countsSingle_cop2=countsSingle(:,6:end);
countsExisting_juv2=countsExisting_juv(:,6:end);
countsExisting_cop2=countsExisting(:,6:end);
countsAll_juv2=countsAll_juv(:,6:end);
countsAll_cop2=countsAll(:,6:end);

clear countsAll countsAll_juv countsSingle countsSingle_juv

% Combine the desired counts
countsSingle_juv = [countsSingle_juv1 countsSingle_juv2];
countsSingle_cop = [countsSingle_cop1 countsSingle_cop2];

countsExisting_juv = [countsExisting_juv1 countsExisting_juv2];
countsExisting_cop = [countsExisting_cop1 countsExisting_cop2];

countsAll_juv = [countsAll_juv1 countsAll_juv2];
countsAll_cop = [countsAll_cop1 countsAll_cop2];

cd ..

%%

%cd C:\Users\SA01TA\Documents\Sealice_NorthMinch\Council_SEPA_MSS_meetings

% Load files:
% MAIN WORKSPACE
% C:\Users\sa01ta\Documents\Sealice_NorthMinch\180516_siteScoping
% EXTRA SINGLE TEST SITES
% C:\Users\sa01ta\Documents\Sealice_NorthMinch\Council_SEPA_MSS_meetings\work180524_3.mat

focalSite1=159;
focalSite2=160;

% Plot density 
% - spatial domain is hard-coded into plotMeshPDensSimple
% - zlim is hard-coded into plotMeshPDensSimple
%plotMeshPDensSimple(mesh.mesh,[],[],'noprint');
%plotMeshPDensSimple(mesh.mesh,mean(countsExisting_juv,2)*28.2/24,[],'noprint');
%plotMeshPDensSimple(mesh.mesh,mean(countsExisting,2)*28.2/24,[],'noprint');
%plotMeshPDensSimple(mesh.mesh,mean(countsSingle_juv_test2,2)*28.2/24,[],'noprint');
%plotMeshPDensSimple(mesh.mesh,mean(countsSingle_test1,2)*28.2/24,[],'noprint');
%plotMeshPDensSimple(mesh.mesh,mean(countsAll_juv_test2,2)*28.2/24,[],'noprint');

% Scaling the calculated value by 400/2300 to change N fish per tonne to
% 400 (not 2300)

[lon,lat]=OS.convertAndTransform(startlocs.Easting,startlocs.Northing);

% counts = countsExisting_cop;
counts = countsSingle_test1+countsSingle_test2;
% counts = countsExisting_cop+...
%     [countsSingle_test1,repmat(mean(countsSingle_test1,2),1,11)]+...
%     [countsSingle_test2,repmat(mean(countsSingle_test2,2),1,11)];
os=0;
[out,z2]=plotMeshPDensSimple(mesh.mesh,mean(counts,2)*(fecundityDaily/hoursPerDay)*(400/2300),[],os,1,'noprint');

%out=plotMeshPDensSimple(mesh.mesh,mean(countsExisting_fhmraScale,2)*28.2/(24),[],'noprint');


% Plot sites
hold on
if oldnumeric==0
    % New site
    %scatter(startlocs.Easting(focalSite1),startlocs.Northing(focalSite1),30,'+','black')
    %scatter(startlocs.Easting(focalSite2),startlocs.Northing(focalSite2),30,'+','black')
    scatter(lon(focalSite1),lat(focalSite1),30,'+','black')
    scatter(lon(focalSite2),lat(focalSite2),30,'+','black')
    
    % Existing sites - active
    existingSites=startlocs(1:156,:);
    %scatter(existingSites.Easting,existingSites.Northing,20,'filled','blue');
    %scatter(existingSites.Easting(existingSites.currentBio~=0),existingSites.Northing(existingSites.currentBio~=0),10,'filled','black')
    scatter(lon(existingSites.currentBio~=0),lat(existingSites.currentBio~=0),10,'filled','black')
    %scatter(startlocs.Easting,startlocs.Northing,1+startlocs.meanBio/5,'filled')
    % Existing sites - inactive
    %hold on
    %scatter(existingSites.Easting(existingSites.currentBio==0),existingSites.Northing(existingSites.currentBio==0),5,'filled','red')
    scatter(lon(existingSites.currentBio==0),lat(existingSites.currentBio==0),5,'filled','red')
    %scatter(existingSites.Easting(existingSites.meanBio==0),existingSites.Easting(existingSites.meanBio==0),15,'filled','red')
else
    scatter(startlocs(:,2),startlocs(:,3),1+startlocs(:,4)/5,'filled')
    scatter(startlocs(startlocs(:,4)==0,2),startlocs(startlocs(:,4)==0,3),15,'filled','blue')
end
% hold on
% for i=1:size(existingSites,1)
%     text(existingSites.Easting(i),existingSites.Easting(i),existingSites.SEPA_Site(i));
% end
%

% for i=1:size(mesh.mesh.uvnode_os,1)
%     text(mesh.mesh.uvnode_os(i,1),mesh.mesh.uvnode_os(i,2),int2str(i));
% end

%colormap(flipud(pink))

colormap(parula(5))
%c=cmocean('deep');
%colormap(c([1,50,100,150,200,250],:))

cb1=colorbar('Position',[.87 .13 .05 .77]);
cb1.Label.String = 'log_{10}(density) (m^{-3})';
title('(a)')
%%
print('-painters','-dpng','-r600','density_copepodid_new_190607.png')
% print('-painters','-dpng','-r600','density_copepodid_existing_190607.png')
% print('-painters','-dpng','-r600','density_copepodid_all_190607.png')
%close
%%
hist(z2)
h1=histogram('BinCounts',~isnan(z2),...
    'Normalization','probability','DisplayStyle','stairs','Linestyle','-','EdgeColor','k','LineWidth',1);
histogram(z2(z2<1))

z2(isnan(z2))=0;
length(z2(z2<1))/length(z2)
length(z2(z2<0.1))/length(z2)
length(z2(z2<0.01))/length(z2)
length(z2(z2<0.001))/length(z2)

%% --------------- plot weekly densities -----------------------------------
weekDates=startDate+(0:28)*7;
% set x and y limits
outdir='C:\Users\SA01TA\Documents\Sealice_NorthMinch\180516_siteScoping\weeklyFigs\'

% Add a loop to work through weeks and plot separately here (for movie, in
% addition to mean/prevalence plot)
for t=1
    plotMeshPDensSimple(mesh.mesh,countsExisting_cop(:,t)*(fecundityDaily/hoursPerDay)*(400/2300),[],os,1,'noprint');
    hold on
    if oldnumeric==0
        % New site
        scatter(startlocs.Easting(159:160),startlocs.Northing(159:160),30,'+','black')
        % Existing sites - active
        existingSites=startlocs(1:156,:);
        scatter(existingSites.Easting(existingSites.currentBio~=0),existingSites.Northing(existingSites.currentBio~=0),10,'filled','black')
        % Existing sites - inactive
        scatter(existingSites.Easting(existingSites.currentBio==0),existingSites.Northing(existingSites.currentBio==0),5,'filled','red')
    end
    
    cb1=colorbar('Position',[.87 .13 .05 .77]);
    cb1.Label.String = 'log_{10}(Density) (m^{-3})';

    colormap(parula(5))
    title(datestr(weekDates(t)))
    
%     print('-painters','-dpng','-r600',[outdir 'density_copepodid_existing_' int2str(t) '.png'])
%     close
end

%% Plot dispersal kernels for sites
% - Density of copepodid particles at distances from source
% All sites (existing)
tic

startDate=datetime(2016,03,01);
endDate=datetime(2016,07,13);

% These first tow take about an hour to calculate for an initial 19 week
% period
[dispDists159,dispDists159_noscale] = dispersalDistances(mesh,startlocs,startDate,endDate,0,24,2,159,1000,50000);
disp('done 1')
[dispDists160,dispDists160_noscale] = dispersalDistances(mesh,startlocs,startDate,endDate,0,24,2,160,1000,50000);
disp('done 2')

% This one took about half an hour
[dispDistsExisting,dispDistsExisting_noscale] = dispersalDistances(mesh,startlocs,startDate,endDate,0,24,2,1:156,1000,50000);
disp('done existing')

toc

%save dispDists_160301_160713 dispDists159 dispDists160 dispDistsExisting dispDists159_noscale dispDists160_noscale dispDistsExisting_noscale

% Want to plot a histogram/kernel density estimate of all distances excluding the focal site, 
% with histogram of the focal site(s) overlaid
%%
% 
% subplot(1,2,1)
% h1=histogram('BinEdges',0:32,'BinCounts',log10(dispDists159(1:32)+1),'DisplayStyle','stairs');
% hold on
% h2=histogram('BinEdges',0:32,'BinCounts',log10(dispDists160(1:32)+1),'DisplayStyle','stairs');
% h3=histogram('BinEdges',0:32,'BinCounts',log10(dispDistsExisting(1:32)+1),'DisplayStyle','stairs');
% legend('Site 1', 'Site 2', 'Existing sites');

%subplot(1,2,1)
h1=histogram('BinEdges',0:32,'BinCounts',dispDists159(1:32),'Normalization','probability','DisplayStyle','stairs','Linestyle','-','EdgeColor','k','LineWidth',1);
hold on
h2=histogram('BinEdges',0:32,'BinCounts',dispDists160(1:32),'Normalization','probability','DisplayStyle','stairs','Linestyle','--','EdgeColor','k','LineWidth',1);
h3=histogram('BinEdges',0:32,'BinCounts',dispDistsExisting(1:32),'Normalization','probability','DisplayStyle','stairs','Linestyle',':','EdgeColor','k','LineWidth',1);
legend('Site 1', 'Site 2', 'Existing sites')
set(gca,'yscale','log')
xlabel('Distance from source (km)')
ylabel('Proportion of copepodid particles')

% subplot(1,2,2)
% h1=histogram('BinEdges',0:25,'BinCounts',log10(dispDists159_noscale(1:25)+1),'Normalization','probability','DisplayStyle','stairs');
% hold on
% h2=histogram('BinEdges',0:25,'BinCounts',log10(dispDists160_noscale(1:25)+1),'Normalization','probability','DisplayStyle','stairs');
% h3=histogram('BinEdges',0:25,'BinCounts',log10(dispDistsExisting_noscale(1:25)+1),'Normalization','probability','DisplayStyle','stairs');
% legend('Site 1', 'Site 2', 'Existing sites');
% set(gca,'yscale','log')

print('-painters','-dpng','-r600','dispKernels_incMort_160301_160713.png')

    
%% ------------------ prevalence plot ---------------------------------------
% if oldnumeric==0
%     locXY=[startlocs.Easting([1,2]),startlocs.Northing([1,2])];
% else
%     locXY=startlocs([1,2],2:3);
% end
% locationsToCount=allsiteIDs([1,2]);

focalSite1=159;
focalSite2=160;

counts = countsExisting_cop+countsSingle_test2;
%counts = countsSingle_test1;

nSimWeeks=18;
[out,z2]=plotMeshPDensSimple(mesh.mesh,sum(counts~=0,2)/nSimWeeks,[],os,0,'noprint');
hold on
% for i=1:size(locXY,1)
%     text(locXY(i,1),locXY(i,2),locationsToCount{i});
% end

% if oldnumeric==0
%     scatter(startlocs.Easting,startlocs.Northing,1+startlocs.meanBio/5,'filled')
%     scatter(startlocs.Easting(startlocs.meanBio==0),startlocs.Easting(startlocs.meanBio==0),15,'filled','blue')
% else 
%     scatter(startlocs(:,2),startlocs(:,3),1+startlocs(:,4)/5,'filled')
%     scatter(startlocs(startlocs(:,4)==0,2),startlocs(startlocs(:,4)==0,3),15,'filled','blue')
% end

if oldnumeric==0
    % New site
    scatter(startlocs.Easting(focalSite1),startlocs.Northing(focalSite1),30,'+','black')
    scatter(startlocs.Easting(focalSite2),startlocs.Northing(focalSite2),30,'+','black')
    % Existing sites - active
    existingSites=startlocs(1:156,:);
    scatter(existingSites.Easting(existingSites.currentBio~=0),existingSites.Northing(existingSites.currentBio~=0),10,'filled','black')
    % Existing sites - inactive
    scatter(existingSites.Easting(existingSites.currentBio==0),existingSites.Northing(existingSites.currentBio==0),5,'filled','red')
end

%map = [1:-0.2:0; 1:-0.2:0; 1:-0.2:0]';
%map = [1:-0.2:0; 1; 1]';
colormap(parula(5))
caxis([0 1]) % sets colorbar limits 
cb1=colorbar('Position',[.87 .13 .05 .77]);
cb1.Label.String = 'Prevalence (proportion of weeks)';

%print('-painters','-dpng','-r600','prevalence_copepodid_testBoth_190501.png')


%% Plot the number of sources at each cell
% plotMeshPDensSimple(mesh.mesh,sourceCountsAll,[],'noprint');
% colorbar
% colormap(parula(5))



%% ----- Load connectivity -----------------------------
% Read data from two different directories and combine

%cd W:\sa01ta\fishfarm_long\20160301_135_160301_134
%cd W:\sa01ta\fishfarm_long\20160301_135_160301_134_maxage336
load W:\sa01ta\fishfarm_long\20160301_135_160301_134_maxage336\connectivity160507.mat
connectivity1=connectivity(:,:,3:end);

%cd W:\sa01ta\fishfarm_long\20160614_135_160614_134_maxage336
%cd W:\sa01ta\fishfarm_long\20160614_135_160614_134_maxage336_dev70h
load W:\sa01ta\fishfarm_long\20160614_135_160614_134_maxage336\connectivity160507.mat
connectivity2=connectivity(:,:,6:end); % files in previous directory go to 13/07/2016

cd ..

clear connectivity

connectivity=cat(3,connectivity1,connectivity2);

%% -------------- Plot connectivity map ----------------------------
FHMRA=shaperead('C:\Users\SA01TA\Documents\Sealice_NorthMinch\Management_areas\FHMRA.shp');
FMA=shaperead('C:\Users\SA01TA\Documents\Sealice_NorthMinch\Management_areas\FMA.shp');
%%
focalSite=159:160;

plotMeshPDensSimple(mesh.mesh,[],[],1,0,'noprint');

mapshow(FMA,'FaceAlpha',0.10)
for i=1:60%size(FMA) % Stop at 60 to prevent adding Gigha area text over axis
    [m,ind]=max(FMA(i).X(1:end-1));
    text(mean(FMA(i).X(ind)-2000),mean(FMA(i).Y(ind)),FMA(i).NAME)
end

connectivity2=connectivity;
connectivity2(157:158,:)=0;
connectivity2(:,157:158)=0;

%startlocs2=startlocs;
%startlocs2(157:158

hold on
if oldnumeric==0
    plotConnections(mean(connectivity2,3),[startlocs.Easting,startlocs.Northing],0.0001,100);
else
    plotConnections(mean(connectivity2,3),startlocs2(:,2:3),0.0001,100);
end

if oldnumeric==0
    % New site
    scatter(startlocs.Easting(focalSite),startlocs.Northing(focalSite),100,'+','black')
    scatter(startlocs.Easting(focalSite),startlocs.Northing(focalSite),100,'o','black')
    % Existing sites - active
    existingSites=startlocs(1:156,:);
    scatter(existingSites.Easting(existingSites.currentBio~=0),existingSites.Northing(existingSites.currentBio~=0),10,'filled','black')
    % Existing sites - inactive
    scatter(existingSites.Easting(existingSites.currentBio==0),existingSites.Northing(existingSites.currentBio==0),5,'filled','red')
end

% if oldnumeric==0
%     scatter(startlocs.Easting,startlocs.Northing,1+startlocs.meanBio/5,'filled','red')
%     scatter(startlocs.Easting(startlocs.meanBio==0),startlocs.Easting(startlocs.meanBio==0),15,'filled','blue')
% else
%     scatter(startlocs(:,2),startlocs(:,3),1+startlocs(:,4)/5,'filled','red')
%     scatter(startlocs(startlocs(:,4)==0,2),startlocs(startlocs(:,4)==0,3),15,'filled','blue')
% end

%print('-painters','-dpng','-r600','connectivity_190509.png')

%%
x = 0:.05:2*pi;
y = sin(x);
z = zeros(size(x));
col = x;  % This is the color, vary with x in this case.
surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);


%%
% Get order S->N
[B,I]=sort(startlocs.Northing);
% TEST1 is 33, TEST2 is 34
I=I(23:44);

names=startlocs.SEPA_Site(I);


% Absolute relative values
subplot(1,2,1)
colormap(flipud(gray))
mx=1;
%mx=max(max(mean(c1(I,I,:),3)));
imagesc(mean(connectivity(I,I,:),3)/mx)
set(gca,'Ydir','Normal')
cb1=colorbar;
cb1.Label.String = 'P(connection)';
xlabel('destination site')
xticks(1:length(names))
xticklabels(names)
ylabel('source site')
title('(a)')
yticks(1:length(names))
yticklabels(names)
xtickangle(90)

% Logarithmic relative values
subplot(1,2,2)
colormap(flipud(gray))
mx=1;
%mx=max(max(mean(connectivity(I,I,:),3)));
imagesc(log(mean(connectivity(I,I,:),3)/mx))
set(gca,'Ydir','Normal')
cb1=colorbar;
cb1.Label.String = 'ln(P(connection))';
xlabel('destination site')
xticks(1:length(names))
xticklabels(names)
ylabel('source site')
title('(b)')
yticks(1:length(names))
yticklabels(names)
xtickangle(90)

%print('-painters','-dpng','-r600','connectivityMean_localToTest.png')


%% ---------- Connectivity - list source and destination areas for a specific site ------------
nhours=168;
nSimWeeks=29;
focalSite=160;

meanConnect=mean(connectivity,3);

% Sources
disp('Sources');
find(meanConnect(:,focalSite))
% Strength of sources
% Need to divide by nhours to get prob of single particle (already divided
% by 5 in calculation) - think this is correct
meanConnect(find(meanConnect(:,focalSite)),focalSite)/nhours 
% FMA of sources
startlocs.FMA(find(meanConnect(:,focalSite)))
% SEPA_Site of sources
startlocs.SEPA_Site(find(meanConnect(:,focalSite)))
% Proportion of weeks in which connection is made
inCon=squeeze(connectivity(find(meanConnect(:,focalSite)),focalSite,:));
inCon~=0;
sum(inCon~=0,2)/nSimWeeks

% Destinations
disp('Destinations');
find(meanConnect(focalSite,:))
% Strength of destinations
meanDest=meanConnect(focalSite,find(meanConnect(focalSite,:)))/nhours
% FMA of destinations
startlocs.FMA(find(meanConnect(focalSite,:)))
% SEPA_Site of destinations
startlocs.SEPA_Site(find(meanConnect(focalSite,:)))
% Proportion of weeks in which connection is made
outCon=squeeze(connectivity(focalSite,find(meanConnect(focalSite,:)),:));
sum(outCon~=0,2)/nSimWeeks

%% ------- Connectivity double check: what about nearby elements? -------------
focalSite=159;
startlocs(focalSite,:);
focalLoc=[startlocs.Easting(focalSite),startlocs.Northing(focalSite)];

% List elements within a certain radius
r=1000;
near=find(pdist2(focalLoc,mesh.mesh.uvnode_os)<r);

% Check for sources of particles going in those elements
for t=3:nweeks
    weekStartDate=startDate+(t-1)*7;
    weekStartDate;
    %endDate=startDate;
    endDate=weekStartDate+6;
    fprintf('Week %d - %s to %s (inclusive)\n',t,datestr(weekStartDate,'yyyymmdd'),datestr(endDate,'yyyymmdd'));

    status=2;
    elementSearch(mesh,weekStartDate,endDate,0,24,status,near);
    
end





%% Check unique site IDs in output
startDateString='20130802';
startDate=datetime(startDateString,'InputFormat','yyyyMMdd');
endDate=startDate+6;

arrivalfiles=dir('arrivals*.out');
% pLocations will contain the files to be used in the element count
pArrivals = {};
daysFound=0;
for date=startDate:endDate
    for i=1:size(arrivalfiles)
        if (strcmp(['arrivals_',datestr(date,'yyyymmdd'),'.out'],arrivalfiles(i).name)==1)
            %fprintf([arrivalfiles(i).name '\n']);
            daysFound=daysFound+1;
            pArrivals{daysFound}=fullfile(pwd,arrivalfiles(i).name);
        end
    end
end

dsArrivals = datastore(pArrivals,'Type','tabulartext','FileExtensions','.out','ReadVariableNames',1);
dsArrivals.Delimiter = ' ';
dsArrivals.NumHeaderLines = 0;
dsArrivals.TextscanFormats = {'%d','%{yyyyMMdd}D','%f','%s','%{yyyyMMdd}D','%f','%s','%f','%f'};
dsArrivals.ReadSize= 500000;
dsArrivals.SelectedVariableNames = {'startLocation'};

u={};
counter=0
reset(dsArrivals)
while hasdata(dsArrivals)
%while counter < 3
    counter=counter+1
    locationDataChunk = read(dsArrivals);
    u1=unique(locationDataChunk.startLocation);
    if ~isempty(u)
        u=unique(vertcat(u,u1));
    else
        u=u1;
    end
end

%% Get max density of lice
x = [mesh.mesh.nodexy_os(mesh.mesh.trinodes(:,1),1) mesh.mesh.nodexy_os(mesh.mesh.trinodes(:,2),1) mesh.mesh.nodexy_os(mesh.mesh.trinodes(:,3),1)];
y = [mesh.mesh.nodexy_os(mesh.mesh.trinodes(:,1),2) mesh.mesh.nodexy_os(mesh.mesh.trinodes(:,2),2) mesh.mesh.nodexy_os(mesh.mesh.trinodes(:,3),2)];
area = 0.5*abs((x(:,1)-x(:,3)).*(y(:,2)-y(:,1))-(x(:,1)-x(:,2)).*(y(:,3)-y(:,1)));

densityExisting_cop = (countsExisting_cop*(fecundityDaily/hoursPerDay)*(400/2300))./repmat(area,1,1);

% Maximum over entire mesh
max(max(densityExisting_cop,[],1))

% Maximum within range of site
focalSite=159;
startlocs(focalSite,:);
focalLoc=[startlocs.Easting(focalSite),startlocs.Northing(focalSite)];
r=1000;
near=find(pdist2(focalLoc,mesh.mesh.uvnode_os)<r);
max(max(densityExisting_cop(near,:),[],1))

focalSite=160;
startlocs(focalSite,:);
focalLoc=[startlocs.Easting(focalSite),startlocs.Northing(focalSite)];
r=1000;
near=find(pdist2(focalLoc,mesh.mesh.uvnode_os)<r);
max(max(densityExisting_cop(near,:),[],1))



