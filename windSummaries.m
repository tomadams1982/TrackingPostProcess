% Read in Dima's wind average data and plot

load('D:\Dima_wind\PWj_TIDES01_jp79244_20130620-1802202.mat')
loadMeshAndSites

% Find nearest centroids to sites
addpath('C:\Users\SA01TA\Documents\code\matlab\particletracking');
dsLocations = datastore('C:\Users\SA01TA\Documents\Sealice_NorthMinch\Kames2018\all_sites_header.dat',...
    'Type','tabulartext','ReadVariableNames',1);
dsLocations.Delimiter = '\t';
dsLocations.NumHeaderLines = 0;
dsLocations.TextscanFormats = {'%s','%d','%d','%f','%f','%f','%f','%s','%s','%s'};
sites = readall(dsLocations);

sites = sites(157:160,:);

for i=1:4
    [nID(i),dist(i)]=nearest_centroid(sites.Easting(i),sites.Northing(i),mesh.mesh.uvnode_os);
end

startdate=datetime(2013,7,1);
dates_month=startdate+calmonths([0:56]);
startdate_wk=datetime(2013,6,20);
dates_wk=startdate_wk+7*[0:242];

siteNames={'Loch Pooltiel','West Jura','Test 1','Test 2'};


%% Make plots for a given site

for i=4
    % calculate direction
    [t,r]=cart2pol(TIDES.u_7d(nID(i),:),TIDES.v_7d(nID(i),:));
    t=pi/2-t; % polar to compass conversion
    t=t+pi; % convert "to" directions to "from"
    t=mod(t,2*pi);

    [t_ma,r_ma]=cart2pol(TIDES.mean_u_ma(nID(i),:),TIDES.mean_v_ma(nID(i),:));
    t_ma=pi/2-t_ma; % polar to compass conversion
    t_ma=t_ma+pi; % convert "to" directions to "from"
    t_ma=mod(t_ma,2*pi);

%     % plot for site 1 (KPOOL1) - weekly
%     subplot(4,1,1)
%     plot(dates_wk,TIDES.u_7d(nID(i),:))
%     title(siteNames{i})
%     xlabel('date');
%     ylabel('u (m/s)');
%     subplot(4,1,2)
%     plot(dates_wk,TIDES.v_7d(nID(i),:))
%     xlabel('date');
%     ylabel('v (m/s)');
%     subplot(4,1,3)
%     plot(dates_wk,t(:))
%     xlabel('date');
%     ylabel('direction from (rad)');
%     subplot(4,1,4)
%     plot(dates_wk,r(:))
%     xlabel('date');
%     ylabel('wind speed (m/s)');
%     print('-painters','-dpng','-r600',['windWeekly_Site' int2str(i) '.png'])
% 
    % plot for site 1 (KPOOL1) - weekly 2016
%     subplot(4,1,1)
%     plot(dates_wk(133:185),TIDES.u_7d(nID(i),133:185))
%     title(siteNames{i})
%     xlabel('date');
%     ylabel('u (m/s)');
%     subplot(4,1,2)
%     plot(dates_wk(133:185),TIDES.v_7d(nID(i),133:185))
%     xlabel('date');
%     ylabel('v (m/s)');
    subplot(2,1,1)
    plot(dates_wk(133:185),t(133:185))
    xlabel('date');
    ylabel('direction from (rad)');
    subplot(2,1,2)
    plot(dates_wk(133:185),r(133:185))
    xlabel('date');
    ylabel('wind speed (m/s)');
    print('-painters','-dpng','-r600',['windWeekly2016_2_Site' int2str(i) '.png'])

    % plot for site 1 (KPOOL1) - monthly
%     subplot(4,1,1)
%     plot(dates_month,TIDES.mean_u_ma(nID(i),:))
%     xlabel('date');
%     ylabel('u (m/s)');
%     title(siteNames{i})
%     subplot(4,1,2)
%     plot(dates_month,TIDES.mean_v_ma(nID(i),:))
%     xlabel('date');
%     ylabel('v (m/s)');
    subplot(2,1,1)
    plot(dates_month,t_ma(:))
    xlabel('date');
    ylabel('direction from (rad)');
    subplot(2,1,2)
    plot(dates_month,r_ma)
    xlabel('date');
    ylabel('wind speed (m/s)');
    print('-painters','-dpng','-r600',['windMonthly_2_Site' int2str(i) '.png'])
end

%%








