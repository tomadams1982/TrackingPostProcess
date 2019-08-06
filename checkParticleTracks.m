
% Run "plotDomains.m", first code block

% Remember - need to rename file extension to "dat" in order that Matlab
% can recognise it
%T1=readtable("C:\Users\SA01TA\Documents\particle_track\locations_20181227.dat");

basedir='C:\Users\SA01TA\Documents';
%basedir='C:\Users\SA01TA\OneDrive - SAMS\Documents\';

%d='C:/Users/SA01TA/Documents/particle_track/output/20190228_ffTest_1domain_noSlip/';
%d='C:/Users/SA01TA/Documents/particle_track/output/20190228_ffTest_1domain_noSlip_euler25/';
%d='C:/Users/SA01TA/Documents/particle_track/output/20190228_ffTest_1domain_backToCentroid_euler25/';
% COMPASS - run first two blocks of COMPSS/plotDomains.m first
%d='C:/Users/SA01TA/Documents/COMPASS/runTest/'; 
%d='C:\Users\SA01TA\Documents\COMPASS\runTest\14day\'
%d=[basedir '\COMPASS\runTest\boundaryTest4\']
d=[basedir '\Sealice_NorthMinch\190517_SSF_Linnhe\testOutput\']
% FVCOM only - run first block of COMPASS/plotDomains.m first
%d='W:/sa01ta/20190328_6_dye_6hrRelease/';


fileList=dir([d '/locations_*']);
Tall=readtable([d,fileList(1).name]);
for f=2:length(fileList)
    Tall=[Tall;readtable([d,fileList(f).name])];
end

%%
% Compa

%ids=1:156; % fishfarms
%ids=1:2644; % COMPASS nephrops test 1
ids=15:16:17439;


hold all
scatter(Tall.x(ids+1),Tall.y(ids+1))
% for part=0:155
%     plot(Tall.x(Tall.ID==part),Tall.y(Tall.ID==part),'b')
% end
for part=1:length(ids)
    if startsWith(Tall.startLocation(Tall.ID==ids(part)),'ROMS')
        plot(Tall.x(Tall.ID==ids(part)),Tall.y(Tall.ID==ids(part)),'r')
    else
        plot(Tall.x(Tall.ID==ids(part)),Tall.y(Tall.ID==ids(part)),'b')
    end
end

%%
print('-painters','-dpng','-r600',[basedir 'COMPASS\' 'boundaryTestTracks_ROMSbuffer.png'])