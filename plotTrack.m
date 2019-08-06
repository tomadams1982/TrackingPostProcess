locations1 = readtable('C:\Users\SA01TA\Documents\particle_track\locations_20171231_2.dat');

plot(locations1.x(locations1.ID==0),locations1.y(locations1.ID==0))
hold on
scatter(-5.6738,57.3783);
