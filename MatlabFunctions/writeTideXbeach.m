function writeTideXbeach(time_stop,water_level,run_fol)
%Creates the tide.txt file needed for xbeach
%   time_stop: Model time in seconds for end of simulation 
%   water_level: Water level for simulation [m]
%   run_fol:  Path to folder where run is occuring 

fname = 'tide.txt';
fid = fopen(fname,'w');

time_stop = time_stop + 420; % Give us a little more time to ensure model runs 
time = 0;

for ii = 1:2
    fprintf(fid,'%1.8e  %1.8e  %1.8e\n',time,water_level,water_level);
    time = time+time_stop;
end
fclose(fid);
if ~exist(run_fol)
    mkdir(run_fol)
    movefile(fname,run_fol)
else
    movefile(fname,run_fol)
end
end



