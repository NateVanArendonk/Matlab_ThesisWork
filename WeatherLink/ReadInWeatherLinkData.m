fname = 'test.txt'; % Name of text file exported from weatherlink program
fol = 'C:\Users\ahooshmand\Desktop\test\'; % Relative path of folder housing text file 
fid = fopen([fol fname]);
% Read in the data - Delimiter is white space and its somewhat variable
% hence the many different types, I really think theres only two different
% sizes but it works so I'll keep them all 
data = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %*[^\n]','headerLines',2','Delimiter',{' ','	','	','	','	','	'});
fclose(fid);

% Here is the order of the data 
% 1. Date
% 2. Time
% 3. AM/PM
% 4. Temp Outside (Instrument Temp)
% 5. Hi Temp 
% 6. Low Temp
% 7. Out Humidity 
% 8. Dew Point 
% 9. Wind Speed 
% 10. Wind Direction 
% 11. Wind Run 
% 12. Hi speed 
% 13. Hi Dir 
% 14. Wind Chill 
% 15. Heat Index  
% 16. THW Index 
% 17. Pressure

dataWant = [1,2,3,9,10,17];
data = data(dataWant);

% Put in a structure 

% NEED TO FIND WAY TO CONVERT FROM ENE TO DEGREES 

fclose('all');

