function writeJONSWAP(Hsig,Tp,fpath,wave_angle)
%Writes jonswap file necessary for Xbeach
%   Given wave height and period, this function creates the necessary
%   jonswap file for xbeach
%   wave_angle = angle of wave attack in nautical convention.  Coming from 

% Name of file
fname = 'jonswap';

% Open File
fid = fopen(fname,'w');

% Write text to file
fprintf(fid,'Hm0        = %.2f\n',Hsig); % Significant Wave Height [m]
fprintf(fid,'fp         = %.2f\n',1/Tp); % Peak frequency of wave spectrum [s-1] 
fprintf(fid,'mainang    = %2.1f\n',wave_angle); % Main wave angle (nautical convention): range = 180 - 360      Note: Nautical is direction where waves come from, measured clockwise from North
fprintf(fid,'gammajsp   = 3.3000\n'); % Peak enhancement factor in JONSWAP expression 
fprintf(fid,'s          = 20.0000\n'); % Direction spreading coefficient
fclose(fid);

% Move File
movefile('jonswap',fpath)

end

