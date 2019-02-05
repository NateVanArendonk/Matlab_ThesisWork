function r = makePixelInstsDemo( )

% Make pixel instruments using pixel toolbox. duplicate of demoInstsFile.

% start by forgetting everything you know
PIXForget;
% and set the station name
PIXSetStation('whidbey');

% sea level
zmsl = 0;

instID = [];

% % vbar insts, two y transects 
% y = [450 700];
% x = [125:25:225];
% 
% for ii=1:length(x)
%     
%     name = ['vbar' num2str(x(ii))];
%     tid = PIXMakeVBAR( x(ii), y(1), x(ii), y(2), zmsl, name);
%     instID = [instID tid];
%     
% end

% runup
% note: normally for a runup line you want to keep a fixed X and Y.
% This requires knowing the Z for each point on the line. For a fixed
% station you may have done the survey transects or have other information
% that allows you to know that. For a UAV operation, you probably do not,
% Thus this is the "no Z" call.

% These are the requested instruments.... 
%       insts(1).xyz = [25 -15 0; 130 -15 0]; 
%       insts(2).xyz = [18 -25.1 0; 130 -25.1 0]; 
%       insts(3).xyz = [0 -35.2 0; 130 -35.2 0]; 
%       insts(4).xyz = [0 -45.3 0; 130 -45.3 0]; 
%       insts(5).xyz = [0 -55.4 0; 130 -55.4 0]; 
%       insts(6).xyz = [0 -65.3 0; 130 -65.5 0];  
      
      
xshoreOFF = [100 100 100 100 100 100]'      
xshoreON = [25 18 0 0 0 0]'
ys = [-15 -25.1 -35.2 -45.3 -55.4 -65.3]';
rot = 0;
zmsl = 0;
for ii=1:length(ys)
    xshore = xshoreOFF(ii);
    xdune = xshoreON(ii);
    y = ys(ii);
    [xr, yr, zr] = runupPIXArray( xshore, xdune, y, zmsl, rot );
    name = ['runup' num2str(fix(y))];
    tid = PIXCreateInstrument( name, 'runup', PIXFixedZ+PIXUniqueOnly );
    PIXAddPoints( tid, xr, yr, zr, name );
    instID = [instID tid];
end

% % cBathy, 5 meter points
% dx = 5;
% dy = 5;
% x1 = 80;
% x2 = 400;
% y1 = 450;
% y2 = 900;
% 
% tid = PIXCreateInstrument( 'mBW', 'matrix', PIXFixedZ+PIXDeBayerStack );
% PIXAddMatrix( tid, x1, y1, x2, y2, zmsl, dx, dy );
% instID = [instID tid];

% % stability slices
% tid = PIXCreateInstrument( 'x300Slice', 'line', PIXFixedZ+PIXUniqueOnly );
% PIXAddLine( tid, 300, 500, 300, 540, 7, .1, 'x300Slice' );
% instID = [instID tid];
% tid = PIXCreateInstrument( 'y517Slice', 'line', PIXFixedZ+PIXUniqueOnly );
% PIXAddLine( tid, 100, 520, 115, 520, 7, .1, 'y517Slice' );
% instID = [instID tid];
% 
% 
% make a package of all insts
pid = PIXCreatePackage('whidbey', instID);
e = matlab2Epoch(now);

% build the initial r
r = PIXCreateR( pid, e, zmsl);

end

