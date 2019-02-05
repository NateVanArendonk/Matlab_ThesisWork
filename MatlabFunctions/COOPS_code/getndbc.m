function spec = getndbc(url,tstart,tend,plotflag)
% GETNDBC - Download 2D wave spectral data from NDBC
%
% GETNDNBC(URL,TSTART,TEND,PLOTFLAG) uses the Java-NetCDF
%   toolbox to automatically download 2D wave energy and 
%   direction data from NDBC OPeNDAP servers. Directional 
%   spectra are calculated using the Maximum Likelihood 
%   Method (Earle et al.,1999).
%
% INPUT AND SYTAX:
%   SDATA = GETNDBC(URL.NC) - retrieves all the spectral wave
%       data contained in the netCDF file URL.NC. WARNING: 
%       Downloading an entire year of directional data may take 
%       a while!
%
%   SDATA = GETNDBC(URL.NC,TSTART,TEND) - retrieves information
%       between the times given in TSTART and TEND. TSTART and
%       TEND can be a string (datestr) input or a datenum.
%
%   SDATA = GETNDBC(URL.NC,TSTART,TEND,1) - produces a movie
%       of output. Default value [0] produces no graphical output.
%
% OUTPUT:
%   Data are returned in a structure with the following terms:
%           'mtime' - matlab datenum
%           'frequency' - frequency (Hz)
%           'direciton' - direction (nautical convention)
%           'dspec' - directional spectra (m*m/Hz/deg)
%           'hsig'  - significant wave height (m)
%           'tpeak' - peak wave period (s)
%           'mwd'   - peak wave direction (degrees, nautical convention
%                     specifying the direction waves are coming from)
%
% EXAMPLES:
%   1) Download a month of data between Feb and Mar for station 
%      46026 and plot results.
%
%       url=['http://dods.ndbc.noaa.gov/thredds/',...
%            'dodsC/data/swden/',...
%            '46026/46026w2010.nc'];
%       spec = getndbc(url,'02-01-2010','03-01-2010',1)
%
%   2) Download the last week of available data for station 46029 and 
%      produce a plot of the results.
%
%      url=['http://dods.ndbc.noaa.gov/thredds/',...
%           'dodsC/data/swden/',...      %the 9999 is used for 
%           '46029/46029wb9999.nc'];     %recent (unverified)data                                     
%      spec = getndbc(url,now-7,[],1)
%
% OPeNDAP NDBC  
%   A list of NDBC stations with available directional data can be
%   found here: http://dods.ndbc.noaa.gov/thredds/catalog/data/swden/
%
% JAVA-NetCDF Tools
%   The Java-NetCDF toolbox can be downloaded here: 
%   http://sourceforge.net/apps/trac/njtbx/wiki/DownloadNjtbx-current
%
% REFERENCE
%   Earle, M.D., Steele, K.E., Wang, D.W.C, 1999. Use of advanced
%       directional wave spectra analysis methods. Ocean Engineering 
%       26: 1421-1434.

% A. Stevens
% astevens@usgs.gov
% 1/24/2011

%open the netCDF object and get variables
info=ncinfo(url);

vars={'spectral_wave_density';...
    'mean_wave_dir';...
    'principal_wave_dir';...
    'wave_spectrum_r1';...
    'wave_spectrum_r2'};

%deal with the time inputs
tsecs=double(ncread(url,'time'));
times=datenum(1970,1,1)+(tsecs./86400);
if ~exist('tstart','var')
    is=1;
else
    if ischar(tstart)
        tstart=datenum(tstart);
    end
    if isempty(tstart)
        is=1;
    else
        is=find(times>=tstart,1,'first');
    end
end
if ~exist('tend','var')
    ie=numel(times);
else
    if ischar(tend)
        tend=datenum(tend);
    end
    if isempty(tend)
        ie=numel(times);
    else
        ie=find(times<=tend,1,'last');
    end
end
if isempty(is) || isempty(ie)
    str=sprintf(['No Information in specified file ',...
        'between %s and %s.'],datestr(tstart),datestr(tend));
    error(str); %#ok
end

%how many frequency bins
dims= arrayfun(@(x)(x.Name),info.Dimensions,'un',0);
fdim=strcmp('frequency',dims);
nbins=info.Dimensions(fdim).Length;

count=(ie-is)+1;
%get the raw data
rdata=cellfun(@(x)(double(squeeze(ncread(...
    url,x,[1 1 1 is],[1 1 nbins count])))'),...
    vars,'un',0);
 
%start building the output struct
metavals={'station';...
    'comment';...
    'location'};
globalatts=cellfun(@(x)(ncreadatt(url,'/',x)),....
    metavals,'un',0);

spec=struct('station',globalatts{1},...
    'comment',globalatts{2},...
    'location',globalatts{3},...
    'mtime',times(is:ie),...
    'frequency',double(ncread(url,'frequency')));

%calculate the directional spectra using
%maximum likelihood method
dires = 10; % directional bin size (degrees)
spec.direction = (0+dires/2:dires:360-(dires/2));
spec.dspec=zeros(length(spec.frequency),...
    length(spec.direction),length(spec.mtime));
for i=1:(ie-is)
    for j=1:length(spec.frequency)
        [dirdismlm,junk]=mlmndbc(rdata{4}(i,j),... 
            rdata{5}(i,j),...
            rdata{2}(i,j),...
            rdata{3}(i,j),spec.direction);%#ok
        spec.dspec(j,:,i)=(rdata{1}(i,j).*...
            dirdismlm)./dires;
    end
end

%calculate wave parameters
if length(spec.frequency)==47
    bwidths=[0.01;repmat(0.005,13,1);...
        repmat(0.01,26,1);repmat(0.02,7,1)];
elseif length(spec.frequency)==64 %cdip
    bwidths=[repmat(0.005,15,1);0.0075;...
        repmat(0.01,48,1)];
else %for unverified data, the binning is different 
    bwidths=[repmat(0.005,13,1);...
        repmat(0.01,26,1);repmat(0.02,7,1)];
end

spec.hsig=cellfun(@(x)(4.*sqrt(sum(x'.*bwidths))),...
    num2cell(rdata{1},2));
[junk,fi]=cellfun(@(x)(max(x)),... 
    num2cell(rdata{1},2));%#ok
spec.tpeak=1./spec.frequency(fi);
spec.mwd=cellfun(@(x,y)(x(y)),...
    num2cell(rdata{2},2),...
    num2cell(fi));

%plot
if ~exist('plotflag','var')
    plotflag=0;
end

if plotflag
    dirr=[spec.direction,spec.direction(1)].*(pi/180);
    dmat=repmat(dirr,length(spec.frequency),1);
    fmat=repmat(spec.frequency,1,length(dirr));
    dspec=[spec.dspec(:,:,1),spec.dspec(:,1,1)];
    [x,y]=pol2cart(dmat,fmat);
    
    radials=[min(spec.frequency),0.1:0.1:0.4,max(spec.frequency)];
    ry=linspace(0,2*pi,100);
    [rx,ry]=cellfun(@(x)(pol2cart(ry,repmat(x,1,length(ry)))),...
        num2cell(radials),'un',0);
    [rtx,rty]=pol2cart(ones(1,length(radials)).*(pi/4),radials);
    
    
    meridians=(0:45:360-45).*(pi/180);
    [mtx,mty]=pol2cart(meridians(1:2:end),...
        ones(1,4)*max(spec.frequency));
    my=linspace(min(spec.frequency),max(spec.frequency),50);
    [mx,my2]=cellfun(@(x)(pol2cart(repmat(x,1,length(my)),my)),...
        num2cell(meridians),'un',0);
    
    [dx,dy]=pol2cart(repmat(spec.mwd(1).*(pi/180),1,length(my)),my);
    
    minx=floor(min(spec.mtime));
    maxx=ceil(max(spec.mtime));
    dint=round((maxx-minx)/4);
    
    fig=figure;
    set(gcf,'render','zb')
    ax(1)=subplot(321);
    plot(spec.mtime,spec.hsig,'k-','linewi',2)
    
    ax(2)=subplot(323);
    plot(spec.mtime,spec.tpeak,'k-','linewi',2)
    
    ax(3)=subplot(325);
    plot(spec.mtime,spec.mwd,'k.')
    
    set(ax(1:2),'xticklabel',[]);
    set(ax,'xlim',[minx maxx],...
        'xtick',(minx:dint:maxx),...
        'nextplot','add');
    labs={'H_{sig}(m)';...
        'T_{peak}(s)';...
        'Direc.(deg)'};
    cellfun(@(x,y)(set(get(x,'ylabel'),...
        'string',['\bf\it',y])),...
        num2cell(ax)',labs)
    
    yl=cellfun(@(x)(get(x,'ylim')),...
        num2cell(ax)','un',0);
    lh1=cellfun(@(x,y)(line([spec.mtime(1),...
        spec.mtime(1)],x,'parent',y,...
        'linewi',2,'color','r')),...
        yl,num2cell(ax)','un',0);
    
    datetick(ax(3),'x',19,'keepticks','keeplimits')
    
    subplot(122)
    ph=pcolor(x,y,dspec);
    shading flat
    hold on
    lh=plot(dx,dy,'r-','linewidth',2); 
    
    %figure out a reasonable clim
    cdata=sort(spec.dspec(spec.dspec~=0));
    cval=interp1(linspace(0,100,numel(cdata)),...
                     cdata,98);
    
    cellfun(@(x,y)(plot(x,y,'k')),rx,ry)
    cellfun(@(x,y)(plot(x,y,'k')),mx,my2)
    set(gca,'view',[90 -90],...
        'da',[1 1 1],...
        'visible','off',...
        'clim',[0 cval])
    colormap(flipud(hot))
    
    text(rtx,rty,cellfun(@(x)(sprintf('%0.2f',x)),...
        num2cell(radials),'un',0),...
        'horizontalalign','center',...
        'verticalalign','bot',...
        'fontsize',8,...
        'fontang','it',...
        'fontweight','bo',...
        'rotation',-45)
    th=text(mtx,mty,{'N','E','S','W'},...
        'fontsize',8,...
        'fontang','it',...
        'fontweight','bo');
    valign={'bot','m','top','m'};
    halign={'c','l','c','r'};
    set(th,{'horizontalalign'},halign',...
        {'VerticalAlignment'},valign')
    
    c1=colorbar('horiz');
    set(get(c1,'xlabel'),'string',...
        '\bf\itSpectral Density (m^2/Hz/deg)')
    
    set(fig,'units','normalized');
    xp=0.01;
    yp=0.01;
    ah=uicontrol(gcf,'style','text',...
        'units','normalized',...
        'position',[xp yp 0.5 0.05]);
     set(ah,'string','Click in figure window to begin animation.',...
         'backgroundcolor',get(fig,'color'))
    
    waitforbuttonpress
    for i=1:length(spec.mtime)
        set(ah,'visible','off')
        dspec=[spec.dspec(:,:,i),spec.dspec(:,1,i)];
        [dx,dy]=pol2cart(repmat(spec.mwd(i).*(pi/180),1,length(my)),my);
        cellfun(@(x)(set(x,'xdata',[spec.mtime(i) spec.mtime(i)])),lh1)
        set(lh,'xdata',dx,'ydata',dy);
        set(ph,'cdata',dspec)
        pause(0.1)
    end
    
end



function [dirdismlm,wavedir]=mlmndbc(r1,r2,alpha1,alpha2,wavedir)
% This program is to use the MLM to compute the wave directional
% distribution function based on the given r1,r2,alph1, and alph2.
%
% firt convert the r1,r2,alp1,alp2, to 
% obtain c22,c33,c12,c13,q12,q13. (c11 was set up as one)
%
% because in this program the co and quad are used to obtain direction 
% distribution fuction for a given frequency the direction distributioin 
% will be normalized to unity later so we can set the xk and c11 to be 
% unity and will not affect the end result. because this unity set up for 
% c11 and xk, the c22, c33, c23, q12 and q13 obtained here are only good 
% for this calculation should not be considered for other purposes.

xk=1; % wave number

% ref numbers are equation numbers in M.D. Earle et. al. (1999)
c11=1;                        % (29)
tmp1=(xk*xk)/2.;              % (30)
tmp6=cos(2*alpha2*pi/180);    
tmp7=sin(2*alpha2*pi/180);
tmp8=sin(alpha1*pi/180);
tmp9=cos(alpha1*pi/180);
c22=c11*tmp1*(1-r2*tmp6);     % (31)
c33=c11*tmp1*(1+r2*tmp6);     % (32)
c23=c11*tmp1*r2*tmp7;         % (33)
q12=-c11*xk*r1*tmp8;          % (34)
q13=-c11*xk*r1*tmp9;          % (35)
xk=((c22+c33)/c11).^0.5;      % (18)
d=c11*c22*c33-c11*(c23.^2)-c33*(q12^2)-c22*(q13^2)+2.*q12*q13*c23;
if(abs(d) >= 0.0000000001)
   temp=2*(c22*c33-c23^2)/d;
   a0=temp+(c11*c22+c11*c33-q12^2-q13^2)*(xk^2)/d;
   a1=-2*(q12*c33-q13*c23)*xk/d;
   b1=2.*(q12*c23-q13*c22)*xk/d;
   a2=0.5*(c11*c33-c11*c22-q13^2+q12^2)*(xk^2)/d;
   b2=-(c11*c23-q12*q13)*(xk^2)/d;
   %[a0 a1 a2 b1 b2]
   ndd=270-wavedir;
   denim=0.5*a0+a1*cos(ndd*pi/180)+b1*sin(ndd*pi/180)+...
      a2*cos(2*ndd*pi/180)+b2*sin(2*ndd*pi/180); % (6)
   dsprd=abs(denim.^(-1));
   sum1=sum(dsprd);
   dirdismlm=dsprd/sum1;
else
   dirdismlm=wavedir*0;
end


