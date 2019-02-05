function meta = importascii()
% The basic layout of this GUI was made with the help of guidegetter,
% available on the File Exchange at Mathworks.com

hf2 = figure('units','normalized',...
    'position',[0.332 0.569 0.197 0.351],...
    'menubar','none','name','importascii',...
    'numbertitle','off','color',[0.94 0.94 0.94]);


uicontrol(hf2,'style','pushbutton','units','normalized',...
    'position',[0.122 0.773 0.179 0.114],...
    'string','Open','backgroundcolor',[0.94 0.94 0.94],...
    'callback',@importascii_getfile);
gd2.filestr = uicontrol(hf2,'style','text','units','normalized',...
    'position',[0.37 0.8 0.499 0.0579],...
    'string','No file selected',...
    'horizontalalignment','left',...
    'backgroundcolor',[0.94 0.94 0.94]);


uicontrol(hf2,'style','text','units','normalized',...
    'position',[0.0596 0.650 0.323 0.0601],...
    'string','Grid Type','backgroundcolor',[0.94 0.94 0.94]);
gd2.gtype = uicontrol(hf2,'style','popupmenu','units','normalized',...
    'position',[0.432 0.646 0.375 0.0735],...
    'string',{'Cartesian';'Spherical'},'backgroundcolor',[1 1 1],...
    'callback',@importascii_gtype);

uicontrol(hf2,'style','text','units','normalized',...
    'position',[0.0596 0.579 0.323 0.0601],...
    'string','Projection','backgroundcolor',[0.94 0.94 0.94]);
gd2.prj = uicontrol(hf2,'style','edit','units','normalized',...
    'position',[0.432 0.577 0.471 0.0713],...
    'string','Undefined','backgroundcolor',[1 1 1]);

uicontrol(hf2,'style','text','units','normalized',...
    'position',[0.0521 0.492 0.323 0.0601],...
    'string','Zone','backgroundcolor',[0.94 0.94 0.94]);
gd2.zonestr = uicontrol(hf2,'style','edit','units','normalized',...
    'position',[0.432 0.492 0.471 0.0713],...
    'string','Undefined','backgroundcolor',[1 1 1]);

uicontrol(hf2,'style','text','units','normalized',...
    'position',[0.0571 0.403 0.323 0.0601],...
    'string','Horiz. Datum','backgroundcolor',[0.94 0.94 0.94]);
gd2.hdatum = uicontrol(hf2,'style','edit','units','normalized',...
    'position',[0.434 0.405 0.471 0.0713],...
    'string','Undefined','backgroundcolor',[1 1 1]);

uicontrol(hf2,'style','text','units','normalized',...
    'position',[0.0546 0.314 0.323 0.0601],...
    'string','Vert. Datum','backgroundcolor',[0.94 0.94 0.94]);
gd2.vdatum = uicontrol(hf2,'style','edit','units','normalized',...
    'position',[0.432 0.316 0.471 0.0713],...
    'string','Undefined','backgroundcolor',[1 1 1]);


uicontrol(hf2,'style','text','units','normalized',...
    'position',[0.0546 0.229 0.323 0.0601],...
    'string','Horiz. Units','backgroundcolor',[0.94 0.94 0.94]);
gd2.hunits = uicontrol(hf2,'style','popupmenu','units',...
    'normalized','position',[0.429 0.227 0.375 0.0735],...
    'string',{'m';'ft';'US Survey ft'},'backgroundcolor',[1 1 1]);

uicontrol(hf2,'style','text','units','normalized',...
    'position',[0.0521 0.156 0.323 0.0601],...
    'string','Vert. Units','backgroundcolor',[0.94 0.94 0.94]);
gd2.vunits = uicontrol(hf2,'style','popupmenu','units','normalized',...
    'position',[0.432 0.154 0.375 0.0735],...
    'string',{'m';'ft';'US Survey ft'},'backgroundcolor',[1 1 1]);


gd2.doimport=uicontrol(hf2,'style','pushbutton','units','normalized',...
    'position',[0.618 0.0312 0.251 0.0846],...
    'string','Import','backgroundcolor',[0.94 0.94 0.94],...
    'enable','off','callback',@importascii_done);


guidata(hf2,gd2)

uiwait
gd2=guidata(hf2);

meta.pathname=gd2.guipath;
meta.filename=gd2.filename;
meta.projection=get(gd2.prj,'string');
meta.zone=get(gd2.zonestr,'string');
meta.horiz_datum=get(gd2.hdatum,'string');
meta.vert_datum=get(gd2.vdatum,'string');

val=get(gd2.gtype,'val');
switch val
    case 1
        meta.ptype='cartesian';
        val2=get(gd2.hunits,'value');
        switch val2
            case 1 
                meta.horiz_units='m';
            case 2
                meta.horiz_units='ft';
            case 3
                meta.horiz_units='US ft';
        end
    case 2
        meta.ptype='spherical';
        meta.horiz_units='degrees';
end

val3=get(gd2.vunits,'val');
switch val3
    case 1
        meta.vert_units='m';
    case 2
        meta.vert_units='ft';
    case 3
        meta.vert_units='US ft';
end

close(hf2)