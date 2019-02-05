clearvars
addpath C:\Functions_Matlab
stns = dir('station_data/*.mat');
completed_names = ['z'];


% Loop through stations
for ii = 1:length(stns)
    % Grab the name of the current station
    cur_name = stns(ii).name;
    % Get rid of any spaces
    under_ind = strfind(cur_name,'_');
    if ~isempty(under_ind)
        name2find = cur_name(1:under_ind-1);
    else 
        period_ind = strfind(cur_name,'.');
        name2find = cur_name(1:period_ind-1);
    end
    
    % Find other stations if they exist and record their indice
    if ~contains(completed_names,name2find) % If the station of interest hasn't already been concatenated
        inds = [ii]; % house all the same station indices
        for jj = 1:length(stns)
            % Get the .mat file name
            [~,temp_name,~] = fileparts(stns(jj).name);
            % if it is the same name as the current station we are looking
            % at, save the indice
            if contains(temp_name,name2find) && jj ~= ii
                inds(end+1) = jj;
            end
        end
        
        % Load each station into a structure 
        for jj = 1:length(inds)
            fname = sprintf('station_data/%s',stns(inds(jj)).name);
            S(jj) = load(fname);
        end
        
        % Now add time vector as a datenum
        for jj = 1:length(S)
            S(jj).time = datenum(S(jj).yr, S(jj).mo, S(jj).da, S(jj).hr, S(jj).mn, 30);
        end
        
        % Concatenate time series 
        T = S(1);
        if length(S) > 1
            for jj = 2:length(S)
                T.usaf = [T.usaf; S(jj).usaf];
                T.wban = [T.wban; S(jj).wban];
                T.yr = [T.yr; S(jj).yr];
                T.mo = [T.mo; S(jj).mo];
                T.da = [T.da; S(jj).da];
                T.hr = [T.hr; S(jj).hr];
                T.mn = [T.mn; S(jj).mn];
                T.wnddir = [T.wnddir; S(jj).wnddir];
                T.spd = [T.spd; S(jj).spd];
                T.gust = [T.gust; S(jj).gust];
                T.slp = [T.slp; S(jj).slp];
                T.alt = [T.alt; S(jj).alt];
                T.stp = [T.stp; S(jj).stp];
                T.time = [T.time; S(jj).time];
                T.lon = [T.lon; S(jj).lon];
                T.lat = [T.lat; S(jj).lat];
                T.elev = [T.elev; S(jj).elev];
            end
        end
        
        % Sort the time series to get it in chronological order
        [~,I] = sort(T.time,'ascend');
        T.usaf = T.usaf(I);
        T.wban = T.wban(I);
        T.yr = T.yr(I);
        T.mo = T.mo(I);
        T.da = T.da(I);
        T.hr = T.hr(I);
        T.mn = T.mn(I);
        T.wnddir = T.wnddir(I);
        T.spd = T.spd(I);
        T.gust = T.gust(I);
        T.slp = T.slp(I);
        T.alt = T.alt(I);
        T.stp = T.stp(I);
        T.time = T.time(I);
        
        % Now Plot
        clf
        plot(T.time,T.spd)
        datetick()
        xlabel('Time [years]')
        ylabel('Wind Speed [mph]')
        title(name2find)
        fig_name = sprintf('%s_Speed_PostConcate.png',name2find);
        printFig(gcf,fig_name,[11 11],'png',300)
        movefile(fig_name,'PostConcatenatePlots')
        completed_names = [completed_names name2find];
        
        % Save the File 
        save(name2find,'-struct','T');
        save_nm = sprintf('%s.mat',name2find);
        movefile(save_nm,'combined_station_data')
        clear S T I    
    end
    fprintf('%s is completed moving on\n',name2find)
end