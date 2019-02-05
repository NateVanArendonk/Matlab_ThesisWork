
clearvars
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
    
    % Find other stations if they exist
    if ~contains(completed_names,name2find) 
        inds = [ii]; % to house all the same station indices
        for jj = 1:length(stns)
            % Get the .mat file name
            [~,temp_name,~] = fileparts(stns(jj).name);
            % if it is the same name as the current station we are looking
            % at, save the indice
            if contains(temp_name,name2find) && jj ~= ii
                inds(end+1) = jj;
            end
        end
        % Plot each station
        clf
        title(name2find)
        for jj = 1:length(inds)
            fname = sprintf('station_data/%s',stns(inds(jj)).name);
            S = load(fname);
            S.datenum = datenum(S.yr, S.mo, S.da, S.hr, S.mn, 30);
            plot(S.datenum,S.spd)
            %             plot(S.lon,S.lat,'*','MarkerSize',9)
            hold on
        end
        ylabel('Wind Speed [mph]')
        xlabel('Time [years]')
        datetick()
        title(name2find)
        fig_name = sprintf('%s_Speed_PreProcessed.png',name2find);
        printFig(gcf,fig_name,[11 11],'png',300)
        movefile(fig_name,'PrePlots')
        completed_names = [completed_names name2find];
        pause
    end
end
