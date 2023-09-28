function x = ImportTS(outdir,groupname)
    
    base= outdir;
    gn= groupname;
    if ~exist(fullfile(base,'TimeMats'),'dir')
        mkdir(fullfile(base,'TimeMats'))
        mkdir(fullfile(base,'TimeMats',gn))
    end
    if ~exist(fullfile(base,'TimeMats',gn),'dir')
        mkdir(fullfile(base,'TimeMats',gn))
    end
    prompt = questdlg('Please Select directory with time series files','Select',...
        'Select Dir','Cancel','Cancel');
    switch prompt
        case 'Select Dir'
            TSdir = uigetdir(fullfile(base));
            if TSdir == 0
                return
            end
        case 'Cancel'
            return
    end
    Import;
    x=1;

    
    function Import
        subfiledir=dir(fullfile(TSdir,'sub*'));
        savedir=fullfile(outdir, 'TimeMats', gn,'\');
        
        for i = 1:length(subfiledir)
            subfiles(i) = convertCharsToStrings(subfiledir(i).name);
        end
        for i = subfiles
            c1= convertStringsToChars(i);
            c2 = str2num(c1(10));
        end
        for j = subfiles
            timemats= readmatrix(fullfile(TSdir, convertStringsToChars(j)));
            for k = 1:size(timemats,2)
                timemats(:,k)=normalize(detrend(timemats(:,k),3));
            end
            Global_ts = mean(timemats,2);
            j2=strsplit(j,'.');
            sub =j2{1};
            filename = [savedir sub(1:end-1) '_TS.mat'];
            save(filename,'timemats','Global_ts')
            disp([sub(1:10) ' Completed Import'])
            clearvars timemats
        end
    end
    
        
    
end

