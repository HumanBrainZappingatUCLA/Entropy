function x = BENchange_calc(basedir,groupname)

    subfiledir=dir(fullfile(basedir,'BEN', groupname, 'BENoutputs\sub*'));
    filedir=fullfile(basedir,'BEN', groupname, 'BENoutputs');
    savefile= groupname;
    
    savedir= fullfile(basedir, 'BEN',groupname,'Correlations');
    for i = 1:length(subfiledir)
        subfiles(i) = convertCharsToStrings(subfiledir(i).name);
    end
    for i = subfiles
        c1= convertStringsToChars(i);
        c2 = str2num(c1(end-4));
        %c2 = str2num(c1(10));
        c3a = [c1(1:end-5) num2str(c2+1) c1(end-3:end)];
        c3b = [c1(1:end-5) num2str(c2-1) c1(end-3:end)];
        if mod(c2,2) == 0
            if ~strcmp(subfiles,convertCharsToStrings(c3b))
                ind = find(subfiles == i);
                subfiles(ind) = [];
            end
        elseif mod(c2,2) == 1
            if ~strcmp(subfiles,convertCharsToStrings(c3a))
                ind = find(subfiles == i);
                subfiles(ind) = [];
            end
        end
    end


    k1=1;
    k2=1;
    for n = subfiles
        n2 = convertStringsToChars(n);
        n3 = str2num(n2(end-4));
        load(fullfile(filedir, n2))

        p=1;
        ROInum=size(SENroi,1);
        for i=1:ROInum
            BEN_array(1,i)=SENroi(i,1);
        end
        if mod(n3,2) == 0
            BEN_After(k2,:) = BEN_array;
            BENgs_After(k2,:)=SENglobal;
            k2=k2+1;
        else
            BEN_Before(k1,:) = BEN_array;
            BENgs_Before(k1,:)=SENglobal;
            BENsubs(k1,1)=convertCharsToStrings(n2(1:9));
            k1=k1+1;
        end

        clearvars SENroi SENglobal BEN_array

    end
    
    BENchange=BEN_After-BEN_Before;
    filename=fullfile(savedir, [savefile '_BENchange.mat']);
    save(filename,'BENsubs','BENchange')
    x=0;
end