function x = RunEntropyROI(basedir,groupname,TSdir1,TSdir2,Transform_mat)
    
    
    subdir = dir(fullfile(basedir, 'TimeMats', groupname, 'sub*'));
    subfiledir1=dir(fullfile(TSdir1, 'sub*'));
    subfiledir2=dir(fullfile(TSdir2, 'sub*'));
    savedir = fullfile(basedir, 'BEN', groupname);
    rprime=0.3;
    %workdir = fullfile(basedir, 'TimeMats', groupname);
    

    for i = 1:length(subfiledir1)
        subs1(i)= convertCharsToStrings(subfiledir1(i).name);
        %subs(1)= convertCharsToStrings(subdir(i).name);
    end
    for i = 1:length(subfiledir2)
        subs2(i)= convertCharsToStrings(subfiledir2(i).name);
        %subs(1)= convertCharsToStrings(subdir(i).name);
    end
    for j = subs1
        j2= convertStringsToChars(j);
        load(fullfile(TSdir1, j2));
        dat = timemats(:,Transform_mat)';
        inmaskdatfile=fullfile(savedir,'braindats', [j2(1:end-6) '_braindat.dat']); %path for braindat.dat creation
        std_d=std(dat,0,2);
        dat=dat';                         % reorganize data to be time series by voxels. This is important.
        loc=1:size(dat,2);
        r=rprime;    % if the BEN map is nearly zero everywhere, increase r to be 0.6 here
        ndat=dat;     
        tlen=size(ndat,1);
        sum=0;
        %% Region SEN
        for dim=[3]
            for scale=0
                tlen=size(ndat,1);
                fprintf('\t\t\t saving inmask data for calc SEN.\n');
                fid=fopen(inmaskdatfile,'wb');
                fwrite(fid, length(loc), 'int');   % number of timecourses
                fwrite(fid, tlen, 'int');
        
                fwrite(fid, std_d, 'float');
                fwrite(fid, ndat, 'float');
                fclose(fid);
                if isunix
                    senfile=fullfile(savedir,'SENdats',[j2(1:end-6) '_SEN.dat']);
                    str=['!SampEn -d ' num2str(dim) ' -r ' num2str(r) ' -i ' inmaskdatfile  ' -o ' senfile];
                    eval(str);
                else
                    senfile=fullfile(savedir,'SENdats',[j2(1:end-6) '_SEN.dat']);
                    str=['SampEn.exe -d ' num2str(dim) ' -r ' num2str(r) ' -i ' inmaskdatfile  ' -o ' senfile];
                    system(str);
                end
                fid=fopen(senfile,'rb');
                slen=fread(fid, 1, 'int');
                dim=fread(fid, 1, 'int');
                ratio=fread(fid, 1, 'float');
                sSEN=fread(fid, slen, 'float');
                fclose(fid);

            end

        end

        SENroi = sSEN;
        clearvars -except SENroi savedir j2 j subs1 subs2 TSdir1 TSdir2 Global_ts Transform_mat rprime
        dat2 = Global_ts';
        inmaskdatfile=fullfile(savedir,'braindats', [j2(1:end-6) '_braindatGlobal.dat']); %path for braindat.dat creation
        std_d=std(dat2,0,2);
        dat2=dat2';                         % reorganize data to be time series by voxels. This is important.
        loc=1:size(dat2,2);
        r=rprime;    % if the BEN map is nearly zero everywhere, increase r to be 0.6 here
        ndat=dat2;     
        tlen=size(ndat,1);
        sum=0;       
        for dim=[3]
            for scale=0
                tlen=size(ndat,1);
                fprintf('\t\t\t saving inmask data for calc SEN.\n');
                fid=fopen(inmaskdatfile,'wb');
                fwrite(fid, length(loc), 'int');   % number of timecourses
                fwrite(fid, tlen, 'int');
        
                fwrite(fid, std_d, 'float');
                fwrite(fid, ndat, 'float');
                fclose(fid);
                if isunix
                    senfile=fullfile(savedir,'SENdats',[j2(1:end-6) '_SENGlobal.dat']);
                    str=['!SampEn -d ' num2str(dim) ' -r ' num2str(r) ' -i ' inmaskdatfile  ' -o ' senfile];
                    eval(str);
                else
                    senfile=fullfile(savedir,'SENdats',[j2(1:end-6) '_SENGlobal.dat']);
                    str=['SampEn.exe -d ' num2str(dim) ' -r ' num2str(r) ' -i ' inmaskdatfile  ' -o ' senfile];
                    system(str);
                end
                fid=fopen(senfile,'rb');
                slen=fread(fid, 1, 'int');
                dim=fread(fid, 1, 'int');
                ratio=fread(fid, 1, 'float');
                sSEN=fread(fid, slen, 'float');
                fclose(fid);

            end

        end
        SENglobal = sSEN;
        sub = split(j,'_');
        filename = fullfile(savedir, 'BENoutputs', [sub{1} '_BENoutputs1.mat']);
        save(filename,'SENroi','SENglobal');
        disp([sub{1} ' Completed SampEn'])
        clearvars -except savedir j subs1 subs2 TSdir1 TSdir2 subfiledir1 subfiledir2 Transform_mat rprime
    end
clear j
    for j = subs2
        j2= convertStringsToChars(j);
        load(fullfile(TSdir2, j2));
        dat = timemats(:,Transform_mat)';
        inmaskdatfile=fullfile(savedir,'braindats', [j2(1:end-6) '_braindat2.dat']); %path for braindat.dat creation
        std_d=std(dat,0,2);
        dat=dat';                         % reorganize data to be time series by voxels. This is important.
        loc=1:size(dat,2);
        r=rprime;    % if the BEN map is nearly zero everywhere, increase r to be 0.6 here
        ndat=dat;     
        tlen=size(ndat,1);
        sum=0;
        %% Region SEN
        for dim=[3]
            for scale=0
                tlen=size(ndat,1);
                fprintf('\t\t\t saving inmask data for calc SEN.\n');
                fid=fopen(inmaskdatfile,'wb');
                fwrite(fid, length(loc), 'int');   % number of timecourses
                fwrite(fid, tlen, 'int');
        
                fwrite(fid, std_d, 'float');
                fwrite(fid, ndat, 'float');
                fclose(fid);
                if isunix
                    senfile=fullfile(savedir,'SENdats',[j2(1:end-6) '_SEN2.dat']);
                    str=['!SampEn -d ' num2str(dim) ' -r ' num2str(r) ' -i ' inmaskdatfile  ' -o ' senfile];
                    eval(str);
                else
                    senfile=fullfile(savedir,'SENdats',[j2(1:end-6) '_SEN2.dat']);
                    str=['SampEn.exe -d ' num2str(dim) ' -r ' num2str(r) ' -i ' inmaskdatfile  ' -o ' senfile];
                    system(str);
                end
                fid=fopen(senfile,'rb');
                slen=fread(fid, 1, 'int');
                dim=fread(fid, 1, 'int');
                ratio=fread(fid, 1, 'float');
                sSEN=fread(fid, slen, 'float');
                fclose(fid);

            end

        end

        SENroi = sSEN;
        clearvars -except SENroi savedir j2 j subs1 subs2 TSdir1 TSdir2 Global_ts Transform_mat rprime
        dat2 = Global_ts';
        inmaskdatfile=fullfile(savedir,'braindats', [j2(1:end-6) '_braindatGlobal2.dat']); %path for braindat.dat creation
        std_d=std(dat2,0,2);
        dat2=dat2';                         % reorganize data to be time series by voxels. This is important.
        loc=1:size(dat2,2);
        r=rprime;    % if the BEN map is nearly zero everywhere, increase r to be 0.6 here
        ndat=dat2;     
        tlen=size(ndat,1);
        sum=0;       
        for dim=[3]
            for scale=0
                tlen=size(ndat,1);
                fprintf('\t\t\t saving inmask data for calc SEN.\n');
                fid=fopen(inmaskdatfile,'wb');
                fwrite(fid, length(loc), 'int');   % number of timecourses
                fwrite(fid, tlen, 'int');
        
                fwrite(fid, std_d, 'float');
                fwrite(fid, ndat, 'float');
                fclose(fid);
                if isunix
                    senfile=fullfile(savedir,'SENdats',[j2(1:end-6) '_SENGlobal2.dat']);
                    str=['!SampEn -d ' num2str(dim) ' -r ' num2str(r) ' -i ' inmaskdatfile  ' -o ' senfile];
                    eval(str);
                else
                    senfile=fullfile(savedir,'SENdats',[j2(1:end-6) '_SENGlobal2.dat']);
                    str=['SampEn.exe -d ' num2str(dim) ' -r ' num2str(r) ' -i ' inmaskdatfile  ' -o ' senfile];
                    system(str);
                end
                fid=fopen(senfile,'rb');
                slen=fread(fid, 1, 'int');
                dim=fread(fid, 1, 'int');
                ratio=fread(fid, 1, 'float');
                sSEN=fread(fid, slen, 'float');
                fclose(fid);

            end

        end
        SENglobal = sSEN;
        sub = split(j,'_');
        filename = fullfile(savedir, 'BENoutputs', [sub{1} '_BENoutputs2.mat']);
        save(filename,'SENroi','SENglobal');
        disp([sub{1} ' Completed SampEn'])
        clearvars -except savedir j subs1 subs2 TSdir1 TSdir2 subfiledir1 subfiledir2 Transform_mat rprime
    end
    x=[];
end