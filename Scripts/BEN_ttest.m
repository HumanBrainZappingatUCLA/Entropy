function x = BEN_ttest(basedir,groupname)

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
    uncorr=1;
    corr=1;
    whichTM= size(BEN_After,2);
    if whichTM == 419
        load('ROI_TransformMatrix419.mat')
    elseif whichTM == 235
        load('ROI_TransformMatrix235.mat')
    end
    for i = 1:length(BEN_After)
        BENb=BEN_Before(:,i);
        BENa=BEN_After(:,i);
        row = find(isnan(BENb));
        BENb(row)=[];
        BENa(row)=[];
        [h,p,ci,stats] = ttest(BENa',BENb');
        eff=meanEffectSize(BENb',BENa',Paired=true,Effect="robustcohen");

        Allresults(i,1)=h;
        Allresults(i,2)=p;
        Allresults(i,3)=stats.tstat;
        Allresults(i,4)=stats.df;
        Allresults(i,5)=mean(BENa);
        Allresults(i,6)=std(BENa);
        Allresults(i,7)=mean(BENb);
        Allresults(i,8)=std(BENb);
        Databackup{1,i}=Transform_matBEN(i);
        Databackup{2,i}=[BENb BENa];
        CohenEffect(i,1)= eff.Effect;
        CohenEffect(i,2)= eff.ConfidenceIntervals(1);
        CohenEffect(i,3)= eff.ConfidenceIntervals(2);
        if h == 1
            if p < 0.05
                SigUncorrect05(uncorr,1)=Transform_matBEN(i);
                SigUncorrect05(uncorr,2)=h;
                SigUncorrect05(uncorr,3)=p;
                SigUncorrect05(uncorr,4)=stats.tstat;
                SigUncorrect05(uncorr,5)=stats.df;
                SigUncorrect05(uncorr,6)=mean(BENa);
                SigUncorrect05(uncorr,7)=std(BENa);
                SigUncorrect05(uncorr,8)=mean(BENb);
                SigUncorrect05(uncorr,9)=std(BENb);
                uncorr=uncorr+1;
            end
            if p < (0.05/ROInum)
                Sig05(corr,1)=Transform_matBEN(i);
                Sig05(corr,2)=h;
                Sig05(corr,3)=p;
                Sig05(corr,4)=stats.tstat;
                Sig05(corr,5)=stats.df;
                Sig05(corr,6)=mean(BENa);
                Sig05(corr,7)=std(BENa);
                Sig05(corr,8)=mean(BENb);
                Sig05(corr,9)=std(BENb);
                corr=corr+1;
            end
        end
    end

    [fdr,Q] = mafdr(Allresults(:,2));
    [fdr2] = mafdr(Allresults(:,2),'BHFDR','true');

    ind = find(fdr2 < 0.05);
    if ~isempty(ind)
        SigMatrix05(ind) = fdr2(ind);
        for i=1:length(ind)
            Sig05Storey(i,1)=Transform_matBEN(ind(i));
            Sig05Storey(i,2)=fdr2(ind(i));
            Sig05Storey(i,3)=Allresults(ind(i),3);
        end
            
    end
    if exist('Sig05Storey')
        k=1;
        for i=1:size(Sig05Storey,1)
            if Sig05Storey(i,1) > 91 && Sig05Storey(i,1) < 201 || Sig05Storey(i,1) > 294
                %if Sig05Storey(i,2) < 0.02
                    Sig05Care(k,:)=Sig05Storey(i,:);
                    k=k+1;
                %end
            end
        end
    end
    if exist('SigUncorrect05')
        nets=SigUncorrect05(:,1);
        n11 = find(nets < 32);
        n12 = find(nets >= 201 & nets <= 230);
        n21 = find(nets >= 32 & nets <= 68);
        n22 = find(nets >= 231 & nets <= 270);
        n31 = find(nets >= 69 & nets <= 91);
        n32 = find(nets >= 271 & nets <= 293);
        n41 = find(nets >= 92 & nets <= 113);
        n42 = find(nets >= 294 & nets <= 318);
        n51 = find(nets >= 114 & nets <= 126);
        n52 = find(nets >= 319 & nets <= 331);
        n61 = find(nets >= 127 & nets <= 148);
        n62 = find(nets >= 332 & nets <= 361);
        n71 = find(nets >= 149 & nets <= 200);
        n72 = find(nets >= 362 & nets <= 400);
        n81 = find(nets >= 401 & nets <= 419);
        NetworkEngagement{1,1}='Visual: ';
        NetworkEngagement{1,2}=sum(length(n11)+length(n12));
        NetworkEngagement{2,1}='SomMotor: ';
        NetworkEngagement{2,2}=sum(length(n21)+length(n22));
        NetworkEngagement{3,1}='DAN: ';
        NetworkEngagement{3,2}=sum(length(n31)+length(n32));
        NetworkEngagement{4,1}='SN: ';
        NetworkEngagement{4,2}=sum(length(n41)+length(n42));
        NetworkEngagement{5,1}='Limbic: ';
        NetworkEngagement{5,2}=sum(length(n51)+length(n52));
        NetworkEngagement{6,1}='Exec: ';
        NetworkEngagement{6,2}=sum(length(n61)+length(n62));
        NetworkEngagement{7,1}='DMN: ';
        NetworkEngagement{7,2}=sum(length(n71)+length(n72));
        NetworkEngagement{8,1}='Subcort: ';
        NetworkEngagement{8,2}=sum(length(n81));
        
    end
    filename=fullfile(savedir, [savefile 'SigMatrix05Storey_SampEn.mat']);
    if exist('Sig05')
        if exist('SigMatrix05')
            filename=fullfile(savedir, [savefile 'SigMatrix05Storey_SampEn.mat']);
            save(filename,'Sig05','SigMatrix05','Allresults','SigUncorrect05','Sig05Storey','BENsubs','Databackup','NetworkEngagement','CohenEffect')
            x=1;
        else
            filename=fullfile(savedir, [savefile 'SigMatrix05Bonf_SampEn.mat']);
            save(filename,'Sig05','Allresults','SigUncorrect05','BENsubs','Databackup','NetworkEngagement','CohenEffect')
            x=1;
        end
    elseif exist('SigMatrix05')
        filename=fullfile(savedir, [savefile 'SigMatrix05Storey_SampEn.mat']);
        save(filename,'SigMatrix05','Allresults','SigUncorrect05','Sig05Storey','BENsubs','Databackup','NetworkEngagement','CohenEffect')
        x=1;
    elseif exist('SigUncorrect05')
        filename=fullfile(savedir, [savefile 'SigMatrixUncorr_SampEn.mat']);
        save(filename,'SigUncorrect05','Allresults','BENsubs','Databackup','NetworkEngagement','CohenEffect')
        x=2;
    else
        x=0;
    end

end