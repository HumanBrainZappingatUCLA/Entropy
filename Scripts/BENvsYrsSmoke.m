function x = BENvsYrsSmoke(basedir,groupname,Transform_matBEN)
    load('C:\Users\tjjordan\Documents\UCLAResearch\TS_Analysis\Behavior\Yrs_Smoked.mat')
    x=[];
    subfiledir=dir(fullfile(basedir,'BEN',groupname,'BENoutputs','sub*'));
    filedir=fullfile(basedir,'BEN',groupname,'BENoutputs');
    
    savedir= fullfile(basedir,'BEN',groupname, 'Correlations');
    saveplots=fullfile(savedir,'Plots');
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
        for i=1:size(SENroi,1)
            BEN_array(1,i)=SENroi(i,1);
        end
        if mod(n3,2) == 0
            BEN_After(k2,:) = BEN_array;
            BENgs_After(k2,:)=SENglobal;
            bensubs(k2)=str2num(n2(5:9));
            k2=k2+1;
        else
            BEN_Before(k1,:) = BEN_array;
            BENgs_Before(k1,:)=SENglobal;
            k1=k1+1;
        end
    
        clearvars SENroi SENglobal BEN_array
    
    end
    Yrs_smoke=cell2mat(Yrsmoke);
    k=1;
    for n = 1:2:length(subfiles)
        n2 = convertStringsToChars(subfiles(n));
        n3 = str2num(n2(5:9));
        if ismember(n3,Yrs_smoke(:,1))
            ind=find(Yrs_smoke==n3);
            Smoketime(k,:)= Yrs_smoke(ind(1),:);
            Smokesubs(k)=n3;
            k=k+1;
        end
    end
    
    [sharedval1,ind1]=intersect(bensubs,Smokesubs);
    [sharedval2,ind2]=intersect(Smokesubs,bensubs);
    %load('ROI_TransformMatrix235.mat')
    BENA_corr=BEN_After(ind1,:);
    BENB_corr=BEN_Before(ind1,:);
    Smoker_corr=Smoketime(ind2,2);
    [row,col]=find(isnan(BENB_corr));
    for i = 1:length(row);
        BENB_corr(row(i),col(i))=0;
    end
    % BENB_corr(isnan(BENB_corr))=0;
    % BENB_corr(isinf(BENB_corr))=max(BENB_corr(:))+5;
    % BENA_corr(isnan(BENA_corr))=0;
    % BENA_corr(isinf(BENA_corr))=max(BENA_corr(:))+5;
    for i=1:size(BENB_corr,2)
        [rho,pval] = corr(Smoker_corr,BENB_corr(:,i));%,'type','Spearman');
        mdl=fitlm(Smoker_corr,BENB_corr(:,i));
        rhosq=mdl.Rsquared.Ordinary;
        %p = polyfit(BENB_corr(:,i),Smoker_corr,1);
        stats_Before(i,1)=rho;
        stats_Before(i,2)=pval;
        stats_Before(i,3)=rhosq;
    end
    ind3=find(stats_Before(:,2) < 0.05);
    [fdr] = mafdr(stats_Before(:,2));
    ind4 = find(fdr < 0.05);
    if ~isempty(ind4)
        SigcorrB=zeros(1,size(BENB_corr,2));
        SigcorrB(ind4) = stats_Before(ind4,1);
    end
    
    
    for i=1:size(BENA_corr,2)
        [rho,pval] = corr(Smoker_corr,BENA_corr(:,i));%,'type','Spearman');
        mdl=fitlm(Smoker_corr,BENA_corr(:,i));
        rhosq=mdl.Rsquared.Ordinary;
        %p = polyfit(BENB_corr(:,i),Smoker_corr,1);
        stats_After(i,1)=rho;
        stats_After(i,2)=pval;
        stats_After(i,3)=rhosq;
    end
    
    ind3A=find(stats_After(:,2) < 0.05);
    [fdr2] = mafdr(stats_After(:,2));
    ind4A = find(fdr2 < 0.05);
    if ~isempty(ind4A)
        SigcorrA=zeros(1,size(BENA_corr,2));
        SigcorrA(ind4A) = stats_After(ind4A,1);
        
    end
    for i=1:length(ind4A)
        SigcorrA05(i,1)=ind4A(i);
        SigcorrA05(i,2)=stats_After(ind4A(i),1);
        SigcorrA05(i,3)=fdr2(ind4A(i));
    end
    for i=1:length(ind4)
        SigcorrB05(i,1)=ind4(i);
        SigcorrB05(i,2)=stats_Before(ind4(i),1);
        SigcorrB05(i,3)=fdr(ind4(i));
        Sigcomp(i,1)=ind4(i);
        Sigcomp(i,2)=stats_Before(ind4(i),1);
        Sigcomp(i,3)=stats_After(ind4(i),1);
        Sigcomp(i,4)=stats_Before(ind4(i),1)-stats_After(ind4(i),1);
        Sigcomp(i,5)=stats_Before(ind4(i),3);
        Sigcomp(i,6)=stats_After(ind4(i),3);
        if ~ismember(ind4(i),ind4A)
            Sigcomp(i,7)=1;
        else
            Sigcomp(i,7)=0;
        end
    end
%     SigcompStor05=Sigcomp(Sigcomp(:,7)==1,1:6);
%     SigcompStor05rsq=SigcompStor05(SigcompStor05(:,5)>0.1,1:6);

    for i=1:length(ind3)
        Sigcomp2(i,1)=ind3(i);
        Sigcomp2(i,2)=stats_Before(ind3(i),1);
        Sigcomp2(i,3)=stats_After(ind3(i),1);
        Sigcomp2(i,4)=stats_Before(ind3(i),1)-stats_After(ind3(i),1);
        Sigcomp2(i,5)=stats_Before(ind3(i),3);
        Sigcomp2(i,6)=stats_After(ind3(i),3);
        if ~ismember(ind3(i),ind3A)
            Sigcomp2(i,7)=1;
        else
            Sigcomp2(i,7)=0;
        end
    end
    SigcompUncorr=Sigcomp2(Sigcomp2(:,5)==1,1:4);
    
    filename=fullfile(savedir,['BENvsYrsSmoke_PrevsPost_' groupname '.mat']);
    save(filename,'SigcompUncorr','Sigcomp2')%,'Sigcomp','SigcompStor05rsq')
    
    BENA_plot=BENA_corr(:,Sigcomp2(:,1));
    BENB_plot=BENB_corr(:,Sigcomp2(:,1));
    Smoker_plot=Smoker_corr;
    LW = 2;
    LW2 = 4;
    FS = 42;
    DS  =15;
    for i=1:size(BENA_plot,2)
    %     p1 = polyfit(Smoker_corr,BENB_plot(:,i),1);
    %     p2 = polyfit(Smoker_corr,BENA_plot(:,i),1);
    %     a=p1(1);
    %     b=p1(2);
        rvalB=Sigcomp2(i,2);
        rvalA=Sigcomp2(i,3);
        rvalDiff=Sigcomp2(i,4);
        rsqB=Sigcomp2(i,5);
        rsqA=Sigcomp2(i,6);
        
        Smoker_plot=Smoker_corr;
        BENBp=BENB_plot(:,i);
        BENAp=BENA_plot(:,i);
        if any(isinf(BENBp))
            indrm=find(isinf(BENBp));
            BENBp(indrm)=[];
            Smoker_plot(indrm)=[];
            BENAp(indrm)=[];
        end
        if Transform_matBEN(Sigcomp2(i,1)) == 100
            ROI = 'Insula';
        elseif Transform_matBEN(Sigcomp2(i,1)) == 141
            ROI = 'DLPFC';
        else
            ROI = num2str(Transform_matBEN(Sigcomp2(i,1)));
        end
        g(1,1)=gramm('x',Smoker_plot','y',BENBp');
        g(1,1).set_names('x','Smoking Exposure (yrs)','y','Regional SampEn');
        g(1,1).geom_point();
        g(1,1).stat_glm('geom','line');
        g(1,1).set_title(ROI,'Interpreter','latex');
        g(1,1).set_color_options('map',[0 0 0],'lightness',1);
        g(1,1).no_legend;
        g(1,2)=gramm('x',Smoker_plot','y',BENAp');
        g(1,2).set_names('x','Smoking Exposure (yrs)','y','Regional SampEn');
        g(1,2).geom_point();
        g(1,2).stat_glm('geom','line');
        
        g(1,2).set_title(ROI,'Interpreter','latex');
        g(1,2).set_color_options('map','brewer_dark');
        g(1,2).no_legend;
    %     filename = [outdir2 convertStringsToChars(roi{i}) ' - ' convertStringsToChars(roi{i2}) '_GCvsRT_Scatter.mat'];
    %     save(filename,'g');
        g.draw();
        axes(g(1,1).facet_axes_handles);
        g(1,1).results.geom_point_handle.MarkerSize = DS;
        g(1,1).facet_axes_handles.LineWidth =LW;
        g(1,1).facet_axes_handles.FontSize=FS;
        g(1,1).title_axe_handle.Children.FontSize=FS;
        g(1,1).title_axe_handle.Children.FontWeight='bold';
        g(1,1).results.stat_glm.line_handle.LineWidth=LW2;
        axes(g(1,2).facet_axes_handles);
        g(1,2).results.geom_point_handle.MarkerSize = DS;
        g(1,2).facet_axes_handles.LineWidth =LW;
        g(1,2).facet_axes_handles.FontSize=FS;
        g(1,2).title_axe_handle.Children.FontSize=FS;  
        g(1,2).results.stat_glm.line_handle.LineWidth=LW2;
        f = gcf;
        f.Position = [1 1 1900 1000];
        hold on
    %     x = 1:30;
    %     yB= p1(1)*x + p1(2);
    %     plot(x,yB,'LineWidth',2)
        rstrB = num2str(round(rvalB,2));
        sigstatsB = {['r = ' rstrB]};
    %     yA= p2(1)*x + p2(2);
    %     plot(x,yA,'LineWidth',2)
        rstrA = num2str(round(rvalA,2));
        rstrDiff=num2str(round(rvalDiff,2));
        sigstatsA = {['r = ' rstrA]};
        %sigstatsDiff = {['\Deltar = ' rstrDiff]};
        %annotation('textbox',[.2 .25 .25 .25],'String',sigstatsB,'FitBoxToText','on','Position',[0.434,0.859,0.066,0.047],'FontSize',FS); %Above Legend
        annotation('textbox',[.2 .25 .5 .5],'String',sigstatsB,'FitBoxToText','on','Position',[0.4191,0.8825,0.0804,0.063],'FontSize',30); %Above Legend
        annotation('textbox',[.2 .25 .25 .25],'String',sigstatsA,'FitBoxToText','on','Position',[0.921,0.859,0.066,0.047],'FontSize',18); %Above Legend
        %anot= annotation('textbox',[.2 .25 .25 .25],'String',sigstatsDiff,'FitBoxToText','on','Position',[0.5,0.98,0,0],'FontSize',18); %Above Legend
        %anot.FitBoxToText='on';
        filename=fullfile(saveplots,['ROI' convertStringsToChars(num2str(Transform_matBEN(Sigcomp2(i,1)))) '_PrevsPostCorr.png']);
        saveas(gcf, filename);
        close all
        Plotdata2.(['ROI' convertStringsToChars(ROI) '_vsYrsSmoke']).x=Smoker_plot;
        Plotdata2.(['ROI' convertStringsToChars(ROI) '_vsYrsSmoke']).y=BENBp;
        clear g
        x=1;
    end
    filename=fullfile(savedir,'BENplotdataYrsSmoke.mat');
    save(filename,'Plotdata2')
end