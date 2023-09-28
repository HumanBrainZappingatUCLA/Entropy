function x = BENvsBehavior(basedir,groupname,Transform_matBEN)
    load('C:\Users\tjjordan\Documents\UCLAResearch\TS_Analysis\Behavior\DLPFC_Behaviorstats.mat')
    warning('off','MATLAB:handle_graphics:exceptions:SceneNode')
    x=[];
    subfiledir=dir(fullfile(basedir,'BEN',groupname,'BENoutputs','sub*'));
    filedir=fullfile(basedir,'BEN',groupname,'BENoutputs');
    Behave_names{1, 5}(1) = "Craving";
    
    savedir= fullfile(basedir,'BEN',groupname, 'Correlations');
    saveplots=fullfile(savedir,'Plots','Behavior');
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
        ROInum=size(SENroi,1);
        p=1;
        for i=1:ROInum
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
    
    k=1;
    for n = 1:2:length(subfiles)
        n2 = convertStringsToChars(subfiles(n));
        n3 = str2num(n2(5:9));
        if ismember(n3,Behaviorstats(:,1))
            ind=find(Behaviorstats==n3);
            Behave(k,:)= Behaviorstats(ind(1),:);
            Behavesubs(k)=n3;
            k=k+1;
        end
    end
    
    [sharedval1,ind1]=intersect(bensubs,Behavesubs);
    [sharedval2,ind2]=intersect(Behavesubs,bensubs);
    
    BENA_corr=BEN_After(ind1,:);
    BENB_corr=BEN_Before(ind1,:);
    Behave_corr=Behave(ind2,2:end);
%     [row,col]=find(isnan(BENB_corr));
%     for i = 1:length(row)
%         BENB_corr(row(i),col(i))=[];
%         BENA_corr(row(i),col(i))=[];
% 
%     end
    p=1;
    for i=1:size(BENB_corr,2)
        %p=1;
        for j=1:size(Behave_corr,2)
            BENb=BENB_corr(:,i);
            BENa=BENA_corr(:,i);
            Behave_corr2=Behave_corr(:,j);
            row = find(isnan(BENb));
            row2 = find(isnan(Behave_corr2));
            row3 = [row;row2];
            row3 = unique(row3);
            BENb(row3)=[];
            BENa(row3)=[];
            Behave_corr2(row3)=[];
            BEN_diff=BENa-BENb;
            Behave_corr2(isoutlier(BEN_diff))=[];
            BEN_diff(isoutlier(BEN_diff))=[];
            BEN_diff(isoutlier(Behave_corr2))=[];
            Behave_corr2(isoutlier(Behave_corr2))=[];
            [rho,pval] = corr(BEN_diff,Behave_corr2);%,'type','Spearman');
            mdl=fitlm(BEN_diff,Behave_corr2);
            rhosq=mdl.Rsquared.Ordinary;
            %p = polyfit(BENB_corr(:,i),Smoker_corr,1);
            stats_Before(p,1)=i;
            stats_Before(p,2)=j;
            stats_Before(p,3)=rho;
            stats_Before(p,4)=pval;
            stats_Before(p,5)=rhosq;
            p=p+1;
            
        end
    end
    
    p=1;
    for i=1:size(stats_Before,1)
        if stats_Before(i,3) > 0.2
            Corrstats(p,:)=stats_Before(i,:);
            p=p+1;
        end
    end
    sprd=1;
    for i=1:size(Corrstats,1)
        if Corrstats(i,5) > 0.1
            CorrstatsSprd(sprd,:)=Corrstats(i,:);
            sprd=sprd+1;
        end
    end
    fdr=mafdr(Corrstats(:,4));
    ind05=find(fdr<0.05);
    
    Corrstats05=Corrstats(ind05,:);
    Corrstats05(:,4)=fdr(ind05);
    filename=fullfile(savedir,['BENvsBehavior_' groupname '.mat']);
    save(filename,'Corrstats','Corrstats05','ind05','CorrstatsSprd')
    
    
    
    %% Plot ones correlating with Yrs Smoke
     clearvars -except Corrstats05 BENA_corr BENB_corr Behave_corr ind05 groupname savedir fdr Behave_names saveplots Transform_matBEN
    
    %load('ROI_TransformMatrix235.mat')
    q=1;
    for i=1:size(Corrstats05,1)
        %if Corrstats05(i,2) ==5 || Corrstats05(i,2) ==13 || Corrstats05(i,2) ==7 || Corrstats05(i,2) ==9 || Corrstats05(i,2) ==11 || Corrstats05(i,2) ==15 || Corrstats05(i,2) ==17
        if Corrstats05(i,2) ==5 || Corrstats05(i,2) ==6 || Corrstats05(i,2) ==17% && Corrstats05(i,5) >= 0.3
            Corr_plot(q,:)=Corrstats05(i,:);
            q=q+1;
        end
    end

    LW = 2;
    LW2 = 4;
    FS = 48;
    DS  =15;
    for i=1:size(Corr_plot,1)
        BENA_plot=BENA_corr(:,Corr_plot(i,1));
        BENB_plot=BENB_corr(:,Corr_plot(i,1));
        row = find(isnan(BENB_plot));
        BENB_plot(row)=[];
        BENA_plot(row)=[];
        
        BENchange=BENA_plot-BENB_plot;
        Behave_plot=Behave_corr(:,Corr_plot(i,2));
        Behave_plot(row)=[];
        row2 = find(isnan(Behave_plot));
        Behave_plot(row2)=[];
        BENchange(row2)=[];
        rmoutliers = isoutlier(Behave_plot);
        BENchange(rmoutliers)=[];
        Behave_plot(rmoutliers)=[];
        rval=Corr_plot(i,3);
        if any(isinf(BENchange))
            indrm=find(isinf(BENchange));
            BENchange(indrm)=[];
            Behave_plot(indrm)=[];
        end
        if Transform_matBEN(Corr_plot(i,1)) == 100
            ROI = 'L Insula';
        elseif Transform_matBEN(Corr_plot(i,1)) == 141 || Transform_matBEN(Corr_plot(i,1)) == 137 || Transform_matBEN(Corr_plot(i,1)) == 138
            ROI = 'L dlPFC';
        elseif Transform_matBEN(Corr_plot(i,1)) == 133
            ROI = 'L ITG';
        elseif Transform_matBEN(Corr_plot(i,1)) == 314 || Transform_matBEN(Corr_plot(i,1)) == 318
            ROI = 'R SFG';
        else
            ROI = num2str(Transform_matBEN(Corr_plot(i,1)));
        end
        title=[ROI ' Entropy vs ' convertStringsToChars(Behave_names{1,Corr_plot(i,2)})];
        g(1,1)=gramm('x',BENchange','y',Behave_plot');
        behave_axtitle=['Change in ' convertStringsToChars(Behave_names{1,Corr_plot(i,2)})];
        g(1,1).set_names('x','Change in Sample Entropy','y',behave_axtitle);
        g(1,1).geom_point();
        g(1,1).stat_glm('geom','line');
        g(1,1).set_title(title,'Interpreter','latex');
        g(1,1).set_color_options('map',[0 0 0],'lightness',1);
        g(1,1).no_legend;
        g.draw();
        axes(g(1,1).facet_axes_handles);
        g(1,1).results.geom_point_handle.MarkerSize = DS;
        g(1,1).facet_axes_handles.LineWidth =LW;
        g(1,1).facet_axes_handles.FontSize=FS;
        g(1,1).title_axe_handle.Children.FontSize=FS;
        g(1,1).results.stat_glm.line_handle.LineWidth=LW2;
        f = gcf;
        f.Position = [1 1 1900 1000];
        hold on
        rstr = num2str(round(rval,2));
        sigstats = {['r = ' rstr]};
        annotation('textbox',[0.848421052631579 0.278 0.083684210526316 0.0570000000000003],'String',sigstats,'FitBoxToText','off','FontSize',32); %Above Legend
        %anot.FitBoxToText='on';
        behavename= convertStringsToChars(Behave_names{1,Corr_plot(i,2)});
        behavename= behavename(find(~isspace(behavename)));
        filename=fullfile(saveplots,['ROI' convertStringsToChars(num2str(Transform_matBEN(Corr_plot(i,1)))) '_vsBehave' behavename '.png']);
        saveas(gcf, filename);
        
        close all
        Plotdata2.(['ROI' convertStringsToChars(num2str(Transform_matBEN(Corr_plot(i,1)))) '_vsBehave' behavename]).x=BENchange;
        Plotdata2.(['ROI' convertStringsToChars(num2str(Transform_matBEN(Corr_plot(i,1)))) '_vsBehave' behavename]).y=Behave_plot;
        Plotdata=[num2cell(BENchange),num2cell(Behave_plot)];
        excelname=['ROI' convertStringsToChars(num2str(Transform_matBEN(Corr_plot(i,1)))) '_vsBehave' convertStringsToChars(Behave_names{1,Corr_plot(i,2)})];
        writecell(Plotdata,fullfile(savedir,'BENplotdata.xlsx'),'Sheet',excelname(1:end-7))
        clear g Plotdata
    
        x=1;
    end
    filename=fullfile(savedir,'BENplotdata.mat');
    save(filename,'Plotdata2')
end

