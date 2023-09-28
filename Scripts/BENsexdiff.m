clear all
basedir='C:\Users\tjjordan\Documents\UCLAResearch\TS_Analysis';
groupname='AvgAllPrevsDLPFCPost419';
load('C:\Users\tjjordan\Documents\UCLAResearch\TS_Analysis\Behavior\SubJ_Sex.mat')

subfiledir=dir(fullfile(basedir,'BEN',groupname,'BENoutputs','sub*'));
filedir=fullfile(basedir,'BEN',groupname,'BENoutputs');

savedir= fullfile(basedir,'BEN',groupname, 'Correlations');
saveplots=fullfile(savedir,'Plots','Edu');
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
SubSex=cell2mat(SubSex);
k=1;
for n = 1:2:length(subfiles)
    n2 = convertStringsToChars(subfiles(n));
    n3 = str2num(n2(5:9));
    if ismember(n3,SubSex(:,1))
        ind=find(SubSex==n3);
        Smoketime(k,:)= SubSex(ind(1),:);
        Smokesubs(k)=n3;
        k=k+1;
    end
end

[sharedval1,ind1]=intersect(bensubs,Smokesubs);
[sharedval2,ind2]=intersect(Smokesubs,bensubs);

BENB_corr=BEN_Before(ind1,:);
Smoker_corr=Smoketime(ind2,2);
uncorr=1;
corr=1;
for i=1:size(BENB_corr,2)
    n=1;
    m=1;
    for j=1:size(BENB_corr,1)
        if Smoker_corr(j)==1
            BENm(m)=BENB_corr(j,i);
            m=m+1;
        elseif Smoker_corr(j)==2
            BENf(n)=BENB_corr(j,i);
            n=n+1;
        end
    end
    [h,p,ci,stats]=ttest2(BENm,BENf);
    Allresults(i,1)=i;
    Allresults(i,2)=h;
    Allresults(i,3)=p;
    Allresults(i,4)=stats.tstat;
    Allresults(i,5)=mean(BENm);
    Allresults(i,6)=mean(BENf);
    if p < 0.05/419
        Sig05(corr,:)=Allresults(i,:);
        corr=corr+1;
    elseif p < 0.05
        Sig05uncorr(uncorr,:)=Allresults(i,:);
        uncorr=uncorr+1;

    end
    
end

 [fdr] = mafdr(Allresults(:,3));
 indc=find(fdr < 0.05);

 