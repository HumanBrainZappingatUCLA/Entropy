clear all
basedir='C:\Users\tjjordan\Documents\UCLAResearch\TS_Analysis';
groupname='AvgAllPrevsDLPFCPost419';
load('C:\Users\tjjordan\Documents\UCLAResearch\TS_Analysis\Behavior\SubJ_Ethn.mat')

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
SubEthn=cell2mat(SubEthn);
k=1;
for n = 1:2:length(subfiles)
    n2 = convertStringsToChars(subfiles(n));
    n3 = str2num(n2(5:9));
    if ismember(n3,SubEthn(:,1))
        ind=find(SubEthn==n3);
        Smoketime(k,:)= SubEthn(ind(1),:);
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
    p= anova1(BENB_corr(:,i),Smoker_corr,'off');
    Allresults(i,1)=i;
    Allresults(i,2)=p;
    if p < 0.05/419
        Sig05(corr,:)=Allresults(i,:);
        corr=corr+1;
    elseif p < 0.05
        Sig05uncorr(uncorr,:)=Allresults(i,:);
        uncorr=uncorr+1;

    end
    %close all
end

 [fdr] = mafdr(Allresults(:,2));
 indc=find(fdr < 0.05);

g1=1;
g2=1;
g3=1;
g4=1;
g5=1;
g6=1;
g7=1;
g8=1;

for j = 1:size(Smoker_corr,1)
    if Smoker_corr(j) == 1
        AIANben(g1,:)=BENB_corr(j,:);
        g1=g1+1;
    elseif Smoker_corr(j) == 2
        Asianben(g2,:)=BENB_corr(j,:);
        g2=g2+1;
    elseif Smoker_corr(j) == 3
        NHOPben(g3,:)=BENB_corr(j,:);
        g3=g3+1;
    elseif Smoker_corr(j) == 4
        BAAben(g4,:)=BENB_corr(j,:);
        g4=g4+1;
    elseif Smoker_corr(j) == 5
        Whiteben(g5,:)=BENB_corr(j,:);
        g5=g5+1;
    elseif Smoker_corr(j) == 6
        HISPben(g6,:)=BENB_corr(j,:);
        g6=g6+1;
    elseif Smoker_corr(j) == 7
        MTORben(g7,:)=BENB_corr(j,:);
        g7=g7+1;
    elseif Smoker_corr(j) == 8
        UKNben(g8,:)=BENB_corr(j,:);
        g8=g8+1;
    
    end
end
AllgroupBENm=[];
AllgroupBENstd=[];
Whichgroup=[];
if exist("AIANben",'var')
    AllgroupBENm=[ AllgroupBENm ;mean(AIANben,1)];
    AllgroupBENstd=[ AllgroupBENstd ;std(AIANben,1)];
    Whichgroup=[Whichgroup,1];
end
if exist("Asianben",'var')
    AllgroupBENm=[ AllgroupBENm ;mean(Asianben,1)];
    AllgroupBENstd=[ AllgroupBENstd ;std(Asianben,1)];
    Whichgroup=[Whichgroup,2];
end
if exist("NHOPben",'var')
    AllgroupBENm=[ AllgroupBENm ;mean(NHOPben,1)];
    AllgroupBENstd=[ AllgroupBENstd ;zeros(1,419)];
    Whichgroup=[Whichgroup,3];
end
if exist("BAAben",'var')
    AllgroupBENm=[ AllgroupBENm ;mean(BAAben,1)];
    AllgroupBENstd=[ AllgroupBENstd ;std(BAAben,1)];
    Whichgroup=[Whichgroup,4];
end
if exist("Whiteben",'var')
    AllgroupBENm=[ AllgroupBENm ;mean(Whiteben,1)];
    AllgroupBENstd=[ AllgroupBENstd ;std(Whiteben,1)];
    Whichgroup=[Whichgroup,5];
end
if exist("HISPben",'var')
    AllgroupBENm=[ AllgroupBENm ;mean(HISPben,1)];
    AllgroupBENstd=[ AllgroupBENstd ;std(HISPben,1)];
    Whichgroup=[Whichgroup,6];
end
if exist("MTORben",'var')
    AllgroupBENm=[ AllgroupBENm ;mean(MTORben,1)];
    AllgroupBENstd=[ AllgroupBENstd ;std(MTORben,1)];
    Whichgroup=[Whichgroup,7];
end
if exist("UKNben",'var')
    AllgroupBENm=[ AllgroupBENm ;mean(UKNben,1)];
    AllgroupBENstd=[ AllgroupBENstd ;std(UKNben,1)];
    Whichgroup=[Whichgroup,8];
end

nodelist=[ 35 98 99 100 143 234 235 236 302 303 304 305 340 137 138 139 140 141 142];
Nodecell=cell(length(nodelist),size(AllgroupBENm,1)+1);
Nodecell{1,1}="Node Number";
for i = 1:size(Whichgroup,2)
    if Whichgroup(i) == 1
        Nodecell{1,i+1}=convertCharsToStrings(['American Indian/Alaska Native (N=' num2str(size(AIANben,1)) ')']);
    elseif Whichgroup(i) == 2
        Nodecell{1,i+1}=convertCharsToStrings(['Asian (N=' num2str(size(Asianben,1)) ')']);
    elseif Whichgroup(i) == 3
        Nodecell{1,i+1}=convertCharsToStrings(['Native Hawaiian or Other Pacific Islander (N=' num2str(size(NHOPben,1)) ')']);
    elseif Whichgroup(i) == 4
        Nodecell{1,i+1}=convertCharsToStrings(['Black or African American (N=' num2str(size(BAAben,1)) ')']);
    elseif Whichgroup(i) == 5
        Nodecell{1,i+1}=convertCharsToStrings(['White (N=' num2str(size(Whiteben,1)) ')']);
    elseif Whichgroup(i) == 6
        Nodecell{1,i+1}=convertCharsToStrings(['Hispanic (N=' num2str(size(HISPben,1)) ')']);
    elseif Whichgroup(i) == 7
        Nodecell{1,i+1}=convertCharsToStrings(['More than One Race (N=' num2str(size(MTORben,1)) ')']);
    elseif Whichgroup(i) == 8
        Nodecell{1,i+1}=convertCharsToStrings(['Unknown or Prefer Not to Specify (N=' num2str(size(UKNben,1)) ')']);
    end
end
for i = 1:length(nodelist)
    Nodecell{i+1,1}=nodelist(i);
    for j = 1:size(AllgroupBENm,1)
        Nodecell{i+1,j+1}=[num2str(round(AllgroupBENm(j,nodelist(i)),1)) '(' num2str(round(AllgroupBENstd(j,nodelist(i)),1)) ')'];
    end
end

writecell(Nodecell,'C:\Users\tjjordan\Documents\UCLAResearch\TS_Analysis\Behavior\EthnDiff.xlsx')