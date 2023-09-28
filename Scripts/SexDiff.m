clear all
load('C:\Users\tjjordan\Documents\UCLAResearch\TS_Analysis\BEN\AvgAllPrevsDLPFCPost419\Correlations\AvgAllPrevsDLPFCPost419_BENchange.mat')
load('C:\Users\tjjordan\Documents\UCLAResearch\Behavior_MasterSheet\SmokingCleanSheet.mat')
for i = 1:size(SmokingCleanSheet,1)
    sublist(i,1)=convertCharsToStrings(SmokingCleanSheet{i,1});
end
for i = 1:size(BENsubs,1)
    s1 = convertStringsToChars(BENsubs(i));
    subs(i)= convertCharsToStrings(s1(5:end));
end
corr=1;
uncorr=1;
ROInum=size(BENchange,2);
for k =1:ROInum

n=1;
m=1;
for i =1:length(BENsubs)-1
    s1=subs(i);
    ind= find(sublist == s1);
    ind1=ind(1);
    sex=convertCharsToStrings(SmokingCleanSheet{ind1,2});
    if sex == "M"
        Mvals(n)= BENchange(i,k);
        n=n+1;
    elseif sex == "F"
        Fvals(m)= BENchange(i,k);
        m=m+1;
    end
    
end

[h,p,ci,stats] = ttest2(Mvals',Fvals');

Allresults(k,1)=h;
Allresults(k,2)=p;
Allresults(k,3)=stats.tstat;
Allresults(k,4)=mean(Mvals);
Allresults(k,5)=std(Mvals);
Allresults(k,6)=mean(Fvals);
Allresults(k,7)=std(Fvals);
if p < 0.05
    SigUncorrect05(uncorr,1)=k;
    SigUncorrect05(uncorr,2)=h;
    SigUncorrect05(uncorr,3)=p;
    SigUncorrect05(uncorr,4)=stats.tstat;
    SigUncorrect05(uncorr,5)=mean(Mvals);
    SigUncorrect05(uncorr,6)=std(Mvals);
    SigUncorrect05(uncorr,7)=mean(Fvals);
    SigUncorrect05(uncorr,8)=std(Fvals);
    uncorr=uncorr+1;
end
if p < (0.05/ROInum)
    Sig05(corr,1)=k;
    Sig05(corr,2)=h;
    Sig05(corr,3)=p;
    Sig05(corr,4)=stats.tstat;
    Sig05(corr,5)=mean(Mvals);
    Sig05(corr,6)=std(Mvals);
    Sig05(corr,7)=mean(Fvals);
    Sig05(corr,8)=std(Fvals);
    corr=corr+1;
end
clearvars -except ROInum Allresults Sig05 subs sublist SmokingCleanSheet ...
    BENchange BENsubs uncorr corr SigUncorrect05
end