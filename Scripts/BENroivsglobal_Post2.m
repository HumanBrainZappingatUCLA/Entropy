clear all
bendir=fullfile(cd,'..','BEN','G1vsG2','BENoutputs');

getsubs = dir(fullfile(bendir,'*outputs2*'));

for i =1:size(getsubs,1)
    getsubs2(i,1)=convertCharsToStrings(getsubs(i).name);
end

for i = 1:size(getsubs2,1)
    load(fullfile(bendir,getsubs2(i)))
    ROIben(i,:)=SENroi';
    GSben(i,1)=mean(SENroi);
end

load('BENavgcompare_pre.mat', 'GSben')

for i = 1:size(ROIben,2)
    GSdata=GSben(:);
    ROIdata=ROIben(:,i);
    [h,p,ci,stats] = ttest(ROIdata,GSdata);
    Allresults(i,1)=h;
    Allresults(i,2)=p;
    Allresults(i,3)=stats.tstat;
    Allresults(i,4)=mean(GSdata);
    Allresults(i,5)=std(GSdata);
    Allresults(i,6)=mean(ROIdata);
    Allresults(i,7)=std(ROIdata);
    Databackup{1,i}=i;
    Databackup{2,i}=[GSdata ROIdata];
end
[fdr,Q] = mafdr(Allresults(:,2));
ind = find(fdr < 0.05);

for i =1:size(Allresults,1)
    if ismember(i,ind)
        if Allresults(i,3) < 0
            AboveBelowmat(i)=-1;
        elseif Allresults(i,3) > 0
            AboveBelowmat(i)=1;
        else
            AboveBelowmat(i)=0;
        end
    else
        AboveBelowmat(i)=0;
    end
end
load('ROI_TransformMatrix419.mat')
AboveBelowmat2=AboveBelowmat(Transform_mat);
figure
imagesc(AboveBelowmat2)
p0=0;
p1=61.5;
p2=138.5;
p3=184.5;
p4=231.5;
p5=257.5;
p6=309.5;
p7=400.5;
p8=419.5;
LW=1;
line([p1,p1],[p0,p1],'Linewidth',LW,'color','red')
line([p2,p2],[p0,p1],'Linewidth',LW,'color','red')
line([p3,p3],[p0,p1],'Linewidth',LW,'color','red')
line([p4,p4],[p0,p1],'Linewidth',LW,'color','red')
line([p5,p5],[p0,p1],'Linewidth',LW,'color','red')
line([p6,p6],[p0,p1],'Linewidth',LW,'color','red')
line([p7,p7],[p0,p1],'Linewidth',LW,'color','red')
s0=0;
s1=61.5;
s2=138.5;
s3=184.5;
s4=231.5;
s5=257.5;
s6=309.5;
s7=400.5;
s8=419.5;
ticks =[(s0+(s1-s0)/2) (s1+(s2-s1)/2) (s2+(s3-s2)/2) (s3+(s4-s3)/2) (s4+(s5-s4)/2) (s5+(s6-s5)/2) (s6+(s7-s6)/2) (s7+(s8-s7)/2)];
ticklabels = {'Vis','SomMotor','DorAttn','Sal','Lim','Exec','DMN','SubCort'};
xticks(ticks)
xticklabels(ticklabels)

BENavgcompare_post2=AboveBelowmat2;
BENavgcompare_post=AboveBelowmat;
save('BENavgcompare_post.mat','BENavgcompare_post','Allresults','BENavgcompare_post2','Databackup')