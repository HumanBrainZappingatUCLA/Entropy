clear all
load('BENavgcompare_post.mat')
load('BENavgcompare_pre.mat')

for i = 1:size(BENavgcompare_pre,2)
    if BENavgcompare_pre(i) == BENavgcompare_post(i)
        PrevsPost(i)=0;
    elseif BENavgcompare_pre(i) < BENavgcompare_post(i)
        PrevsPost(i)= 1;
    elseif BENavgcompare_pre(i) > BENavgcompare_post(i)
        PrevsPost(i)= -1;
      
    end
end

imagesc(PrevsPost)
colorbar

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

save('C:\Users\tjjordan\Documents\UCLAResearch\TS_Analysis\Atlas_remapping\PrevsPost_avgBENcomp2.mat','PrevsPost')