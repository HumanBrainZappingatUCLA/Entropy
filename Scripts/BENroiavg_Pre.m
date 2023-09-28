clear all
bendir='C:\Users\tjjordan\Documents\UCLAResearch\TS_Analysis\Atlas_remapping\BEN\BENoutputs';

getsubs = dir(fullfile(bendir,'*outputs1*'));

for i =1:size(getsubs,1)
    getsubs2(i,1)=convertCharsToStrings(getsubs(i).name);
end

for i = 1:size(getsubs2,1)
    load(fullfile(bendir,getsubs2(i)))
    ROIben(i,:)=SENroi';
end
ROIben2=mean(ROIben,1);


save('C:\Users\tjjordan\Documents\UCLAResearch\TS_Analysis\Atlas_remapping\BENroiavg_pre.mat','ROIben2')