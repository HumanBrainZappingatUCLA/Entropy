clear all
basedir='C:\Users\tjjordan\Documents\UCLAResearch\TS_Analysis\Behavior';
load(fullfile(basedir,'Behavior_Results.mat'))
sessnum = inputdlg('How many Conditions? 1 - 4','',[1 35]);
sessnum=str2num(sessnum{1});
if sessnum > 4
    disp('Error: Too Many Conditions given')
    return
elseif sessnum < 1
    disp('Error: No conditions given')
    return
end
condgroup=[];
c=1;
while c <= sessnum
    condlist= {'DLPFC','SFG','PPC','v5'};
    condselect = listdlg('PromptString', {'Please select 1 condition?'}, ...
    'ListString',condlist,'SelectionMode','single');
    if isempty(condselect)
        return
    end
    if condselect == 1
        if exist(fullfile(basedir,"DLPFC_Behaviorstats.mat"))
            load(fullfile(basedir,"DLPFC_Behaviorstats.mat"))
            condgroup=[condgroup condselect];
            c=c+1;
        else
            disp("DLPFC Behavior data doesn't exist")
            disp("Pick another Condition")
        end
    elseif condselect == 2
        if exist(fullfile(basedir,"SFG_Behaviorstats.mat"))
            load(fullfile(basedir,"SFG_Behaviorstats.mat"))
            condgroup=[condgroup condselect];
            c=c+1;
        else
            disp("SFG Behavior data doesn't exist")
            disp("Pick another Condition")
        end
    elseif condselect == 3
        if exist(fullfile(basedir,"PPC_Behaviorstats.mat"))
            load(fullfile(basedir,"PPC_Behaviorstats.mat"))
            condgroup=[condgroup condselect];
            c=c+1;
        else
            disp("PPC Behavior data doesn't exist")
            disp("Pick another Condition")
        end
    elseif condselect == 4
        if exist(fullfile(basedir,"v5_Behaviorstats.mat"))
            load(fullfile(basedir,"v5_Behaviorstats.mat"))
            condgroup=[condgroup condselect];
            c=c+1;
        else
            disp("v5 Behavior data doesn't exist")
            disp("Pick another Condition")
        end
    end
    
end
clear Behaviorstats

subs1=DLPFC_behaveChange(:,1);
subs1=subs1';
subs2=v5_behaveChange(:,1);
subs2=subs2';
[subs,ind1,ind2] = intersect(subs1,subs2);
dlpfcchange= DLPFC_behaveChange(ind1,:);
v5change= v5_behaveChange(ind2,:);