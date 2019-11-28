%the proper pipeline for EEG analysis

%paths
clear all
addpath('/Users/mplome/Desktop/eeglab14_1_2b')
addpath('/Users/mplome/Desktop')
addpath('tools');
eeglab;
close
addpath('/Users/mplome/Desktop/Tools')
table_old = readtable("/Users/mplome/Downloads/old.xlsx");
table_young =readtable("/Users/mplome/Downloads/young.xlsx");
oldIDsA = table_old{3:end,1};
oldIDsB= table_old{3:end,2};
all_oldIDs = [oldIDsA; oldIDsB];
youngIDsA = table_young{3:end,1};
youngIDsB= table_young{3:end,2};
all_youngIDs = [youngIDsA; youngIDsB];

raw= '/Users/mplome/Desktop/balanced/EEG_preprocessed'; % path to preprocessed eeg files
etfolder='/Users/mplome/Desktop/balanced/ET'; %MO

d=dir(raw) %what folders are in there (each folder = one subject)

d(1:3)=[] % get rid of the . and .. folders as well as .DS_Store on mac

OLD_OR_YOUNG = {'old', 'yng'};
data_pro_old = {};
data_pro_yng = {};
data_anti_old = {};
data_anti_yng = {};


data_pro_mpe_old = {};
data_pro_mpe_yng = {};
data_anti_mpe_old = {};
data_anti_mpe_yng = {};

for i=1:length(d) %loop over all subjects
    if d(i).isdir
        subjectfolder=dir([d(i).folder filesep d(i).name  ]);
        
        deleteindex=[];
        for ii=1:length(subjectfolder)
            if not(endsWith(subjectfolder(ii).name, '_EEG.mat')) || startsWith(subjectfolder(ii).name,'bip') || startsWith(subjectfolder(ii).name,'red')
                deleteindex(end+1)=ii;
            end
        end
        subjectfolder(deleteindex)=[];
        FullEEG=[];
        for ii=1:length(subjectfolder)
            load ([subjectfolder(ii).folder filesep subjectfolder(ii).name]) % gets loaded as EEG
            fileindex=subjectfolder(ii).name(end-8) %here you need to find the index from thefile (end-someting) indexing
            etfile=  [etfolder filesep d(i).name filesep d(i).name '_AS' fileindex '_ET.mat'] %define string of the complete path to the matching ET file.
             
            %EEG = pop_reref(EEG,[47 83])
            EEG = pop_reref(EEG,[]) %tu jest wazna sprawa, to reref
            
            
            %merge ET into EEG
            ev1=94 %first trigger of eeg and ET file
            ev2=50 % end trigger in eeg and ET file
            EEG=pop_importeyetracker(EEG, etfile,[ev1 ev2], [1:4], {'TIME' 'L_GAZE_X' 'L_GAZE_Y' 'L_AREA'},1,1,0,0,4)
            %sprawdz pupil size, w razie czego zmien
            if ii==1
                FullEEG=EEG;
            else
                FullEEG=pop_mergeset(FullEEG,EEG);
            end
        end

        if isempty(FullEEG)
            continue
        end
        
            
        previous = '';
        for e = 1:length(FullEEG.event)
            if strcmp(FullEEG.event(e).type, 'L_saccade')
                if contains(previous, '10 ')
                    FullEEG.event(e).type = 'L_saccade_10';
                elseif contains(previous, '11 ')
                    FullEEG.event(e).type = 'L_saccade_11';
                end
            end
            
            if ~strcmp(FullEEG.event(e).type, 'L_fixation') ...
                    && ~strcmp(FullEEG.event(e).type, 'L_blink')
                previous = FullEEG.event(e).type;
            end
        end
        
        

        
        id=d(i).name ;
        young=    any(contains(all_youngIDs,id));
        old =     any(contains(all_oldIDs,id));
        
   
        %young means 1, old means 0
        all_ages(i) =    young;
        all_eeg{i} = FullEEG;
        
        
        BEGIN = -0.5;
        END = 0.5;
        %data epoching
        data_pro{i} = pop_epoch(FullEEG, {'11  '}, [BEGIN,END]);
        data_anti{i} = pop_epoch(FullEEG, {'10  '} , [BEGIN,END]);
  
        data_pro{i} = pop_rmbase(data_pro{i},1000*[BEGIN,0]);
        data_pro{i} = pop_rmbase(data_pro{i},1000*[BEGIN,0]);
        
         %I'm interested only in saccades with latencies [100;500]ms ale nie
        %dziala
        data_pro{i} = FilterByEventLatency(data_pro{i}, 'L_saccade_11', 100, 500);
        data_anti{i} = FilterByEventLatency(data_anti{i}, 'L_saccade_10', 100, 500);
        
        %dividing - old/young
        eval(['data_pro_' OLD_OR_YOUNG{young+1} '{end+1} = data_pro{i}']);
        eval(['data_anti_' OLD_OR_YOUNG{young+1} '{end+1} = data_anti{i}']);
        
        data_pro_mpe{i} = my_pop_epoch(FullEEG, '11  ', 'L_saccade_11');
        data_anti_mpe{i} = my_pop_epoch(FullEEG, '10  ', 'L_saccade_10');
        
        data_pro_mpe{i} = TransformData(data_pro_mpe{i});
        data_anti_mpe{i} = TransformData(data_anti_mpe{i});
        
         %dividing - old/young
        eval(['data_pro_mpe_' OLD_OR_YOUNG{young+1} '{end+1} = data_pro_mpe{i}']);
        eval(['data_anti_mpe_' OLD_OR_YOUNG{young+1} '{end+1} = data_anti_mpe{i}']);
        
    end
        
end
%%
csvsz = [length(data_pro), 4];%znow, one moga miec rozne dlugosci
vartypes = {'string', 'logical' 'double', 'double'};
varnames = {'id', 'age', 'PRO', 'ANTI'};
csvdata = table('Size', csvsz, 'VariableTypes', vartypes, 'VariableNames', varnames);
for i=1:length(data_pro)%dlaczego, przeciez one moga miec rozne dlugosci...
    id=d(i).name;
    
    VAR_NAMES = {
        'data_pro'
        'data_anti'
    };
    ITEMS = {
        {data_pro{i}, 'L_saccade_11', 'PRO ', [-100,500]}
       
        {data_anti{i}, 'L_saccade_10', 'ANTI ', [-100,500]}
    };
    I_DATA = 1;
    I_TRIGGER = 2;
    I_TITLE = 3;
    I_LIM = 4;
    
    csvdata(i, 1) = {id};
    csvdata(i, 2) = {all_ages(i)};

    figure;
    sgtitle(id);
    for I=1:length(ITEMS)
        X = ITEMS{I}{I_DATA};
        vmin = ITEMS{I}{I_LIM}(1);
        vmax = ITEMS{I}{I_LIM}(2);
        [Y,ind] = ClipToBounds(X, vmin, vmax);
        
        color = char(ones(length(ind),1) * 'r');
        color(ind) = 'k';

        subplot(4,1,I)
        for k=1:size(X.data,3)
            hold on;
            plot(X.times, X.data(107,:,k), color(k));
            xlim([-500,500])
        end
        title(ITEMS{I}{I_TITLE});
        
        csvdata(i, I + 2) = {sum(ind)/length(ind)};
    end
    path = '/Users/mplome/Desktop/EEG_ANALYSIS/good_and_bad/';
    savefig([path, id '.fig']);
    saveas(gcf,[path, id , '.png']);
    writetable(csvdata, [path 'summary.csv']); 
end
%%
for i=1:length(data_pro_yng)
    [data_pro_yng{i}, pi] =ClipToBounds(data_pro_yng{i}, -100, 500);
    [data_anti_yng{i}, ai] = ClipToBounds(data_anti_yng{i}, -500, 100);
    
    [data_pro_mpe_yng{i}, pi_cell] = clip_tb_cell(data_pro_mpe_yng{i}, -500, 100);
    [data_anti_mpe_yng{i}, ai_cell] = clip_tb_cell(data_anti_mpe_yng{i}, -100, 500);
    
end

for i=1:length(data_pro_old)
    [data_pro_old{i}, pi] =ClipToBounds(data_pro_old{i}, -100, 500);
    [data_anti_old{i}, ai] = ClipToBounds(data_anti_old{i}, -100, 500);
    
    [data_pro_mpe_old{i}, pi_cell] =clip_tb_cell(data_pro_mpe_old{i}, -100, 500);
    [data_anti_mpe_old{i}, ai_cell] = clip_tb_cell(data_anti_mpe_old{i}, -100, 500);
    
%     data_pro_old{i} = Slice(data_pro_old{i}, pri);
%     data_pro_old{i} = Slice(data_pro_old{i}, pli);
%     data_anti_old{i} = Slice(data_anti_old{i}, ari);
%     data_anti_old{i} = Slice(data_anti_old{i}, ali);
% %     
%     data_pro_mpe_old{i} = SliceEpoch_cell(data_pro_mpe_old{i}, pri_cell);
%     data_pro_mpe_old{i} = SliceEpoch_cell(data_pro_mpe_old{i}, pli_cell);
%     data_anti_mpe_old{i} = SliceEpoch_cell(data_anti_mpe_old{i}, ari_cell);
%     data_anti_mpe_old{i} = SliceEpoch_cell(data_anti_mpe_old{i}, ali_cell);
end


% TU WYKRESY
%%
% young

%ANTI
min_len = inf;
for i=1:length(data_anti_mpe_yng)
    min_len = min(min_len, data_anti_mpe_yng{i}.min_trial);
end

tab_anti_yng = zeros(data_anti_mpe_yng{1}.nbchan, min_len, 0);
for i=1:length(data_anti_mpe_yng)%bierze wszysykich mlodych
    tab_anti_yng = cat(3, tab_anti_yng,data_anti_mpe_yng{i}.suffix_data(:,end - min_len +1:end,:));%po 3 wymiarzxe laczy- 3 wymiar to sa triale
    erp_anti_yng(i,:,:) = mean(data_anti_mpe_yng{i}.suffix_data(:,end - min_len +1:end,:),3);
    % Zamien "1:min_len" na "min_len+1:end"
end


%pro
min_len = inf;
for i=1:length(data_pro_mpe_yng)
    min_len = min(min_len, data_pro_mpe_yng{i}.min_trial);
end

tab_pro_yng = zeros(data_pro_mpe_yng{1}.nbchan, min_len, 0);
for i=1:length(data_pro_mpe_yng)%bierze wszysykich mlodych
    tab_pro_yng = cat(3, tab_pro_yng,data_pro_mpe_yng{i}.suffix_data(:,end - min_len +1:end,:));%po 3 wymiarzxe laczy- 3 wymiar to sa triale
    erp_pro_yng(i,:,:) = mean(data_pro_mpe_yng{i}.suffix_data(:,end - min_len +1:end,:),3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%old
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%
%anti 

min_len = inf;
for i=1:length(data_anti_mpe_old)
    min_len = min(min_len, data_anti_mpe_old{i}.min_trial);
end

tab_anti_old = zeros(data_anti_mpe_old{1}.nbchan, min_len, 0);
for i=1:length(data_anti_mpe_old)%bierze wszysykich mlodych
    tab_anti_old = cat(3, tab_anti_old,data_anti_mpe_old{i}.suffix_data(:,end - min_len +1:end,:));%po 3 wymiarzxe laczy- 3 wymiar to sa triale
    erp_anti_old(i,:,:) = mean(data_anti_mpe_old{i}.suffix_data(:,end - min_len +1:end,:),3);
end


%pro
min_len = inf;
for i=1:length(data_pro_mpe_old)
    min_len = min(min_len, data_pro_mpe_old{i}.min_trial);
end

tab_pro_old= zeros(data_pro_mpe_old{1}.nbchan, min_len, 0);
for i=1:length(data_pro_mpe_old)%bierze wszysykich mlodych
   tab_pro_old= cat(3, tab_pro_old,data_pro_mpe_old{i}.suffix_data(:,end - min_len +1:end,:));%po 3 wymiarzxe laczy- 3 wymiar to sa triale
   erp_pro_old(i,:,:) = mean(data_pro_mpe_old{i}.suffix_data(:,end - min_len +1:end,:),3);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%% średnie po kanałach
mean_channel_activity_pro_yng = mean(mean(tab_pro_yng([9 94 104], :, :),3),1);

mean_channel_activity_pro_old = mean(mean(tab_pro_old([9 94 104], :, :),3),1);

mean_channel_activity_anti_yng = mean(mean(tab_anti_yng([9 94 104], :, :),3),1);

mean_channel_activity_anti_old = mean(mean(tab_anti_old([9 94 104], :, :),3),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first do the loop over the subjects, then calcultate erp for each of
nsubj_old=size(data_anti_mpe_old,2)
nsubj_yng=size(data_anti_mpe_yng,2)

%%%%%%%%YOUNG GROUP, group averages
groupaverage_erp_pro_yng = squeeze(mean(mean(erp_pro_yng(:,[19,9,104],:),2),1));
shaded_erp_pro_yng = shadedErrorBar(1:size(erp_pro_yng,3), groupaverage_erp_pro_yng,std(mean(erp_pro_yng(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'r');

groupaverage_erp_anti_yng = squeeze(mean(mean(erp_anti_yng(:,[19,9,104],:),2),1));
shaded_erp_anti_yng = shadedErrorBar(1:size(erp_anti_yng,3), groupaverage_erp_anti_yng,std(mean(erp_anti_yng(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'r');
%%%%%%%%%%%%%%%%%%OLD PEOPLE, group averages
groupaverage_erp_pro_old = squeeze(mean(mean(erp_pro_old(:,[19,9,104],:),2),1));
shadedErrorBar(1:size(erp_pro_old,3), groupaverage_erp_pro_old,std(mean(erp_pro_old(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'r');

groupaverage_erp_anti_old = squeeze(mean(mean(erp_anti_old(:,[19,9,104],:),2),1));
shadedErrorBar(1:size(erp_anti_old,3), groupaverage_erp_anti_old,std(mean(erp_anti_old(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%PLOTS%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
figure
subplot(2,1,1)
shadedErrorBar(1:size(erp_pro_old,3), groupaverage_erp_pro_old,std(mean(erp_pro_old(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'r');
hold on
shadedErrorBar(1:size(erp_pro_yng,3), groupaverage_erp_pro_yng,std(mean(erp_pro_yng(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'b');
legend({'old', 'yng'})
title ('PRO')
hold off

subplot(2,1,2)
shadedErrorBar(1:size(erp_anti_old,3), groupaverage_erp_anti_old,std(mean(erp_anti_old(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'r');
hold on
shadedErrorBar(1:size(erp_anti_yng,3), groupaverage_erp_anti_yng,std(mean(erp_anti_yng(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'b');
legend({'old', 'yng'})
title ('ANTI')
hold off
%%

subplot(2,2,4)
shadedErrorBar(1:size(erp_anti_old,3), groupaverage_erp_anti_old,std(mean(erp_anti_old(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'r');
hold on
shadedErrorBar(1:size(erp_anti_yng,3), groupaverage_erp_anti_yng,std(mean(erp_anti_yng(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'b');
legend({'old', 'yng'})
title ('ANTI LEFT')
hold off


subplot(2,2,2)
shadedErrorBar(1:size(erp_pro_old,3), groupaverage_erp_pro_old,std(mean(erp_pro_old(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'r');
hold on
shadedErrorBar(1:size(erp_pro_yng,3), groupaverage_erp_pro_yng,std(mean(erp_pro_yng(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'b');
legend({'old', 'yng'})
title ('PRO LEFT')
hold off
sgtitle('Old Group vs Young Group')

%%
%%%


%%
figure
subplot(2,2,1)
shadedErrorBar(1:size(erp_pro_yng,3), groupaverage_erp_pro_yng,std(mean(erp_pro_yng(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'r');
hold on
shadedErrorBar(1:size(erp_anti_yng,3), groupaverage_erp_anti_yng,std(mean(erp_anti_yng(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'b');
legend({'PRO RIGHT', 'ANTI LEFT'})
title ('PRO RIGHT vs ANTI LEFT (yng)')
hold off

subplot(2,2,3)
shadedErrorBar(1:size(erp_anti_yng,3), groupaverage_erp_anti_yng,std(mean(erp_anti_yng(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'r');
hold on
shadedErrorBar(1:size(erp_pro_yng,3), groupaverage_erp_pro_yng,std(mean(erp_pro_yng(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'b');
legend({'ANTI RIGHT', 'PRO LEFT'})
title ('ANTI RIGHT vs PRO LEFT (yng)')
hold off

subplot(2,2,2)
shadedErrorBar(1:size(erp_pro_old,3), groupaverage_erp_pro_old,std(mean(erp_pro_old(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'r');
hold on
shadedErrorBar(1:size(erp_anti_old,3), groupaverage_erp_anti_old,std(mean(erp_anti_old(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'b');
legend({'PRO RIGHT', 'ANTI LEFT'})
title ('PRO RIGHT vs ANTI LEFT (old)')
hold off

subplot(2,2,4)
shadedErrorBar(1:size(erp_anti_old,3), groupaverage_erp_anti_old,std(mean(erp_anti_old(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'r');
hold on
shadedErrorBar(1:size(erp_pro_old,3), groupaverage_erp_pro_old,std(mean(erp_pro_old(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'b');
legend({'ANTI RIGHT', 'PRO LEFT'})
title ('ANTI RIGHT vs PRO LEFT (old)')
hold off

sgtitle('PRO VS ANTI')



%%
%%%%BINNED

%%%przede wszytkim upewnij sie czy faktycznie binujemy na zakresie takim
%%%jaki chce nicolas


%deklaracja wariabli
for i=1:length(data_pro_yng)
    data_pro_mpe_yng{i} = BinnedData(data_pro_mpe_yng{i}, 10);%we want to have 10bins
    data_pro_mpe_yng{i} = BinnedData(data_pro_mpe_yng{i}, 10);
    data_anti_mpe_yng{i} = BinnedData(data_anti_mpe_yng{i}, 10);
    data_anti_mpe_yng{i} = BinnedData(data_anti_mpe_yng{i}, 10);
end 
for i=1:length(data_pro_old)
    data_pro_mpe_old{i} = BinnedData(data_pro_mpe_old{i}, 10);%we want to have 10bins
    data_pro_mpe_old{i} = BinnedData(data_pro_mpe_old{i}, 10);
    data_anti_mpe_old{i} = BinnedData(data_anti_mpe_old{i}, 10);
    data_anti_mpe_old{i} = BinnedData(data_anti_mpe_old{i}, 10);
end    
        
      
%YNG


%pro 
for i=1:length(data_pro_mpe_yng)%bierze wszysykich mlodych
    erp_pro_yng_bin(i,:,:) = mean(data_pro_mpe_yng{i}.binned_data, 3);
end
nsubj_yng=size(data_anti_mpe_yng,2);
groupaverage_erp_pro_yng_bin = squeeze(mean(mean(erp_pro_yng_bin(:,[19,9,104],:),2),1));
shaded_erp_pro_yng_bin = shadedErrorBar(1:size(erp_pro_yng_bin,3), groupaverage_erp_pro_yng_bin,std(mean(erp_pro_yng_bin(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'r');



%anti 
for i=1:length(data_anti_mpe_yng)%bierze wszysykich mlodych
    erp_anti_yng_bin(i,:,:) = mean(data_anti_mpe_yng{i}.binned_data, 3);
end
nsubj_yng=size(data_anti_mpe_yng,2);
groupaverage_erp_anti_yng_bin = squeeze(mean(mean(erp_anti_yng_bin(:,[19,9,104],:),2),1));
shaded_erp_anti_yng_bin = shadedErrorBar(1:size(erp_anti_yng_bin,3), groupaverage_erp_anti_yng_bin,std(mean(erp_anti_yng_bin(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'r');

%OLD

%pro 
for i=1:length(data_pro_mpe_old)%bierze wszysykich mlodych
    erp_pro_old_bin(i,:,:) = mean(data_pro_mpe_old{i}.binned_data, 3);
end
nsubj_old=size(data_anti_mpe_old,2);
groupaverage_erp_pro_old_bin = squeeze(mean(mean(erp_pro_old_bin(:,[19,9,104],:),2),1));
shaded_erp_pro_old_bin = shadedErrorBar(1:size(erp_pro_old_bin,3), groupaverage_erp_pro_old_bin,std(mean(erp_pro_old_bin(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'r');



%anti 
for i=1:length(data_anti_mpe_old)%bierze wszysykich mlodych
    erp_anti_old_bin(i,:,:) = mean(data_anti_mpe_old{i}.binned_data, 3);
end
nsubj_old=size(data_anti_mpe_old,2);
groupaverage_erp_anti_old_bin = squeeze(mean(mean(erp_anti_old_bin(:,[19,9,104],:),2),1));
shaded_erp_anti_old_bin = shadedErrorBar(1:size(erp_anti_old_bin,3), groupaverage_erp_anti_old_bin,std(mean(erp_anti_old_bin(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'r');

%%
%binned - plots
figure
subplot(2,2,1)
shaded_erp_pro_old_bin = shadedErrorBar(1:size(erp_pro_old_bin,3), groupaverage_erp_pro_old_bin,std(mean(erp_pro_old_bin(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'b');
shaded_erp_pro_yng_bin = shadedErrorBar(1:size(erp_pro_yng_bin,3), groupaverage_erp_pro_yng_bin,std(mean(erp_pro_yng_bin(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'c');
legend({'old', 'yng'})
title ('PRO RIGHT')
 ylim([-0.7,1.4])
hold off



subplot(2,2,2)
shaded_erp_pro_old_bin = shadedErrorBar(1:size(erp_pro_old_bin,3), groupaverage_erp_pro_old_bin,std(mean(erp_pro_old_bin(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'b');
shaded_erp_pro_yng_bin = shadedErrorBar(1:size(erp_pro_yng_bin,3), groupaverage_erp_pro_yng_bin,std(mean(erp_pro_yng_bin(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'c');
legend({'old', 'yng'})
title ('PRO LEFT')
 ylim([-0.7,1.4])
hold off



subplot(2,2,3)
shaded_erp_anti_old_bin = shadedErrorBar(1:size(erp_anti_old_bin,3), groupaverage_erp_anti_old_bin,std(mean(erp_anti_old_bin(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'r');
shaded_erp_anti_yng_bin = shadedErrorBar(1:size(erp_anti_yng_bin,3), groupaverage_erp_anti_yng_bin,std(mean(erp_anti_yng_bin(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'm');
legend({'old', 'yng'})
title ('ANTI LEFT')
 ylim([-0.7,1.4])
hold off


subplot(2,2,4)
shaded_erp_anti_old_bin = shadedErrorBar(1:size(erp_anti_old_bin,3), groupaverage_erp_anti_old_bin,std(mean(erp_anti_old_bin(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'r');
shaded_erp_anti_yng_bin = shadedErrorBar(1:size(erp_anti_yng_bin,3), groupaverage_erp_anti_yng_bin,std(mean(erp_anti_yng_bin(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'm');
legend({'old', 'yng'})
title ('ANTI RIGHT')
 ylim([-0.7,1.4])
sgtitle('10 bins')
%nowa figura, mlodzi
figure
subplot(2,2,1)
shaded_erp_pro_yng_bin = shadedErrorBar(1:size(erp_pro_yng_bin,3), groupaverage_erp_pro_yng_bin,std(mean(erp_pro_yng_bin(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'r');
shaded_erp_anti_yng_bin = shadedErrorBar(1:size(erp_anti_yng_bin,3), groupaverage_erp_anti_yng_bin,std(mean(erp_anti_yng_bin(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'm');
legend({'pro right', 'anti left'})
title ('ANTI LEFT  vs PRO RIGHT yng')
 ylim([-0.7,1.4])
 hold off


subplot(2,2,2)
shaded_erp_pro_yng_bin = shadedErrorBar(1:size(erp_pro_yng_bin,3), groupaverage_erp_pro_yng_bin ,std(mean(erp_pro_yng_bin (:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'r');
shaded_erp_anti_yng_bin = shadedErrorBar(1:size(erp_anti_yng_bin,3), groupaverage_erp_anti_yng_bin,std(mean(erp_anti_yng_bin(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'm');
legend({'pro left', 'anti right'})
title ('ANTI RIGHT VS. PRO LEFT yng')
 ylim([-0.7,1.4])
hold off

%nowa fgura, staarzy
%dodaj ylims! bo sa inne so far
subplot(2,2,3)
shaded_erp_pro_old_bin = shadedErrorBar(1:size(erp_pro_old_bin,3), groupaverage_erp_pro_old_bin,std(mean(erp_pro_old_bin(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'r');
shaded_erp_anti_old_bin = shadedErrorBar(1:size(erp_anti_old_bin,3), groupaverage_erp_anti_old_bin,std(mean(erp_anti_old_bin(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'm');
legend({'pro right', 'anti left'})
title ('ANTI LEFT  vs PRO RIGHT old')
 ylim([-0.7,1.4])
hold off


subplot(2,2,4)
shaded_erp_pro_old_bin = shadedErrorBar(1:size(erp_pro_old_bin,3), groupaverage_erp_pro_old_bin ,std(mean(erp_pro_old_bin (:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'r');
shaded_erp_anti_old_bin = shadedErrorBar(1:size(erp_anti_old_bin,3), groupaverage_erp_anti_old_bin,std(mean(erp_anti_old_bin(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'm');
legend({'pro left', 'anti right'})
title ('ANTI RIGHT VS. PRO LEFT old')
 ylim([-0.7,1.4])
hold off
sgtitle('10 bins')
% Koniec %%


%%
%Sinle trial analysis
ITEMS = {
    {data_pro_yng, 'L_saccade_11', 'PRO  YNG'}
    {data_anti_yng, 'L_saccade_10', 'ANTI  YNG'}
    
	{data_pro_old, 'L_saccade_11', 'PRO OLD'}
	{data_anti_old, 'L_saccade_10', 'ANTI OLD'}
	
};
I_DATA = 1;
I_TRIGGER = 2;
I_TITLE = 3;

for I=1:length(ITEMS)
    figure;
    %new method starts here
    X = MergeSets(ITEMS{I}{I_DATA}{:});%
  
    X = SortByLatencyOf(X, ITEMS{I}{I_TRIGGER})%L_saccade_10,11,12,13
    X.data_sorted;
    image_sc_data = squeeze(X.data_sorted(19,:,:))';
    %imagesc(image_sc_data);
    colormap('jet');
    caxis([-25,25]);


    new_to_image  = [];
    rts = [];
    for i = 26: size(image_sc_data,1) - 25
       new_to_image(end + 1, :) =  mean(image_sc_data(i-25:i+25,:),1);
       %if isempty(X.epoch_sorted(i).latency)
       %    rts(end + 1) =0;
       %else
       rts(end + 1) = X.epoch_sorted(i).latency/2;
       %end
    end


    %adding rhe reaction time plot
    imagesc(new_to_image(75:end,151:end))
    hold on
    plot( rts,1:length(rts),'black','LineWidth',3)
    title (ITEMS{I}{I_TITLE})
end
    
function [X] = MergeSets(X, varargin)
    varmat = cell2mat(varargin);
    X.event = horzcat(X.event, varmat.event);
    X.trials = sum(vertcat(X.trials, varmat.trials));
    X.xmin = min(vertcat(X.xmin, varmat.xmin));
    X.xmax = max(vertcat(X.xmax, varmat.xmax));
    X.data = cat(3, X.data, varmat.data);
    X.epoch = horzcat(X.epoch, varmat.epoch);
end


%here sorting
function [dataset_out] = SortByLatencyOf(dataset, event)
    dataset_out = EpochLatencyOf(dataset, event);
    [B,I] = sort([dataset_out.epoch.latency]);
    dataset_out.epoch_sorted = dataset_out.epoch(I);
    dataset_out.data_sorted = dataset_out.data(:,:,I);
end

function [dataset] = EpochLatencyOf(dataset, event)
    for i=1:length(dataset.epoch)
        sacc_idx = strcmp(dataset.epoch(i).eventtype, event);
        if sum(sacc_idx) == 1
            dataset.epoch(i).latency = dataset.epoch(i).eventlatency{sacc_idx};
        else
            dataset.epoch(i).latency = -1;
        end
    end
end


