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
data_pro_right_old = {};
data_pro_right_yng = {};
data_anti_left_old = {};
data_anti_left_yng = {};
data_pro_left_old = {};
data_pro_left_yng = {};
data_anti_right_old = {};
data_anti_right_yng = {};

data_pro_right_mpe_old = {};
data_pro_right_mpe_yng = {};
data_anti_left_mpe_old = {};
data_anti_left_mpe_yng = {};
data_pro_left_mpe_old = {};
data_pro_left_mpe_yng = {};
data_anti_right_mpe_old = {};
data_anti_right_mpe_yng = {};

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
        
            
        
        countblocks = 0;
        previous = '';
        for e = 1:length(FullEEG.event)
            if contains (FullEEG.event(e).type,'94')
                countblocks = countblocks+1;
            end
            if countblocks == 2 || countblocks == 3 || countblocks == 4
                if contains (FullEEG.event(e).type,'10')
                    FullEEG.event(e).type = '12 ';
                elseif contains (FullEEG.event(e).type,'11')
                    FullEEG.event(e).type = '13 ';
               
                end
                if contains (FullEEG.event(e).type,'40') 
                    FullEEG.event(e).type = '41 ';
                end
            end

            if strcmp(FullEEG.event(e).type, 'L_saccade')
                if contains(previous, '10 ')
                    FullEEG.event(e).type = 'L_saccade_10';
                elseif contains(previous, '11 ')
                    FullEEG.event(e).type = 'L_saccade_11';
                elseif contains(previous, '12 ')
                    FullEEG.event(e).type = 'L_saccade_12';
                elseif contains(previous, '13 ')
                    FullEEG.event(e).type = 'L_saccade_13';
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
        
         %%%%przed dzieleniem na epochi
         % eegfilt(13-30)
        FullEEG.data = eegfilt(FullEEG.data, 500, 13, 30);


        BEGIN = -0.5;
        END = 0.5;
        %data epoching

        data_pro_right_mpe{i} = my_pop_epoch(FullEEG, '11  ', 'L_saccade_11');
        data_pro_left_mpe{i} = my_pop_epoch(FullEEG, '10  ', 'L_saccade_10');
        data_anti_right_mpe{i} = my_pop_epoch(FullEEG, '13  ', 'L_saccade_13');
        data_anti_left_mpe{i} = my_pop_epoch(FullEEG, '12  ', 'L_saccade_12');
        
        data_pro_right_mpe{i} = TransformData(data_pro_right_mpe{i});
        data_pro_left_mpe{i} = TransformData(data_pro_left_mpe{i});
        data_anti_right_mpe{i} = TransformData(data_anti_right_mpe{i});
        data_anti_left_mpe{i} = TransformData(data_anti_left_mpe{i});
        
         %dividing - old/young
        eval(['data_pro_right_mpe_' OLD_OR_YOUNG{young+1} '{end+1} = data_pro_right_mpe{i}']);
        eval(['data_pro_left_mpe_' OLD_OR_YOUNG{young+1} '{end+1} = data_pro_left_mpe{i}'])
        eval(['data_anti_right_mpe_' OLD_OR_YOUNG{young+1} '{end+1} = data_anti_right_mpe{i}']);
        eval(['data_anti_left_mpe_' OLD_OR_YOUNG{young+1} '{end+1} = data_anti_left_mpe{i}']);
        
    end
        
end

%%
for i=1:length(data_pro_right_mpe_yng)
    
    [data_pro_right_mpe_yng{i}, pri_cell] =clip_tb_cell(data_pro_right_mpe_yng{i}, -100, 500);
    [data_pro_left_mpe_yng{i}, pli_cell] = clip_tb_cell(data_pro_left_mpe_yng{i}, -500, 100);
    [data_anti_right_mpe_yng{i}, ari_cell] = clip_tb_cell(data_anti_right_mpe_yng{i}, -500, 100);
    [data_anti_left_mpe_yng{i}, ali_cell] = clip_tb_cell(data_anti_left_mpe_yng{i}, -100, 500);
    
end

for i=1:length(data_pro_right_mpe_old)
    
    [data_pro_right_mpe_old{i}, pri_cell] =clip_tb_cell(data_pro_right_mpe_old{i}, -100, 500);
    [data_pro_left_mpe_old{i}, pli_cell] = clip_tb_cell(data_pro_left_mpe_old{i}, -500, 100);
    [data_anti_right_mpe_old{i}, ari_cell] = clip_tb_cell(data_anti_right_mpe_old{i}, -500, 100);
    [data_anti_left_mpe_old{i}, ali_cell] = clip_tb_cell(data_anti_left_mpe_old{i}, -100, 500);
    
end
%within subjects, over all conditions
for i=1:length(data_pro_right_mpe_old)
    data_pro_right_mpe_old{i} = Hilbobo(data_pro_right_mpe_old{i});
    data_pro_left_mpe_old{i} = Hilbobo(data_pro_left_mpe_old{i});
    data_anti_right_mpe_old{i} = Hilbobo(data_anti_right_mpe_old{i});
    data_anti_left_mpe_old{i} = Hilbobo(data_anti_left_mpe_old{i});
end
for i=1:length(data_pro_right_mpe_yng)
    data_pro_right_mpe_yng{i} = Hilbobo(data_pro_right_mpe_yng{i});
    data_pro_left_mpe_yng{i} = Hilbobo(data_pro_left_mpe_yng{i});
    data_anti_right_mpe_yng{i} = Hilbobo(data_anti_right_mpe_yng{i});
    data_anti_left_mpe_yng{i} = Hilbobo(data_anti_left_mpe_yng{i});
end

