% =========================================================================
% This script create a new experiment folder and save itself into that 
% folder 
% =========================================================================

% get dates
date_string = datetime_str;
year_ = year(meas.DateNumber);
month_ = month(meas.DateNumber);
day_ = day(meas.DateNumber);
% generate the experiment folder name
% exp_folder_name = [datestr(datetime, 'yyyy_mm_dd__HH_MM_SS')];
exp_folder_name = datestr(meas.DateNumber, 'yyyy_mm_dd__HH_MM_SS');

folder_path = fullfile(save_path, num2str(year_), num2str(month_), num2str(day_), exp_folder_name);

% create folder (and all folders above)
[~, ~, ~]   = mkdir(folder_path);
[~, ~, ~]  =  mkdir(fullfile(folder_path, 'plots'));

% copy the whole folder DorB_operational to local D hardrive
print_wp_to_txt;
copyfile(main_path, folder_path)

% delete obsolete folder and temporary files
if isfolder(fullfile(folder_path, 'obsolete'))
    warning('off','MATLAB:DELETE:DirectoryDeletion');
    delete(fullfile(folder_path, 'obsolete','*'));
    try
        rmdir(fullfile(folder_path, 'obsolete'));
    end
    warning('on','MATLAB:DELETE:DirectoryDeletion');
end

tmpfiles = [dir([folder_path , '\*.asv']) ; dir([folder_path , '\*\*.asv'])];
for ind=1:length(tmpfiles)
   delete(fullfile(tmpfiles(ind).folder,tmpfiles(ind).name))
end

% copy the whole folder minime_operational to git folder
% git_path = 'C:\Users\ASAFMAN\Documents\GitHub\minime_operational';
% copyfile(main_path, git_path)