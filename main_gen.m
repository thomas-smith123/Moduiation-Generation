clear all;
close all;

%% defination
fs = 64e9;
resample_fs = 5e9;
max_signal_per_frame = 4;
resolution = [512,512];
window_length = resolution(2);
overlap = 128;
path = './data_gen';
snr_ = [-10,-5,0,5,10,15];
number = 600;
check = true;
thread = 4;
%% gen
parfor i=1:thread
    gen = keysight_signal_gen(fs,resample_fs,max_signal_per_frame,...
        resolution,window_length,overlap,char(strcat(path,'_',string(i),'/')),snr_,number,check);
    gen.call();
end
%%
