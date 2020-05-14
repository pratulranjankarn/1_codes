clc;
close all;
clear all;
fid = fopen('delayfile_final.txt','r');
A = textscan(fid,'%s %f %f %f %s %f %f %f %f %f %f','headerlines',1);
fclose(fid);
n = find(strcmp(A{5},'APE'));
orig_time = A{1}(n);
T0 = [2000 01 01 00 00 00];
T_i = [2011 01 01 00 00 00];
T_f = [2012 04 30 23 59 59];
j = 1;
k = 1;
l = 1;
m = 1;
for i=1:length(orig_time)
    if strcmp(orig_time{i,1}(19),'.')
        sec = strcat('0',orig_time{i,1}(18));
    else
        sec = orig_time{i,1}(18:19);
    end
    year = str2double(orig_time{i,1}(1:4));
    mnth = str2double(orig_time{i,1}(6:7));
    date = str2double(orig_time{i,1}(9:10));
    hwr = str2double(orig_time{i,1}(12:13));
    mnt = str2double(orig_time{i,1}(15:16));
    scnd = str2double(sec);
    T = [year mnth date hwr mnt scnd];
    if etime(T,T0) < etime(T_i,T0)
        orig_b_E_orig{j,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_E.sacii');
        orig_b_E_24{j,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_E_2-4.sacii');
        orig_b_E_48{j,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_E_4-8.sacii');
        orig_b_N_orig{j,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_N.sacii');
        orig_b_N_24{j,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_N_2-4.sacii');
        orig_b_N_48{j,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_N_4-8.sacii');
        if isnan(A{10}(n(i)))
            delay_bf(j,:) = [string(orig_time{i,1}) A{8}(n(i)) A{9}(n(i)) string(0)];
        else
            orig_b_E_816{m,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_E_8-16.sacii');
        
            orig_b_N_816{m,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_N_8-16.sacii');
            m = m + 1;
        
            delay_bf(j,:) = [string(orig_time{i,1}) A{8}(n(i)) A{9}(n(i)) A{10}(n(i))];
        end
        j = j + 1;
    elseif etime(T,T0) > etime(T_f,T0)
        orig_a_E_orig{k,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_E.sacii');
        orig_a_E_24{k,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_E_2-4.sacii');
        orig_a_E_48{k,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_E_4-8.sacii');
        orig_a_E_816{k,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_E_8-16.sacii');
        orig_a_N_orig{k,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_N.sacii');
        orig_a_N_24{k,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_N_2-4.sacii');
        orig_a_N_48{k,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_N_4-8.sacii');
        orig_a_N_816{k,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_N_8-16.sacii');
        delay_af(k,:) = [string(orig_time{i,1}) A{8}(n(i)) A{9}(n(i)) A{10}(n(i))];
        k = k + 1;
    else
        orig_c_E_orig{l,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_E.sacii');
        orig_c_E_24{l,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_E_2-4.sacii');
        orig_c_E_48{l,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_E_4-8.sacii');
        orig_c_E_816{l,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_E_8-16.sacii');
        orig_c_N_orig{l,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_N.sacii');
        orig_c_N_24{l,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_N_2-4.sacii');
        orig_c_N_48{l,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_N_4-8.sacii');
        orig_c_N_816{l,1} = strcat('/home/user/COMB_EV/Event_wav_final/',...
            orig_time{i,1}(1:10),'H',orig_time{i,1}(12:13),'M',...
            orig_time{i,1}(15:16),'S',sec,'/N_APE_N_8-16.sacii');
        delay_cn(l,:) = [string(orig_time{i,1}) A{8}(n(i)) A{9}(n(i)) A{10}(n(i))];
        l = l + 1;
    end
end
orig_b_E_orig = string(orig_b_E_orig);
orig_b_E_24 = string(orig_b_E_24);
orig_b_E_48 = string(orig_b_E_48);
orig_b_E_816 = string(orig_b_E_816);
orig_b_N_orig = string(orig_b_N_orig);
orig_b_N_24 = string(orig_b_N_24);
orig_b_N_48 = string(orig_b_N_48);
orig_b_N_816 = string(orig_b_N_816);
orig_a_E_orig = string(orig_a_E_orig);
orig_a_E_24 = string(orig_a_E_24);
orig_a_E_48 = string(orig_a_E_48);
orig_a_E_816 = string(orig_a_E_816);
orig_a_N_orig = string(orig_a_N_orig);
orig_a_N_24 = string(orig_a_N_24);
orig_a_N_48 = string(orig_a_N_48);
orig_a_N_816 = string(orig_a_N_816);
orig_c_E_orig = string(orig_c_E_orig);
orig_c_E_24 = string(orig_c_E_24);
orig_c_E_48 = string(orig_c_E_48);
orig_c_E_816 = string(orig_c_E_816);
orig_c_N_orig = string(orig_c_N_orig);
orig_c_N_24 = string(orig_c_N_24);
orig_c_N_48 = string(orig_c_N_48);
orig_c_N_816 = string(orig_c_N_816);

fid = fopen('delay_APE_bf','w');
fprintf(fid,'%s %s %s %s\n',delay_bf');
fclose(fid);

fid = fopen('delay_APE_af','w');
fprintf(fid,'%s %s %s %s\n',delay_af');
fclose(fid);

fid = fopen('delay_APE_cn','w');
fprintf(fid,'%s %s %s %s\n',delay_cn');
fclose(fid);

fid = fopen('bf_crscr_E_orig','w');
fprintf(fid,'%s\n',orig_b_E_orig);
fclose(fid);

fid = fopen('bf_crscr_E_2-4','w');
fprintf(fid,'%s\n',orig_b_E_24);
fclose(fid);

fid = fopen('bf_crscr_E_4-8','w');
fprintf(fid,'%s\n',orig_b_E_48);
fclose(fid);

fid = fopen('bf_crscr_E_8-16','w');
fprintf(fid,'%s\n',orig_b_E_816);
fclose(fid);

fid = fopen('bf_crscr_N_orig','w');
fprintf(fid,'%s\n',orig_b_N_orig);
fclose(fid);

fid = fopen('bf_crscr_N_2-4','w');
fprintf(fid,'%s\n',orig_b_N_24);
fclose(fid);

fid = fopen('bf_crscr_N_4-8','w');
fprintf(fid,'%s\n',orig_b_N_48);
fclose(fid);

fid = fopen('bf_crscr_N_8-16','w');
fprintf(fid,'%s\n',orig_b_N_816);
fclose(fid);

fid = fopen('af_crscr_E_orig','w');
fprintf(fid,'%s\n',orig_a_E_orig);
fclose(fid);

fid = fopen('af_crscr_E_2-4','w');
fprintf(fid,'%s\n',orig_a_E_24);
fclose(fid);

fid = fopen('af_crscr_E_4-8','w');
fprintf(fid,'%s\n',orig_a_E_48);
fclose(fid);

fid = fopen('af_crscr_E_8-16','w');
fprintf(fid,'%s\n',orig_a_E_816);
fclose(fid);

fid = fopen('af_crscr_N_orig','w');
fprintf(fid,'%s\n',orig_a_N_orig);
fclose(fid);

fid = fopen('af_crscr_N_2-4','w');
fprintf(fid,'%s\n',orig_a_N_24);
fclose(fid);

fid = fopen('af_crscr_N_4-8','w');
fprintf(fid,'%s\n',orig_a_N_48);
fclose(fid);

fid = fopen('af_crscr_N_8-16','w');
fprintf(fid,'%s\n',orig_a_N_816);
fclose(fid);

fid = fopen('cn_crscr_E_orig','w');
fprintf(fid,'%s\n',orig_c_E_orig);
fclose(fid);

fid = fopen('cn_crscr_E_2-4','w');
fprintf(fid,'%s\n',orig_c_E_24);
fclose(fid);

fid = fopen('cn_crscr_E_4-8','w');
fprintf(fid,'%s\n',orig_c_E_48);
fclose(fid);

fid = fopen('cn_crscr_E_8-16','w');
fprintf(fid,'%s\n',orig_c_E_816);
fclose(fid);

fid = fopen('cn_crscr_N_orig','w');
fprintf(fid,'%s\n',orig_c_N_orig);
fclose(fid);

fid = fopen('cn_crscr_N_2-4','w');
fprintf(fid,'%s\n',orig_c_N_24);
fclose(fid);

fid = fopen('cn_crscr_N_4-8','w');
fprintf(fid,'%s\n',orig_c_N_48);
fclose(fid);

fid = fopen('cn_crscr_N_8-16','w');
fprintf(fid,'%s\n',orig_c_N_816);
fclose(fid);