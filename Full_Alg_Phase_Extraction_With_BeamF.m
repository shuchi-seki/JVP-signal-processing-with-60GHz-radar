clear all;
close all;
clc;

load('BGT60TR13Crecord.mat')



Frames = BGT60TR13Crecord;


Radar_Parameter.Sampling_Frequency_kHz = 1000;
Radar_Parameter.Lower_RF_Frequency_kHz = 5.8e+07;
Radar_Parameter.Upper_RF_Frequency_kHz = 6.35e+07;
Radar_Parameter.Chirp_Time_sec = 0.000132987;
Radar_Parameter.Samples_per_Chirp = 128;
Radar_Parameter.Chirps_per_Frame = 8;
Radar_Parameter.Samples_per_Frame = 3072;
Radar_Parameter.Frame_Period_sec = 0.0040644;


%% Singal Processing
c = 3e8; % Speed of light (m/s)
fc = (Radar_Parameter.Lower_RF_Frequency_kHz+Radar_Parameter.Upper_RF_Frequency_kHz)/2*1000; % Center frequency (Hz)
CRR = 1/Radar_Parameter.Chirp_Time_sec; % Chirp repetition rate (Hz)
FRR=1/Radar_Parameter.Frame_Period_sec;% Frame repetition rate (Hz)
BW = (Radar_Parameter.Upper_RF_Frequency_kHz-Radar_Parameter.Lower_RF_Frequency_kHz)*1000; % Bandwidth (Hz)
gamma = BW*CRR; % Chirp rate (Hz/s)
range_res = c/(2*BW);
max_range = range_res*fix(Radar_Parameter.Sampling_Frequency_kHz*1e3/CRR)/2;
tau = 1/FRR; % Slow time (s)

Vec_alpha_degree = -90:90;
Vec_alpha = pi*Vec_alpha_degree/180; % DoA interval

Samps_p_Frame = Radar_Parameter.Samples_per_Frame;
Chips_p_Frame = Radar_Parameter.Chirps_per_Frame;
Samps_p_Chirp = Radar_Parameter.Samples_per_Chirp;

Nb_Chirp_Mean_FFT = input("Please Enter Number of used Chirps per Frame:\n"); % <= Chips_p_Frame

NumRXAntenna = 3;

%% Range-FFT Processing

raw_data_matrix= Frames(1:floor( length(Frames)/(Samps_p_Frame+1) ) * (Samps_p_Frame+1)); % Take multiple frames with length Samps_p_Frame+1 (NaN value at the end of each frame)
raw_data_matrix_SubFrame = reshape(raw_data_matrix,Samps_p_Frame+1,[]); 
raw_data_matrix_SubFrame = raw_data_matrix_SubFrame(1:end-1,:); % Create raw data matrix of size (Nb of Samps_p_Frame)x(Nb of Frames)

raw_data_matrix_SubFrame_Ant1 = raw_data_matrix_SubFrame(1:(Samps_p_Chirp*Chips_p_Frame),:); % Extracting raw data from each Antenna
raw_data_matrix_SubFrame_Ant2 = raw_data_matrix_SubFrame((Samps_p_Chirp*Chips_p_Frame)+1:2*(Samps_p_Chirp*Chips_p_Frame),:);
raw_data_matrix_SubFrame_Ant3 = raw_data_matrix_SubFrame(2*(Samps_p_Chirp*Chips_p_Frame)+1:end,:);

%raw_data_matrix_SubFrame(isnan(raw_data_matrix_SubFrame)) = 0;

raw_data_matrix_SubFrame_Ant1 = raw_data_matrix_SubFrame_Ant1 - mean(raw_data_matrix_SubFrame_Ant1); % Centering data
raw_data_matrix_SubFrame_Ant2 = raw_data_matrix_SubFrame_Ant2 - mean(raw_data_matrix_SubFrame_Ant2);
raw_data_matrix_SubFrame_Ant3 = raw_data_matrix_SubFrame_Ant3 - mean(raw_data_matrix_SubFrame_Ant3);


I_Ant1 = [];
Q_Ant1 = [];

I_Ant2 = [];
Q_Ant2 = [];

I_Ant3 = [];
Q_Ant3 = [];

temp_index = [];


i_start = 1 + input(['Please Enter the Index of Target Bin in FFT-Range [<= ', num2str((Samps_p_Chirp/2)-1) ,'] : \n']);

h = waitbar(0,'Please wait...');
for i_seg = 1:length(raw_data_matrix_SubFrame_Ant1(1,:))

    zero_padding_factor = 1; % For FFT

    %=====> Antenna 1 Data

    raw_data_matrix_chirp_seg_Ant1 = raw_data_matrix_SubFrame_Ant1(:,i_seg);
    raw_data_matrix_chirp_seg_Ant1 = reshape(raw_data_matrix_chirp_seg_Ant1,Samps_p_Chirp*zero_padding_factor,[]);
    w = window(@blackmanharris,Samps_p_Chirp*zero_padding_factor);
    raw_data_matrix_chirp_seg_Win_Ant1 = raw_data_matrix_chirp_seg_Ant1.*w;
    raw_data_matrix_chirp_seg_Win_FFT_Ant1 = fft(raw_data_matrix_chirp_seg_Win_Ant1,Samps_p_Chirp*zero_padding_factor);
    
%     figure
%     plot(abs(raw_data_matrix_chirp_seg_Win_FFT_Ant1))
%     pause

    %raw_data_matrix_chirp_seg_Win_FFT_Ant1_mean = mean(raw_data_matrix_chirp_seg_Win_FFT_Ant1.');
    if Nb_Chirp_Mean_FFT ==1
        raw_data_matrix_chirp_seg_Win_FFT_Ant1_mean = raw_data_matrix_chirp_seg_Win_FFT_Ant1(:,Chips_p_Frame/2).';
    elseif Nb_Chirp_Mean_FFT ~= 1
        raw_data_matrix_chirp_seg_Win_FFT_Ant1_mean = mean(raw_data_matrix_chirp_seg_Win_FFT_Ant1(:,1:Nb_Chirp_Mean_FFT).');
    end


    filter_pic_Ant1 = zeros(1,Samps_p_Chirp*zero_padding_factor);
%     i_start = 1 + input("Please Enter the Index of Bin [1 or 2]: \n");
    temp_Ant1 = abs(raw_data_matrix_chirp_seg_Win_FFT_Ant1_mean(1,i_start:Samps_p_Chirp/2));
    [Val_Ant1, Max_index_Ant1] = max(temp_Ant1);
    %Max_index_Ant1
    Max_index_Ant1 = 1;
    filter_pic_Ant1(1,Max_index_Ant1+i_start-1) = 1;

    %temp_index = [temp_index, Val_Ant1];

    raw_data_matrix_chirp_seg_Win_FFT_mean_Filtered_Ant1 = raw_data_matrix_chirp_seg_Win_FFT_Ant1_mean .* filter_pic_Ant1;


    I_signal_Ant1 = real(ifft(raw_data_matrix_chirp_seg_Win_FFT_mean_Filtered_Ant1));
    Q_signal_Ant1 = imag(ifft(raw_data_matrix_chirp_seg_Win_FFT_mean_Filtered_Ant1));

    I_Ant1 = [I_Ant1; I_signal_Ant1.'];
    Q_Ant1 = [Q_Ant1; Q_signal_Ant1.'];
    
    %======================================================================


    %=====> Antenna 2 Data

    raw_data_matrix_chirp_seg_Ant2 = raw_data_matrix_SubFrame_Ant2(:,i_seg);
    raw_data_matrix_chirp_seg_Ant2 = reshape(raw_data_matrix_chirp_seg_Ant2,Samps_p_Chirp*zero_padding_factor,[]);
    w = window(@blackmanharris,Samps_p_Chirp*zero_padding_factor);
    raw_data_matrix_chirp_seg_Win_Ant2 = raw_data_matrix_chirp_seg_Ant2.*w;
    raw_data_matrix_chirp_seg_Win_FFT_Ant2 = fft(raw_data_matrix_chirp_seg_Win_Ant2,Samps_p_Chirp*zero_padding_factor);
    
%     figure
%     plot(abs(raw_data_matrix_chirp_seg_Win_FFT_Ant2))
%     pause

    %raw_data_matrix_chirp_seg_Win_FFT_Ant2_mean = mean(raw_data_matrix_chirp_seg_Win_FFT_Ant2.');
    if Nb_Chirp_Mean_FFT ==1
        raw_data_matrix_chirp_seg_Win_FFT_Ant2_mean = raw_data_matrix_chirp_seg_Win_FFT_Ant2(:,Chips_p_Frame/2).';
    elseif Nb_Chirp_Mean_FFT ~= 1
        raw_data_matrix_chirp_seg_Win_FFT_Ant2_mean = mean(raw_data_matrix_chirp_seg_Win_FFT_Ant2(:,1:Nb_Chirp_Mean_FFT).');
    end

    filter_pic_Ant2 = zeros(1,Samps_p_Chirp*zero_padding_factor);
    %i_start = 2;%*zero_padding_factor;
    temp_Ant2 = abs(raw_data_matrix_chirp_seg_Win_FFT_Ant2_mean(1,i_start:Samps_p_Chirp/2));
    [Val_Ant2, Max_index_Ant2] = max(temp_Ant2);
    %Max_index_Ant2
    Max_index_Ant2 = 1;
    filter_pic_Ant2(1,Max_index_Ant2+i_start-1) = 1;

    %temp_index = [temp_index, Val_Ant2];

    raw_data_matrix_chirp_seg_Win_FFT_mean_Filtered_Ant2 = raw_data_matrix_chirp_seg_Win_FFT_Ant2_mean .* filter_pic_Ant2;


    I_signal_Ant2 = real(ifft(raw_data_matrix_chirp_seg_Win_FFT_mean_Filtered_Ant2));
    Q_signal_Ant2 = imag(ifft(raw_data_matrix_chirp_seg_Win_FFT_mean_Filtered_Ant2));

    I_Ant2 = [I_Ant2; I_signal_Ant2.'];
    Q_Ant2 = [Q_Ant2; Q_signal_Ant2.'];
    
    %======================================================================


     %=====> Antenna 3 Data

    raw_data_matrix_chirp_seg_Ant3 = raw_data_matrix_SubFrame_Ant3(:,i_seg);
    raw_data_matrix_chirp_seg_Ant3 = reshape(raw_data_matrix_chirp_seg_Ant3,Samps_p_Chirp*zero_padding_factor,[]);
    w = window(@blackmanharris,Samps_p_Chirp*zero_padding_factor);
    raw_data_matrix_chirp_seg_Win_Ant3 = raw_data_matrix_chirp_seg_Ant3.*w;
    raw_data_matrix_chirp_seg_Win_FFT_Ant3 = fft(raw_data_matrix_chirp_seg_Win_Ant3,Samps_p_Chirp*zero_padding_factor);
    
%     figure
%     plot(abs(raw_data_matrix_chirp_seg_Win_FFT_Ant3))
%     pause

    %raw_data_matrix_chirp_seg_Win_FFT_Ant3_mean = mean(raw_data_matrix_chirp_seg_Win_FFT_Ant3.');
    if Nb_Chirp_Mean_FFT ==1
        raw_data_matrix_chirp_seg_Win_FFT_Ant3_mean = raw_data_matrix_chirp_seg_Win_FFT_Ant3(:,Chips_p_Frame/2).';
    elseif Nb_Chirp_Mean_FFT ~= 1
        raw_data_matrix_chirp_seg_Win_FFT_Ant3_mean = mean(raw_data_matrix_chirp_seg_Win_FFT_Ant3(:,1:Nb_Chirp_Mean_FFT).');
    end

    filter_pic_Ant3 = zeros(1,Samps_p_Chirp*zero_padding_factor);
    %i_start = 2;%*zero_padding_factor;
    temp_Ant3 = abs(raw_data_matrix_chirp_seg_Win_FFT_Ant3_mean(1,i_start:Samps_p_Chirp/2));
    [Val_Ant3, Max_index_Ant3] = max(temp_Ant3);
    %Max_index_Ant3
    Max_index_Ant3 = 1;
    filter_pic_Ant3(1,Max_index_Ant3+i_start-1) = 1;

    %temp_index = [temp_index, Val_Ant3];

    raw_data_matrix_chirp_seg_Win_FFT_mean_Filtered_Ant3 = raw_data_matrix_chirp_seg_Win_FFT_Ant3_mean .* filter_pic_Ant3;


    I_signal_Ant3 = real(ifft(raw_data_matrix_chirp_seg_Win_FFT_mean_Filtered_Ant3));
    Q_signal_Ant3 = imag(ifft(raw_data_matrix_chirp_seg_Win_FFT_mean_Filtered_Ant3));

    I_Ant3 = [I_Ant2; I_signal_Ant3.'];
    Q_Ant3 = [Q_Ant2; Q_signal_Ant3.'];
    
    %======================================================================



    waitbar(i_seg/length(raw_data_matrix_SubFrame_Ant1(1,:)));

end

close(h)

%% DoA Estimation
% I_Ant1 = I_Ant1 - mean(I_Ant1);
% Q_Ant1 = Q_Ant1 - mean(Q_Ant1);
% I_Ant1 = I_Ant1./max(abs(I_Ant1));
% Q_Ant1 = Q_Ant1./max(abs(Q_Ant1));
% 
% I_Ant2 = I_Ant2 - mean(I_Ant2);
% Q_Ant2 = Q_Ant2 - mean(Q_Ant2);
% I_Ant2 = I_Ant2./max(abs(I_Ant2));
% Q_Ant2 = Q_Ant2./max(abs(Q_Ant2));
% 
% I_Ant3 = I_Ant3 - mean(I_Ant3);
% Q_Ant3 = Q_Ant3 - mean(Q_Ant3);
I_Ant3 = I_Ant3(1:length(I_Ant1));
Q_Ant3 = Q_Ant3(1:length(I_Ant1));
% I_Ant3 = I_Ant3./max(abs(I_Ant3));
% Q_Ant3 = Q_Ant3./max(abs(Q_Ant3));

Ant_confg = input("Please Choose Antennas Configuration :\n 1- For Antenna 1;\n 2- For Antenna 2;\n 3- For Antenna 3;\n 4- For Antennas [1 2];\n 5- For Antennas [1 3];\n 6- For Antennas [2 3]: \n");

if Ant_confg == 1
    Vec_IQ = [complex(I_Ant1,Q_Ant1).';complex(I_Ant1,Q_Ant1).'];
elseif Ant_confg == 2
    Vec_IQ = [complex(I_Ant2,Q_Ant2).';complex(I_Ant2,Q_Ant2).'];
elseif Ant_confg == 3
    Vec_IQ = [complex(I_Ant3,Q_Ant3).';complex(I_Ant3,Q_Ant3).'];
elseif Ant_confg == 4
    Vec_IQ = [complex(I_Ant1,Q_Ant1).';complex(I_Ant2,Q_Ant2).'];
elseif Ant_confg == 5
    Vec_IQ = [complex(I_Ant1,Q_Ant1).';complex(I_Ant3,Q_Ant3).'];
elseif Ant_confg == 6
    Vec_IQ = [complex(I_Ant2,Q_Ant2).';complex(I_Ant3,Q_Ant3).'];
end



Vec_IQ_p_Frame = reshape(Vec_IQ,2,[],length(raw_data_matrix_SubFrame_Ant1(1,:)));

h1 = waitbar(0,'Please wait...DoA Estimation');
for i_SubFrame = 1:length(Vec_IQ_p_Frame(1,1,:))

    Rx = inv(length(Vec_IQ_p_Frame(1,:,1))) * Vec_IQ_p_Frame(:,:,i_SubFrame) * Vec_IQ_p_Frame(:,:,i_SubFrame)';
    [eigVecMx,eigVMx] = eig(Rx);
    eigV_Vec = diag(eigVMx);
    [Val_eigV, Max_index_eigV] = max(eigV_Vec);
    A = eigVecMx(:,Max_index_eigV);

    
    for i_alpha = 1:length(Vec_alpha)

        alpha = Vec_alpha(i_alpha);
        a_alpha = [1; exp(-1i*pi*sin(alpha))];

        a_proj = A * A' * a_alpha;
        
%         w_proj = (Rx \ a_proj) / (a_proj' * (Rx \ a_proj));
%         psd(i_alpha) = real(w_proj' * Rx * w_proj);

        psd(i_alpha) = real(a_proj' * Rx * a_proj);

        
    end

    

    [Val_Doa, Max_index_Doa] = max(psd);

    DoA(i_SubFrame) = Vec_alpha(Max_index_Doa);

    waitbar(i_SubFrame/length(Vec_IQ_p_Frame(1,1,:)));
    %plot(Vec_alpha_degree,psd);pause

end
close(h1)

DoA_mean = mean(DoA);

temp_abr = [];
theta = [];


for i_SubFrame = 1:length(Vec_IQ_p_Frame(1,1,:))

    %a_alpha_max = [1; exp(-1i*pi*sin(DoA(i_SubFrame)))];
    a_alpha_max = [1; exp(-1i*pi*sin(DoA_mean))];

    Vec_IQ_p_Frame_BF = a_alpha_max' * Vec_IQ_p_Frame(:,:,i_SubFrame);
    
    I = real(Vec_IQ_p_Frame_BF).';
    Q = imag(Vec_IQ_p_Frame_BF).';

%     I = real(Vec_IQ_p_Frame(1,:,i_SubFrame)).';
%     Q = imag(Vec_IQ_p_Frame(1,:,i_SubFrame)).';
    

%     for i_Frame = 1:length(I(1,:))
    
        XY = [I, Q];
        Par = TaubinSVD(XY);
        %Par = TaubinNTN(XY);
        
        if i_SubFrame == 1
            temp_abr = [Par(1),Par(2),Par(3)];
            abr = temp_abr;
            
        elseif i_SubFrame <= 30
            temp_abr = [[Par(1),Par(2),Par(3)];temp_abr];
            abr = mean(temp_abr);
            
            
        elseif i_SubFrame > 30
            temp_abr = [[Par(1),Par(2),Par(3)]; temp_abr(2:end,:)];
            abr = mean(temp_abr);
        end
        
    %     [Par(1),Par(2),Par(3)]
    %     pause
       
        IQ_compensated = (1/abr(3)) * complex(I+abr(1),Q+abr(2));
        %IQ_compensated = 1 * complex(I,Q);
%         figure
%         scatter(real(IQ_compensated),imag(IQ_compensated))
%         hold on
%         scatter(I(:,i_SubFrame),Q(:,i_SubFrame),'x')
%         pause
        temp_theta = unwrap(angle(IQ_compensated));
        theta = [theta; mean(temp_theta)];
        %theta = [theta; temp_theta(end)];
        %theta = [theta; abr(3)];
    
%     end
end




% R_vec = length()


diff_theta = diff(unwrap(theta));

[y2,d2] = bandstop([diff_theta],[0.9 1.1],inv(Radar_Parameter.Frame_Period_sec/1),ImpulseResponse="auto",Steepness=0.5);
% [y1,d1] = bandpass(y2,[0.05 4],inv(Radar_Parameter.Frame_Period_sec/1),ImpulseResponse="auto",Steepness=0.95);

[y1,d1] = bandpass([y2],[0.75 4],inv(Radar_Parameter.Frame_Period_sec/1),ImpulseResponse="auto",Steepness=0.95);
%[y1,d1] = bandpass([diff_theta],[0.75 4],inv(Radar_Parameter.Frame_Period_sec/1),ImpulseResponse="auto",Steepness=0.95);

%[y1,d1] = bandpass([diff_theta],[39 45],inv(Radar_Parameter.Frame_Period_sec/1),ImpulseResponse="auto",Steepness=0.95);
%[y1,d1] = bandpass([temp_index],[0.75 4],inv(Radar_Parameter.Frame_Period_sec/3),ImpulseResponse="auto",Steepness=0.95);

figure
lambda = c/fc;
tt = (1e6)*(lambda/2)/(2*pi);
% y1 = tt*y1;
% y1 = y1./max(y1(1:1e2));

%>>>>>>>>>
Group_Delay = 765; % samples
y1 = y1(Group_Delay:end-Group_Delay+1);

time = [0:length(y1)-1]*Radar_Parameter.Frame_Period_sec/1;
y1 = tt*y1;
y1 = y1./max(abs(y1));
subplot(2,1,1),plot(time,delayseq(y1,-2)),grid
title('Data : sd-ppg-n-jvp-45.raw')
xlabel('Time [sec]');ylabel('Normalized Amplitude')

disp('Please select the reference JVP pulse...')
[x_pulse,y_pulse] = ginput(2);
i_pulse = ceil(x_pulse / Radar_Parameter.Frame_Period_sec);

figure
Nb_pulses = input('Please Enter Number of Best JVP Pulses to Detect: \n');
[istart,istop,dist] = findsignal(y1,y1(i_pulse(1):i_pulse(2)),'MaxNumSegments',Nb_pulses);
findsignal(y1,y1(i_pulse(1):i_pulse(2)),'MaxNumSegments',Nb_pulses);
Peak_Points = [];
for i_dpulses = 1:length(istart)
    [Val_Max, Index_Max] = max(y1(istart(i_dpulses):istop(i_dpulses)));
    Peak_Points = [Peak_Points, istart(i_dpulses)+Index_Max];
end

figure
lambda = c/fc;
tt = (1e6)*(lambda/2)/(2*pi);
time = [0:length(y1)-1]*Radar_Parameter.Frame_Period_sec/1;
y1 = tt*y1;
y1 = y1./max(abs(y1));
subplot(2,1,1),plot(time,delayseq(y1,-2))
hold on
plot(Peak_Points*Radar_Parameter.Frame_Period_sec,y1(Peak_Points),'xr','LineWidth',2),grid
xlabel('Time [sec]');ylabel('Normalized Amplitude'),
legend('JVP Signal','Detected Peaks')
title('Data : sd-ppg-n-jvp-45.raw')

figure
NFFT = 4096;
ax_Freq = inv(Radar_Parameter.Frame_Period_sec/1)*[-NFFT/2:(NFFT/2)-1]/NFFT;
plot(ax_Freq,abs(fftshift(fft(y1,NFFT)))/NFFT),grid
xlabel('Frequency [Hz]'),ylabel('Amplitude')
title('Spectrum of Filtered JVP Signal')


figure

plot(Vec_alpha_degree,psd./max(psd),'r'),grid
xlabel('DoA [degrees]'),ylabel('Normilized PSD')

figure

plot(180*DoA/pi,'r'),grid
xlabel('Frame index'),ylabel('DoA [degrees]')

figure

plot(Radar_Parameter.Sampling_Frequency_kHz*[0:length(temp_Ant1)-1]/(2*length(temp_Ant1)),temp_Ant1/Samps_p_Chirp)
xlabel('Frequency [KHz]'),ylabel('Normalized amplitude'),grid
title('FFT-Range')

figure
temp_Frames = Frames(1:Samps_p_Frame*1);
temp_Frames(isnan(temp_Frames))=0;
pspectrum(temp_Frames,Radar_Parameter.Sampling_Frequency_kHz,'spectrogram','OverlapPercent',99,'FrequencyLimits',[0 Radar_Parameter.Sampling_Frequency_kHz/2],'TimeResolution',100*Radar_Parameter.Frame_Period_sec)
view([85 15])
title('Spectrogram FFT-Range')

figure

pspectrum(temp_Frames,Radar_Parameter.Sampling_Frequency_kHz,'persistence','OverlapPercent',99,'FrequencyLimits',[0 Radar_Parameter.Sampling_Frequency_kHz/2],'TimeResolution',100*Radar_Parameter.Frame_Period_sec)
title('Spectrogram FFT-Range')

%==========================================================================
%==========================================================================

