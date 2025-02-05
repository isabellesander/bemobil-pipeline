function [qrs_amp_raw,qrs_i_raw,delay,plots]=pan_tompkins(ecg,fs,gr,f1,f2)

if ~isvector(ecg)
    error('ecg must be a row or column vector');
end
if nargin < 3
    gr = 1;   % on default the function always plots
end
ecg = ecg(:); % vectorize

%% ======================= Initialize =============================== %
delay = 0;
skip = 0;                                                                  % becomes one when a T wave is detected
m_selected_RR = 0;
mean_RR = 0;
ser_back = 0;
ax = zeros(1,6);
if ~exist('f1','var')
    f1=5;                                                                      % cuttoff low frequency to get rid of baseline wander
end
if ~exist('f2','var')
    f2=15;                                                                     % cuttoff frequency to discard high frequency noise
end

%% ============ Noise cancelation(Filtering)( 5-15 Hz) =============== %%
if fs == 200
    % ------------------ remove the mean of Signal -----------------------%
    ecg = ecg - mean(ecg);
    %% ==== Low Pass Filter  H(z) = ((1 - z^(-6))^2)/(1 - z^(-1))^2 ==== %%
    Wn = 12*2/fs;
    N = 3;                                                                  % order of 3 less processing
    [a,b] = butter(N,Wn,'low');                                             % bandpass filtering
    ecg_l = filtfilt(a,b,ecg);
    ecg_l = ecg_l/ max(abs(ecg_l));
    %% ======================= start figure ============================= %%
    
    if gr
        plots(1) = figure('color','w');
        ax(1) = subplot(321);plot(ecg);axis tight;title('Raw signal');
        ax(2)=subplot(322);plot(ecg_l);axis tight;title('Low pass filtered');
    end
    %% ==== High Pass filter H(z) = (-1+32z^(-16)+z^(-32))/(1+z^(-1)) ==== %%
    Wn = 5*2/fs;
    N = 3;                                                                  % order of 3 less processing
    [a,b] = butter(N,Wn,'high');                                            % bandpass filtering
    ecg_h = filtfilt(a,b,ecg_l);
    ecg_h = ecg_h/ max(abs(ecg_h));
    if gr
        ax(3)=subplot(323);plot(ecg_h);axis tight;title('High Pass Filtered');
    end
else
    %%  bandpass filter for Noise cancelation of other sampling frequencies(Filtering)
    Wn=[f1 f2]*2/fs;                                                           % cutt off based on fs
    N = 3;                                                                     % order of 3 less processing
    [a,b] = butter(N,Wn);                                                      % bandpass filtering
    ecg_h = filtfilt(a,b,ecg);
    ecg_h = ecg_h/ max( abs(ecg_h));
    if gr
        plots(1) = figure('color','w');
        ax(1) = subplot(3,2,[1 2]);plot(ecg);axis tight;title('Raw Signal');
        ax(3)=subplot(323);plot(ecg_h);axis tight;title('Band Pass Filtered');
    end
end
%% ==================== derivative filter ========================== %%
% ------ H(z) = (1/8T)(-z^(-2) - 2z^(-1) + 2z + z^(2)) --------- %
if fs ~= 200
    int_c = (5-1)/(fs*1/40);
    b = interp1(1:5,[1 2 0 -2 -1].*(1/8)*fs,1:int_c:5);
else
    b = [1 2 0 -2 -1].*(1/8)*fs;
end

ecg_d = filtfilt(b,1,ecg_h);
ecg_d = ecg_d/max(ecg_d);

if gr
    ax(4)=subplot(324);plot(ecg_d);
    axis tight;
    title('Filtered with the derivative filter');
end
%% ========== Squaring nonlinearly enhance the dominant peaks ========== %%
ecg_s = ecg_d.^2;
if gr
    ax(5)=subplot(325);
    plot(ecg_s);
    axis tight;
    title('Squared');
end

%% ============  Moving average ================== %%
ecg_m = conv(ecg_s ,ones(1 ,round(0.150*fs))/round(0.150*fs));
delay = delay + round(0.150*fs)/2;

if gr
    ax(6)=subplot(326);plot(ecg_m);
    axis tight;
    title('Averaged with 30 samples length,Black noise,Green Adaptive Threshold,RED Sig Level,Red circles QRS adaptive threshold');
    axis tight;
end

%% ===================== Fiducial Marks ============================== %%
[pks,locs] = findpeaks(ecg_m,'MINPEAKDISTANCE',round(0.2*fs));

%% =================== Initialize Some Other Parameters =============== %%
LLp = length(pks);
qrs_c = zeros(1,LLp);           % amplitude of R
qrs_i = zeros(1,LLp);           % index
qrs_i_raw = zeros(1,LLp);       % amplitude of R
qrs_amp_raw= zeros(1,LLp);      % Index
nois_c = zeros(1,LLp);
nois_i = zeros(1,LLp);
SIGL_buf = zeros(1,LLp);
NOISL_buf = zeros(1,LLp);
SIGL_buf1 = zeros(1,LLp);
NOISL_buf1 = zeros(1,LLp);
THRS_buf1 = zeros(1,LLp);
THRS_buf = zeros(1,LLp);

%% Initialize the training phase (2 seconds of the signal) to determine the THR_SIG and THR_NOISE
THR_SIG = max(ecg_m(1:2*fs))*1/3;
THR_NOISE = mean(ecg_m(1:2*fs))*1/2;
SIG_LEV= THR_SIG;
NOISE_LEV = THR_NOISE;
THR_SIG1 = max(ecg_h(1:2*fs))*1/3;
THR_NOISE1 = mean(ecg_h(1:2*fs))*1/2;
SIG_LEV1 = THR_SIG1;
NOISE_LEV1 = THR_NOISE1;

%% ============ Thresholding and desicion rule ============= %%
Beat_C = 0;
Beat_C1 = 0;
Noise_Count = 0;
for i = 1 : LLp
    %% ===== locate the corresponding peak in the filtered signal === %%
    if locs(i)-round(0.150*fs)>= 1 && locs(i)<= length(ecg_h)
        [y_i,x_i] = max(ecg_h(locs(i)-round(0.150*fs):locs(i)));
    else
        if i == 1
            [y_i,x_i] = max(ecg_h(1:locs(i)));
            ser_back = 1;
        elseif locs(i)>= length(ecg_h)
            [y_i,x_i] = max(ecg_h(locs(i)-round(0.150*fs):end));
        end
    end
    %% ================= update the heart_rate ==================== %%
    if Beat_C >= 9
        diffRR = diff(qrs_i(Beat_C-8:Beat_C));
        mean_RR = mean(diffRR);
        comp =qrs_i(Beat_C)-qrs_i(Beat_C-1);
        
        if comp <= 0.92*mean_RR || comp >= 1.16*mean_RR
            THR_SIG = 0.5*(THR_SIG);
            THR_SIG1 = 0.5*(THR_SIG1);
        else
            m_selected_RR = mean_RR;
        end
        
    end
    
    if m_selected_RR
        test_m = m_selected_RR;
    elseif mean_RR && m_selected_RR == 0
        test_m = mean_RR;
    else
        test_m = 0;
    end
    
    if test_m && Beat_C > 0
        if (locs(i) - qrs_i(Beat_C)) >= round(1.66*test_m)
            [pks_temp,locs_temp] = max(ecg_m(qrs_i(Beat_C)+ round(0.200*fs):locs(i)-round(0.200*fs)));
            locs_temp = qrs_i(Beat_C)+ round(0.200*fs) + locs_temp -1;
            
            if pks_temp > THR_NOISE
                Beat_C = Beat_C + 1;
                qrs_c(Beat_C) = pks_temp;
                qrs_i(Beat_C) = locs_temp;
                
                [y_i_t,x_i_t] = max(ecg_h(locs_temp-round(0.150*fs):locs_temp));
                if y_i_t > THR_NOISE1
                    Beat_C1 = Beat_C1 + 1;
                    qrs_i_raw(Beat_C1) = locs_temp-round(0.150*fs) + x_i_t - 1; % save index of raw signal
                    qrs_amp_raw(Beat_C1) = y_i_t; % save amplitude of raw signal
                    SIG_LEV1 = 0.25*y_i_t + 0.75*SIG_LEV1; % when found with the second thres
                end
                SIG_LEV = 0.25*pks_temp + 0.75*SIG_LEV ; % when found with the second threshold
            end
        end
    end
    
    %% == find noise and QRS peaks == %%
    if pks(i) >= THR_SIG
        if Beat_C >= 3
            diffRR = diff(qrs_i(Beat_C-2:Beat_C));
            mean_RR = mean(diffRR);
        end
        
        if Beat_C > 0 && (locs(i) - qrs_i(Beat_C)) <= round(0.36*fs) % case of T wave
            Slope1 = mean(diff(ecg_m(locs(i)-round(0.075*fs):locs(i))));
            Slope2 = mean(diff(ecg_m(qrs_i(Beat_C)-round(0.075*fs):qrs_i(Beat_C))));
            if abs(Slope1) <= abs(0.5*(Slope2)) % slope less then 0.5 of previous R
                nois_c(Noise_Count) = pks(i);
                nois_i(Noise_Count) = locs(i);
                Noise_Count = Noise_Count + 1;
                skip = 1;
                NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV;
            else
                skip = 0;
            end
        end
        
        if skip == 0
            Beat_C = Beat_C + 1;
            Beat_C1 = Beat_C1 + 1; % Ensure Beat_C1 is incremented here
            qrs_c(Beat_C) = pks(i);
            qrs_i(Beat_C) = locs(i);
            
            %% bandpass filter
            qrs_i_raw(Beat_C1) = locs(i)-round(0.150*fs) + x_i - 1; % save index of raw signal
            qrs_amp_raw(Beat_C1) = y_i; % save amplitude of raw signal
            SIG_LEV1 = 0.125*y_i + 0.875*SIG_LEV1;
            
            SIG_LEV = 0.125*pks(i) + 0.875*SIG_LEV ;
        end
        
    elseif (THR_NOISE <= pks(i)) && (pks(i) < THR_SIG)
        NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV;
        NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
    elseif pks(i) < THR_NOISE
        nois_c(Noise_Count) = pks(i);
        nois_i(Noise_Count) = locs(i);
        Noise_Count = Noise_Count + 1;
        NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV;
        NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
    end
    
    % adjust threshold with SNR
    if NOISE_LEV ~= 0 || SIG_LEV ~= 0
        THR_SIG = NOISE_LEV + 0.25*(abs(SIG_LEV - NOISE_LEV));
        THR_NOISE = 0.5*(THR_SIG);
        THR_SIG1 = NOISE_LEV1 + 0.25*(abs(SIG_LEV1 - NOISE_LEV1));
        THR_NOISE1 = 0.5*(THR_SIG1);
    end
    
    SIGL_buf(i) = SIG_LEV;
    NOISL_buf(i) = NOISE_LEV;
    THRS_buf(i) = THR_SIG;
    
    SIGL_buf1(i) = SIG_LEV1;
    NOISL_buf1(i) = NOISE_LEV1;
    THRS_buf1(i) = THR_SIG1;
end

if gr
    subplot(326);
    hold on;
    plot(SIGL_buf,'r--','LineWidth',2);
    plot(NOISL_buf,'k--','LineWidth',2);
    plot(THRS_buf,'g--','LineWidth',2);
end

qrs_i_raw = qrs_i_raw(qrs_i_raw~=0);
qrs_amp_raw = qrs_amp_raw(qrs_amp_raw~=0);

% Interpolation for noisy segments
if Noise_Count > 0
    noisy_segments = nois_i(1:Noise_Count);
    clean_indices = setdiff(1:Beat_C1, noisy_segments);
    clean_R_locs = qrs_i_raw(clean_indices);
    
    % Interpolating the noisy R peaks
    noisy_R_locs = qrs_i_raw(noisy_segments);
    noisy_R_locs_interp = interp1(clean_indices, clean_R_locs, noisy_segments, 'linear', 'extrap');
    
    % Replacing noisy R peaks with interpolated values
    qrs_i_raw(noisy_segments) = noisy_R_locs_interp;
end

if gr
    linkaxes(ax,'x');
end
end
