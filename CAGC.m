clear,clc;
% CAGC : C/A code generate and capture
% option ----------------------------------------------------------
Local_Prn    = 10;              % 待生成C/A码的卫星序列号，其C/A码储存在Local_CAcode中
Contrast_Prn = 10;              % 与Local进行C/A码互相关的卫星序列号，其C/A码储存在Contrast_CAcode中
Sampling_Freq= 16368000;        % 采样率:16.368MHz
Intermediate_Freq = 4092000;    % 理论中频:4.092MHz
Step_Size    = 500;              % 搜索步长:Hz
Search_Scope = 5000;            % 搜索范围:±Hz
CA_Mode      = -1;              % C/A码表现形式，0:1或0 / -1：1或-1
IfPaint      = 0;               % 是否画互相关图，1：是/0：否
load gps_data_20ms.mat;         % 打开文件
% C/A code generate 使用抽头的组合选择生成C/A------------------------
Local_CAcode = CAcodeGenerate(Local_Prn,CA_Mode);
Contrast_CAcode = CAcodeGenerate(Contrast_Prn,CA_Mode);

% Cross Correlation------------------------------------------------
Cor = Cross_Correlation(Local_CAcode,Contrast_CAcode);
% painting---------------------------------------------------------
if(IfPaint)
    x=-511:1:511;
    plot(x,Cor,'.-');      % 画折线图
    grid on;               % 添加网格
    title('互相关图');
    xlabel('偏移量');
    ylabel('值');
end
% write------------------------------------------------------------
% fid = fopen('output.bin', 'wb');
% fwrite(fid, gps_dat, 'int16');
% capture----------------------------------------------------------
SamplesPerCode = round(1023*Sampling_Freq/1023000);
incoming_1ms_IF=gps_dat(1:SamplesPerCode);
for i=1:32
    results(i) = acquire_gps_signal(incoming_1ms_IF, Sampling_Freq, i);
end



SamplesPerCode = round(1023*Sampling_Freq/1023000);
signal1=gps_dat(1:SamplesPerCode);
signal2=gps_dat(SamplesPerCode+1:2*SamplesPerCode);
signal0DC = gps_dat - mean(gps_dat);
% 计算采样周期ts
ts = 1 / Sampling_Freq;
% 载波相位点
phasePoints = (0 : (SamplesPerCode-1)) * 2 * pi * ts;
% 对于给定带宽，设定捕获带宽的个数 500Hz步进 总共2kHz
numberOfFrqBins = round(20 * 2) + 1;
% Generate all C/A codes and sample them according to the sampling freq.
CACodesTable = zeros(32,1023);
for i=1:1:32
    CACodesTable(i,:) = CAcodeGenerate(i,CA_Mode);%makeCaTable.m产生32颗指定卫星的CA码 32*samplesPerCode矩阵
end
% ====================生成用于存放结果的向量======================
% ------------------- 二维区域上所有搜索单元的计算结果 ---------------------------
% 生成一个零矩阵用于存放单颗卫星的相关计算结果，行是频率数，列是码相关值
results = zeros(numberOfFrqBins, SamplesPerCode);   % 电子信息学院20×10方阵中每个人的身高数据
% 产生一个存放载波偏移的向量，它反应出总共有几个不同的载波频率，用于后续的遍历搜索
frqBins = zeros(1, numberOfFrqBins);
% ------------------------- 经过比较之后，最终的捕获结果 -----------------------------
% 创建一个向量存放32颗卫星的载波频移结果
carrFreq     = zeros(1, 32);  % 32个学院里个子最高的人分别位于哪一排
% 创建一个向量存放32颗卫星的码相位
codePhase    = zeros(1, 32);  % 32个学院里个子最高的人分别位于哪一列
% 相关峰值比
peakMetric   = zeros(1, 32);
fprintf('(');
% ======================== 开始捕获过程 ==========================
%-------采用并行码相位捕获：所有的码一次全部搜索完，只需要遍历频率
 
% 进入外层循环，每次从卫星列表中取不同的PRN
for PRN = 1:1:32
 
    CACodeFreqDom = conj(fft(CACodesTable(PRN, :)));%某颗PRN卫星的本地CA码采样值FFT，再取共轭 
 
    % 进入内层循环，遍历CA码相关的带宽为500Hz的频率搜索区域，一共有numberOfFrqBins个
    for frqBinIndex = 1:numberOfFrqBins
        
        %找到每个频率区域的中频载波频率 
        frqBins(frqBinIndex) = Intermediate_Freq - ...
                               (20/2) * 1000 + ...
                               0.5e3 * (frqBinIndex - 1);
  %------产生正余弦信号作为本地载波，频率是上一步得到的本次循环的当前值，相位则是之前的载波相位点-------
 
        sinCarr = sin(frqBins(frqBinIndex) * phasePoints);
        cosCarr = cos(frqBins(frqBinIndex) * phasePoints);
  
  %-----------------对连续的两个1ms中频数据signal1、2进行载波剥离----------------------
        I1      = sinCarr .* signal1;   %同相支路
        Q1      = cosCarr .* signal1;   %正交支路
        I2      = sinCarr .* signal2;
        Q2      = cosCarr .* signal2;
 
  %---------------载波剥离后的1ms数据做FFT-----------------
        IQfreqDom1 = fft(I1 + 1i*Q1);
        IQfreqDom2 = fft(I2 + 1i*Q2);
 
  % 频域上，对 （载波FFT的结果） 和 （CA码FFT＋共轭后的结果） 进行相乘，时域卷积判断最大值
        convCodeIQ1 = IQfreqDom1 .* CACodeFreqDom; %改为转置
        convCodeIQ2 = IQfreqDom2 .* CACodeFreqDom; %改为转置
        
        %先IFFT，再取模，存储卷积相关结果，这个结果是一系列数值，也就是在当前频率上每个码的相关值
        acqRes1 = abs(ifft(convCodeIQ1)) .^ 2;
        acqRes2 = abs(ifft(convCodeIQ2)) .^ 2;
 
  % 寻找卷积相关最大值，复现码的相位值 完成每个频率域的最大值搜索 
        %计算出的相关值分别存放在之前生成的frqBinIndex*samplesPerCode的result矩阵中
 
        if (max(acqRes1) > max(acqRes2))        % 找两个CA信号中更靠谱的一个，将它的所有相关值存放在当前频率对应的这一行
            results(frqBinIndex, :) = acqRes1;  
        else
            results(frqBinIndex, :) = acqRes2;
        end
        
    end % 循环结束，两层循环两个end；别忘记frqBinIndex是遍历存储用的一个“指针”
end
% function---------------------------------------------------------
function CAcode = CAcodeGenerate(Prn,CA_Mode)
    CAcode = zeros(1023,1);
    G1=ones(10,1);%状态
    G2=ones(10,1);
    G1_Code=-2*ones(1023,1);%序列
    G2_Code=-2*ones(1023,1);
    Prn_Selection = [
        2, 6;   % PRN1
        3, 7;   % PRN2
        4, 8;   % PRN3
        5, 9;   % PRN4
        1, 9;   % PRN5
        2, 10;  % PRN6
        1, 8;   % PRN7
        2, 9;   % PRN8
        3, 10;  % PRN9
        2, 3;   % PRN10
        3, 4;   % PRN11
        5, 6;   % PRN12
        6, 7;   % PRN13
        7, 8;   % PRN14
        8, 9;   % PRN15
        9, 10;  % PRN16
        1, 4;   % PRN17
        2, 5;   % PRN18
        3, 6;   % PRN19
        4, 7;   % PRN20
        5, 8;   % PRN21
        6, 9;   % PRN22
        1, 3;   % PRN23
        4, 6;   % PRN24
        5, 7;   % PRN25
        6, 8;   % PRN26
        7, 9;   % PRN27
        8, 10;  % PRN28
        1, 6;   % PRN29
        2, 7;   % PRN30
        3, 8;   % PRN31
        4, 9;   % PRN32
    ];
    if ((Prn>32) || (Prn<1))
	    error('Prn must be a number between 1 and 32\n')
    end
    for i=1:1023
        %异或运算
        G1_temp=XOR(G1(3),G1(10));
        G2_temp=XOR(G2(2),G2(3));
        G2_temp=XOR(G2_temp,G2(6));
        G2_temp=XOR(G2_temp,G2(8));
        G2_temp=XOR(G2_temp,G2(9));
        G2_temp=XOR(G2_temp,G2(10));
        %输出
        G1_Code(i)=G1(10);
        G2_Code(i)=XOR(G2(Prn_Selection(Prn,1)),G2(Prn_Selection(Prn,2)));
        %寄存器移位
        for j=10:-1:2
            G1(j)=G1(j-1);
            G2(j)=G2(j-1);
        end
        G1(1)=G1_temp;
        G2(1)=G2_temp;
    end
    %G1，G2组合生成C/A码
    for i = 1:1023
        CAcode(i)=XOR(G1_Code(i),G2_Code(i));
    end
    if(CA_Mode==-1)
        CAcode = Zero2MinusOne(CAcode);
    end
end

function Result = XOR(one,two)
    if(one*two<0)
        error('XOR error\n')
    end
    if(one~=two)
        Result=1;
    else
        Result=0;
    end
end

function New_CA = Zero2MinusOne(CA)
    New_CA=zeros(1023,1);
    for i=1:1:1023
        if(CA(i)==0)
            New_CA(i)=-1;
        else
            New_CA(i)=1;
        end
    end
end

function Result = Cross_Correlation(Local,Contrast)
    Local=Zero2MinusOne(Local);
    Contrast=Zero2MinusOne(Contrast);
    Result = zeros(1023,1);
    for offset = -511:1:511
        for i=1:1:1023
            if(i+offset>1023)
                Cur=i+offset-1023;
            elseif(i+offset<=0)
                Cur=i+offset+1023;
            else
                Cur=i+offset;
            end
            Result(offset+511+1) = Result(offset+511+1) + Local(Cur)*Contrast(i);
        end
        Result(offset+511+1) = Result(offset+511+1) / 1023;
    end
end

function [acquisition_results] = acquire_gps_signal(incoming_1msIF, fs, prn)
    % Parameters
    fIF = 4.092e6;  % IF frequency in Hz
    fd_range = -10e3:250:10e3;  % Doppler frequency range
    code_step = 0.5;  % Code step in chips
    code_range = 0:code_step:1023;  % Code phase range
    samples_per_ms = length(incoming_1msIF);
    Ts = 1/fs;  % Sampling period
    
    % Calculate actual samples per chip
    samples_per_chip = samples_per_ms/1023;  % should be ~5.585 samples/chip
    
    % Ensure incoming_1msIF is a row vector
    incoming_1msIF = incoming_1msIF(:).'; 
    
    % Time vector for 1ms
    t = (0:samples_per_ms-1) * Ts;
    
    % Initialize results matrix
    results = zeros(length(fd_range), length(code_range));
    
    % Generate base C/A code
    code_chips = CAcodeGenerate(prn,-1);  % Get 1023 chips
    
    % Create sampled code sequence using direct calculation
    chip_index = floor(t * 1023 / 1e-3);  % Convert time to chip index
    chip_index = mod(chip_index, 1023);    % Wrap around for indices >= 1023
    base_code = code_chips(chip_index + 1); % +1 because MATLAB is 1-based indexing
    
    % Verify dimensions
    if length(base_code) ~= samples_per_ms
        error('Base code length (%d) does not match signal length (%d)', ...
            length(base_code), samples_per_ms);
    end
    
    % Progress indicator
    fprintf('Starting acquisition for PRN %d\n', prn);
    total_iterations = length(fd_range);
    
    % Loop through all Doppler frequencies and code phases
    for fd_idx = 1:length(fd_range)
        fd = fd_range(fd_idx);
        
        % Update progress
        if mod(fd_idx, 10) == 0
            fprintf('Progress: %.1f%%\n', 100*fd_idx/total_iterations);
        end
        
        % Generate carrier signal
        carrier_i = cos(2*pi*(fIF + fd)*t);
        carrier_q = sin(2*pi*(fIF + fd)*t);
        
        % Mix incoming signal with carrier
        signal_i = incoming_1msIF .* carrier_i;
        signal_q = incoming_1msIF .* carrier_q;
        
        for code_idx = 1:length(code_range)
            code_phase = code_range(code_idx);
            
            % Calculate circular shift for code phase
            shift_samples = round(code_phase * samples_per_chip);
            code = circshift(base_code, shift_samples);
            
            % Correlate I and Q channels
            I = sum(signal_i .* code);
            Q = sum(signal_q .* code);
            
            % Calculate correlation power
            results(fd_idx, code_idx) = sum(I.^2 + Q.^2);
        end
    end
    
    % Find maximum correlation and its location
    [max_val, idx] = max(results(:));
    [doppler_idx, code_idx] = ind2sub(size(results), idx);
    best_doppler = fd_range(doppler_idx);
    best_code_phase = code_range(code_idx);
    
    fprintf('\nAcquisition Results for PRN %d:\n', prn);
    fprintf('Best Doppler: %.1f Hz\n', best_doppler);
    fprintf('Best Code Phase: %.1f chips\n', best_code_phase);
    fprintf('Peak to Mean Ratio: %.1f\n', max_val/mean(results(:)));
    
    % Package results
    acquisition_results.correlation_matrix = results;
    acquisition_results.doppler_axis = fd_range;
    acquisition_results.code_axis = code_range;
    acquisition_results.best_doppler = best_doppler;
    acquisition_results.best_code_phase = best_code_phase;
    acquisition_results.peak_value = max_val;
    
    % Plot 3D acquisition function
    figure;
    surf(code_range, fd_range/1e3, results);
    xlabel('Code Phase (chips)');
    ylabel('Doppler Frequency (kHz)');
    zlabel('Correlation Power');
    title(sprintf('GPS Signal Acquisition Results for PRN %d', prn));
    colorbar;
    view(45, 45);
    
   
end