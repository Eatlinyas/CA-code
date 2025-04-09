clear,clc;
close all;  % 安全关闭所有绘图窗口
% CAGC : C/A code generate and capture
% option ----------------------------------------------------------
Opt.Local_Prn    = 30;              % 待生成C/A码的卫星序列号，其C/A码储存在Local_CAcode中
Opt.Contrast_Prn = 30;              % 与Local进行C/A码互相关的卫星序列号，其C/A码储存在Contrast_CAcode中
Opt.Sampling_Freq= 16.368e6;        % 采样率:16.368MHz
Opt.Intermediate_Freq = 4.092e6;    % 理论中频:4.092MHz
Opt.Code_Step    = 1;               % 伪距码搜索步长：chips
Opt.Freq_Step    = 50;              % 多普勒频率搜索步长:Hz
Opt.Search_Range = 10e3;            % 搜索范围:±kHz
Opt.IfShowSR     = 2;               % 是否显示搜索结果图（0：不显示任何/ 1：仅显示搜索成功的图/ 2：显示所有）
Opt.IfShowPB     = 1;               % 是否显示progress bar（进度条）
Opt.IfShowDPB    = 0;               % 是否显示detail of progress bar（进度条的详细信息）
Opt.MM_Threshold = 20;              % Max / Mean Threshold，峰值比平均的阈值，大于该阈值认为搜索成功
Opt.CA_Mode      = -0;              % C/A码表现形式，0:1或0 / -1：1或-1
Opt.IfPaint      = 1;               % 是否画互相关图，1：是 / 0：否
% 并行计算   需要并行工具箱（Parallel Computing Toolbox）
Opt.IfPal        = 0;               % 是否使用并行计算，可显著加快计算速度，但此时IfShowSR将为0，无法显示绘图
Opt.Parpoor      = -1;              % 无需自行修改，推荐设置为系统默认值
load gps_data_20ms.mat;             % 打开文件
% C/A code generate 使用抽头的组合选择生成C/A------------------------
Local_CAcode = CAcodeGenerate(Opt.Local_Prn,Opt.CA_Mode);
Contrast_CAcode = CAcodeGenerate(Opt.Contrast_Prn,Opt.CA_Mode);

% Cross Correlation------------------------------------------------
Cor = Cross_Correlation(Local_CAcode,Contrast_CAcode,Opt);

% write------------------------------------------------------------
% fid = fopen('output.bin', 'wb');
% fwrite(fid, gps_dat, 'int16');
% capture----------------------------------------------------------
K=1;  %选取第k个1ms
SamplesPerCode = round(1023*Opt.Sampling_Freq/1023000);
Incoming_1ms_IF=gps_dat((K-1)*SamplesPerCode+1:K*SamplesPerCode);
if(Opt.IfPal)
    % delete(gcp('nocreate'));
    % parpool(Opt.Parpoor)
    parfor i=1:32
        Results(i) = Acquire_GPS_Signal(Incoming_1ms_IF, Opt, i);
    end
else
    for i=1:1:32
        Results(i) = Acquire_GPS_Signal(Incoming_1ms_IF, Opt, i);
    end
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

function Result = Cross_Correlation(Local,Contrast,Opt)
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
    % painting---------------------------------------------------------
    if(Opt.IfPaint)
        x=-511:1:511;
        plot(x,Result,'.-');      % 画折线图
        grid on;               % 添加网格
        title('互相关图');
        xlabel('偏移量');
        ylabel('值');
    end
end

function [acquisition_results] = Acquire_GPS_Signal(Incoming_1ms_IF, Opt, Prn)
    % Parameters
    IFreq = Opt.Intermediate_Freq;  % IF frequency
    Search_Range = -Opt.Search_Range:Opt.Freq_Step:Opt.Search_Range;  % Doppler frequency range
    Code_Range = 0:Opt.Code_Step:1022;
    Samples_Per_ms = length(Incoming_1ms_IF);
    Ts = 1/Opt.Sampling_Freq;

    Samples_Per_Chip = Samples_Per_ms/1023;

    Incoming_1ms_IF = Incoming_1ms_IF(:).';
    t = (0:Samples_Per_ms-1) * Ts;

    results = zeros(length(Search_Range), length(Code_Range));

    % Generate base C/A code
    code_chips = CAcodeGenerate(Prn, Opt.CA_Mode);
    chip_index = floor(t * 1023 / 1e-3);
    chip_index = mod(chip_index, 1023);
    base_code = code_chips(chip_index + 1);

    % Precompute all shifted codes (each row is one shift)
    shifted_codes = zeros(length(Code_Range), Samples_Per_ms);
    for code_idx = 1:length(Code_Range)
        shift_samples = round(Code_Range(code_idx) * Samples_Per_Chip);
        shifted_codes(code_idx, :) = circshift(base_code, shift_samples);
    end

    % Progress indicator
    if(Opt.IfShowPB)
        fprintf('Starting acquisition for PRN %d\n', Prn);
    end
    total_iterations = length(Search_Range);

    for fd_idx = 1:length(Search_Range)
        fd = Search_Range(fd_idx);
        if ((mod(fd_idx, 10) == 0)&&(Opt.IfShowDPB))
            fprintf('Progress: %.1f%%\n', 100 * fd_idx / total_iterations);
        end

        % Carrier
        carrier_i = cos(2*pi*(IFreq + fd)*t);
        carrier_q = sin(2*pi*(IFreq + fd)*t);

        % Mix
        signal_i = Incoming_1ms_IF .* carrier_i;
        signal_q = Incoming_1ms_IF .* carrier_q;

        % Vectorized correlation over all code phases
        I_vec = shifted_codes * signal_i.';  % [code_phase x 1]
        Q_vec = shifted_codes * signal_q.';  % [code_phase x 1]
        results(fd_idx, :) = I_vec'.^2 + Q_vec'.^2;
    end

    % Find best result
    [max_val, idx] = max(results(:));
    [doppler_idx, code_idx] = ind2sub(size(results), idx);
    Best_Doppler = Search_Range(doppler_idx);
    Best_Code_Phase = Code_Range(code_idx);
    results_no_max = results(results ~= max_val);

    if(Opt.IfShowDPB)
        fprintf('\nAcquisition Results for PRN %d:\n', Prn);
        fprintf('Best Doppler: %.1f Hz\n', Best_Doppler);
        fprintf('Best Code Phase: %.1f chips\n', Best_Code_Phase);
        fprintf('Peak to Mean Ratio: %.1f\n', max_val/mean(results(:)));
    end
    % Package result
    % acquisition_results.correlation_matrix = results;
    % acquisition_results.doppler_axis = Search_Range;
    % acquisition_results.code_axis = Code_Range;
    acquisition_results.Best_Doppler = Best_Doppler;
    acquisition_results.Best_Code_Phase = Best_Code_Phase;
    acquisition_results.Peak_Value = max_val;
    acquisition_results.Peak2Mean = max_val/mean(results_no_max);

    % Optional plot
    if(Opt.IfShowSR)
        if(Opt.IfShowSR>1||acquisition_results.Peak2Mean>Opt.MM_Threshold)
            figure();
            surf(Code_Range, Search_Range/1e3, results);
            xlabel('Code Phase (chips)');
            ylabel('Doppler Frequency (kHz)');
            zlabel('Correlation Power');
            title(sprintf('GPS Signal Acquisition Results for PRN %d', Prn));
            colorbar;
            view(45, 45);
        end
    end
end