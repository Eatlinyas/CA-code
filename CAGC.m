clear,clc;
% CAGC : C/A code generate and capture
% option ----------------------------------------------------------
Local_Prn    = 10;              % 待生成C/A码的卫星序列号，其C/A码储存在Local_CAcode中
Contrast_Prn = 10;              % 与Local进行C/A码互相关的卫星序列号，其C/A码储存在Contrast_CAcode中
Sampling_Freq= 16368000;        % 采样率:16.368MHz
Intermediate_Freq = 4092000;    % 理论中频:4.092MHz
Step_Size    = 25;             % 搜索步长:Hz
Search_Scope = 5000;            % 搜索范围:±Hz
CA_Mode      = 00;              % C/A码表现形式，0:1或0/-1：1或-1
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
% capture----------------------------------------------------------
% Result_FFT = fft(gps_dat);
% Amp_FFT = abs(Result_FFT);
% plot(Amp_FFT);

N = 327360;          
Fs = Sampling_Freq;      
df = Fs / N;         % 频率分辨率 = 50 Hz
f = (0:N-1) * df; 

for i=1:1:N
    if(gps_dat(i)>0)
        gps_dat(i)=gps_dat(i)-1;
    else
        gps_dat(i)=gps_dat(i)+1;
    end
end

% 计算 FFT
Y = abs(fft(gps_dat));
Y = Y(1:N/2); % 只取正频率部分
f = f(1:N/2); % 频率轴对应
threshold = max(Y) * 0.1;
% 找到主峰
[peaks, locs] = findpeaks(Y, 'MinPeakHeight', threshold);
freqs = f(locs);  % 主频点
% 绘制
plot(f, Y);
hold on;
plot(freqs, peaks, 'r*'); % 标出峰值
hold off;
% C/A捕获
CA_Freq = max(freqs);
Step = Sampling_Freq/CA_Freq;
t = (0:length(Y)-1)/Sampling_Freq;  % 时间轴
% 下变频，将 4.42 MHz 信号移到基带
Y = Y';
Y_baseband = Y .* exp(-1j * 2 * pi * CA_Freq * t)';  

% 低通滤波，去掉高频分量
LPF = designfilt('lowpassfir', 'FilterOrder', 100, ...
                 'CutoffFrequency', 1e6, ...  % 设定低通截止频率
                 'SampleRate', Sampling_Freq);
Y_filtered = filter(LPF, real(Y_baseband));

% Y_filtered = filtfilt(LPF, real(Y_baseband)); 
% freqz(LPF);  % 查看滤波器频率响应


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