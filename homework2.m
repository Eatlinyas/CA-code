clear,clc;
close all;  % 安全关闭所有绘图窗口
%-----------------------------------------
Gg = 10;            % 载体加加速度
Time = 20e-3;       % 相干积分时间 s
BandWidth = 10;      % 带宽 Hz
Raw_CN = 10:1:50;   % 载噪比 dB-Hz
pi = 3.1415926535;
i = 1;
sigma = zeros(1,41);
sum = zeros(1,41);
%-----------------------------------------
CN = 0.1*exp(Raw_CN);
theta = 0.4828 * Gg /(BandWidth*BandWidth*BandWidth);
while(i<42)
    sigma(i) = 1+1/(2*Time*CN(i));
    sigma(i) = sigma(i) * BandWidth/CN(i);
    sigma(i) = sqrt(sigma(i))*180/pi;
    sum(i) = sigma(i) + theta;
    i=i+1;
end
%------------------------------------------
figure();
plot(Raw_CN,sigma,'.-');
grid on;               % 添加网格
title('热噪声误差');
xlabel('C/N dB-Hz');
ylabel('PLL 热噪声误差');

% figure();
% plot(CN,sum,'.-');
% grid on;               % 添加网格
% title('总误差');
% xlabel('C/N db-Hz');
% ylabel('总误差');