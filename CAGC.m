clear clc;
% CAGC : C/A code generate and capture
% option ---------------------------------------------
Local_Prn    = 1;               % 待生成C/A码的卫星序列号，其C/A码储存在Local_CAcode中
Contrast_Prn = 1;               % 与Local进行C/A码互相关的卫星序列号，其C/A码储存在Contrast_CAcode中
Intermediate_Freq = 4092000;    % 理论中频:4.092MHz
Step_Size    = 1000;            % 搜索步长:1000Hz
Search_Scope = 10000;           % 搜索范围:±10000Hz
% C/A code generate
Local_CAcode = CAcodeGenerate(Local_Prn);
Contrast_CAcode = CAcodeGenerate(Contrast_Prn);
function CAcode = CAcodeGenerate(Prn)
    CAcode = zeros(1023);
    
end