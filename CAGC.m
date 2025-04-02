clear,clc;
% CAGC : C/A code generate and capture
% option ---------------------------------------------
Local_Prn    = 2;               % 待生成C/A码的卫星序列号，其C/A码储存在Local_CAcode中
Contrast_Prn = 1;               % 与Local进行C/A码互相关的卫星序列号，其C/A码储存在Contrast_CAcode中
Intermediate_Freq = 4092000;    % 理论中频:4.092MHz
Step_Size    = 500;            % 搜索步长:Hz
Search_Scope = 0000;           % 搜索范围:±Hz
% C/A code generate
Local_CAcode = CAcodeGenerate(Local_Prn);
Contrast_CAcode = CAcodeGenerate(Contrast_Prn);
function CAcode = CAcodeGenerate(Prn)
    CAcode = zeros(1023,1);
    G1=ones(10,1);%状态
    G2=ones(10,1);
    G1_Code=-1*ones(1023,1);%序列
    G2_Code=-1*ones(1023,1);
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
    Prn_Offset=[
        5;
        6;
        7;
        8;
        17;
        18;
        139;
        140;
        141;
        251;
        252;
        254;
        255;
        256;
        257;
        258;
        469;
        470;
        471;
        472;
        473;
        474;
        509;
        512;
        513;
        514;
        515;
        516;
        859;
        860;
        861;
        862;
    ];
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
        if(i+Prn_Offset(Prn)<=1023)
            G2_Code(i+Prn_Offset(Prn))=XOR(G2(Prn_Selection(Prn,1)),G2(Prn_Selection(Prn,2)));
        else
            G2_Code(i+Prn_Offset(Prn)-1023)=XOR(G2(Prn_Selection(Prn,1)),G2(Prn_Selection(Prn,2)));
        end
        %寄存器移位
        for j=10:-1:2
            G1(j)=G1(j-1);
            G2(j)=G2(j-1);
        end
        G1(1)=G1_temp;
        G2(1)=G2_temp;
    end
    for i = 1:1023
        CAcode(i)=XOR(G1_Code(i),G2_Code(i));
    end
end

function Result = XOR(one,two)
    if(one*two<0)
        error=1
    end
    if(one~=two)
        Result=1;
    else
        Result=0;
    end
end