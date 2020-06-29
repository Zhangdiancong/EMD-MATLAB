%% 
clear all;clf;clc;
%%%%%下载数据%%%%
filename = 'C:\Users\dell\Desktop\数据\1213_Plot_and_Store_Rep_1.0.csv';
M = csvread(filename,1,0);
t=M(1:10000,1);
x1=M(1:10000,2);
%% 
Fs=1000;        %采样频率
Ts=1/Fs;
%%   低通滤波器
fp=500;fs=1000;                 %通带截止频率，阻带截止频率  
rp=30;rs=40;                    %通带、阻带衰减  

wp=2*pi*fp;  ws=2*pi*fs;     %求角频率
[n,wn]=buttord(wp,ws,rp,rs,'s');     %’s’是确定巴特沃斯模拟滤波器阶次和3dB 截止模拟频率  
[z,P,k] = buttap(n);       %设计归一化巴特沃斯模拟低通滤波器，z为极点，p为零点和k为增益  
[bp,ap] = zp2tf(z,P,k);    %转换为Ha(p),bp为分子系数，ap为分母系数  
[bs,as] = lp2lp(bp,ap,wp); %Ha(p)转换为低通Ha(s)并去归一化，bs为分子系数，as为分母系数  
  
[hs,ws]=freqs(bs,as);         %模拟滤波器的幅频响应  
[bz,az]=bilinear(bs,as,Fs);   %对模拟滤波器双线性变换  
[h1,w1]=freqz(bz,az);         %数字滤波器的幅频响应  
m1=filter(bz,az,x1);  
%   ********低通滤波后信号**********
figure  
subplot(2,1,1);  
plot(M(1:10000,1)*1000,M(1:10000,2));axis([450 1000 -5.5*10^-4 5.5*10^-4 ]);
xlabel('t(s)');ylabel('mv');title('原始SEMG信号波形');grid;  
subplot(2,1,2);  
plot(M(1:10000,1)*1000,m1);axis([450 1000 -5.5*10^-4 5.5*10^-4 ]);
xlabel('t(s)');ylabel('mv');title('低通滤波后的时域图形');grid;  


%%    —————–带陷滤波器抑制工频干扰——————-  
%50Hz陷波器：由一个低通滤波器加上一个高通滤波器组成  
%而高通滤波器由一个全通滤波器减去一个低通滤波器构成  
Me=100;               %滤波器阶数  
L=100;                %窗口长度  
beta=100;             %衰减系数  

 
wc1=49/Fs*pi;         %wc1为高通滤波器截止频率，对应51Hz  
wc2=51/Fs*pi;         %wc2为低通滤波器截止频率，对应49Hz  

h=ideal_lp(0.132*pi,Me)-ideal_lp(wc1,Me)+ideal_lp(wc2,Me); %h为陷波器冲击响应  
w=kaiser(L,beta);  
y1=h.*rot90(w);         %y为50Hz陷波器冲击响应序列  
m2=filter(y1,1,m1);     

    %    **********陷波器功率谱，频谱************
 figure  
 subplot(2,1,1);
 plot(abs(h));axis([0 100 0 0.2]);  xlabel('频率(Hz)');ylabel('幅度(mv)');title('陷波器幅度谱');grid;  
    N=512;   
    P=10*log10(abs(fft(y1).^2)/N);  
    f=(0:length(P)-1);  
subplot(2,1,2);
plot(f,P);  xlabel('频率(Hz)');ylabel('功率(dB)');title('陷波器功率谱');grid;     

    %     ************滤波后信号图**********
figure  
subplot (4,1,1); 
plot(M(1:10000,1)*1000,m1); axis([0 6000 -5.5*10^-4 5.5*10^-4 ]);  xlabel('t(s)');ylabel('幅值');title('原始信号');grid;  
subplot(4,1,2);
plot(M(1:10000,1)*1000,m2); axis([0 6000 -5.5*10^-4 5.5*10^-4 ]);  xlabel('t(s)');ylabel('幅值');title('带阻滤波后信号');grid;    
subplot (4,1,3); 
plot(M(1:10000,1)*1000,m1);axis([450 2300 -5.5*10^-4 5.5*10^-4 ]);  xlabel('t(s)');ylabel('幅值');title('原始信号');grid;  
subplot(4,1,4);
plot(M(1:10000,1)*1000,m2);axis([450 2300 -5.5*10^-4 5.5*10^-4 ]);  xlabel('t(s)');ylabel('幅值');title('带阻滤波后信号');grid;  
%% EMD分解
x= m2;
imf = emd(x);
plot_hht(x,imf,1/Fs);

k = 4;
y = imf{k};
N = length(y);

[yenvelope, yfreq, yh, yangle] = HilbertAnalysis(y, 1/Fs);
yModulate = y./yenvelope;
[YMf, f] = FFTAnalysis(yModulate, Ts);
Yf = FFTAnalysis(y, Ts);

figure
subplot(321)
plot(t, y)
title(sprintf('IMF%d', k))
xlabel('Time/s')
ylabel(sprintf('IMF%d', k));

subplot(322)
plot(f, Yf)
title(sprintf('IMF%d的频谱', k))
xlabel('f/Hz')
ylabel('|IMF(f)|');

subplot(323)
plot(t, yenvelope)
title(sprintf('IMF%d的包络', k))
xlabel('Time/s')
ylabel('envelope');

subplot(324)
plot(t(1:end-1), yfreq)
title(sprintf('IMF%d的瞬时频率', k))
xlabel('Time/s')
ylabel('Frequency/Hz');

subplot(325)
plot(t, yModulate)
title(sprintf('IMF%d的调制信号', k))
xlabel('Time/s')
ylabel('modulation');

subplot(326)
plot(f, YMf)
title(sprintf('IMF%d调制信号的频谱', k))
xlabel('f/Hz')
ylabel('|YMf(f)|');

%%   EMD滤波








%%  特征提取
thisdata=m2;
length_t=1000;
delta_t=200;
j=1;
for i=1:delta_t:size(thisdata,1)
    if i+length_t>size(thisdata,1)
        break;
    end
    meandata(j)=(sum(abs(thisdata(i:i+length_t))))/length_t;
    j=j+1;
end
figure
plot(meandata);
xlabel('取样点数/1');
ylabel('肌电电压值/uV');







