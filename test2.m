%% 
clear all;clf;clc;
%%%%%��������%%%%
filename = 'C:\Users\dell\Desktop\����\1213_Plot_and_Store_Rep_1.0.csv';
M = csvread(filename,1,0);
t=M(1:10000,1);
x1=M(1:10000,2);
%% 
Fs=1000;        %����Ƶ��
Ts=1/Fs;
%%   ��ͨ�˲���
fp=500;fs=1000;                 %ͨ����ֹƵ�ʣ������ֹƵ��  
rp=30;rs=40;                    %ͨ�������˥��  

wp=2*pi*fp;  ws=2*pi*fs;     %���Ƶ��
[n,wn]=buttord(wp,ws,rp,rs,'s');     %��s����ȷ��������˹ģ���˲����״κ�3dB ��ֹģ��Ƶ��  
[z,P,k] = buttap(n);       %��ƹ�һ��������˹ģ���ͨ�˲�����zΪ���㣬pΪ����kΪ����  
[bp,ap] = zp2tf(z,P,k);    %ת��ΪHa(p),bpΪ����ϵ����apΪ��ĸϵ��  
[bs,as] = lp2lp(bp,ap,wp); %Ha(p)ת��Ϊ��ͨHa(s)��ȥ��һ����bsΪ����ϵ����asΪ��ĸϵ��  
  
[hs,ws]=freqs(bs,as);         %ģ���˲����ķ�Ƶ��Ӧ  
[bz,az]=bilinear(bs,as,Fs);   %��ģ���˲���˫���Ա任  
[h1,w1]=freqz(bz,az);         %�����˲����ķ�Ƶ��Ӧ  
m1=filter(bz,az,x1);  
%   ********��ͨ�˲����ź�**********
figure  
subplot(2,1,1);  
plot(M(1:10000,1)*1000,M(1:10000,2));axis([450 1000 -5.5*10^-4 5.5*10^-4 ]);
xlabel('t(s)');ylabel('mv');title('ԭʼSEMG�źŲ���');grid;  
subplot(2,1,2);  
plot(M(1:10000,1)*1000,m1);axis([450 1000 -5.5*10^-4 5.5*10^-4 ]);
xlabel('t(s)');ylabel('mv');title('��ͨ�˲����ʱ��ͼ��');grid;  


%%    �����������C�����˲������ƹ�Ƶ���š�����������-  
%50Hz�ݲ�������һ����ͨ�˲�������һ����ͨ�˲������  
%����ͨ�˲�����һ��ȫͨ�˲�����ȥһ����ͨ�˲�������  
Me=100;               %�˲�������  
L=100;                %���ڳ���  
beta=100;             %˥��ϵ��  

 
wc1=49/Fs*pi;         %wc1Ϊ��ͨ�˲�����ֹƵ�ʣ���Ӧ51Hz  
wc2=51/Fs*pi;         %wc2Ϊ��ͨ�˲�����ֹƵ�ʣ���Ӧ49Hz  

h=ideal_lp(0.132*pi,Me)-ideal_lp(wc1,Me)+ideal_lp(wc2,Me); %hΪ�ݲ��������Ӧ  
w=kaiser(L,beta);  
y1=h.*rot90(w);         %yΪ50Hz�ݲ��������Ӧ����  
m2=filter(y1,1,m1);     

    %    **********�ݲ��������ף�Ƶ��************
 figure  
 subplot(2,1,1);
 plot(abs(h));axis([0 100 0 0.2]);  xlabel('Ƶ��(Hz)');ylabel('����(mv)');title('�ݲ���������');grid;  
    N=512;   
    P=10*log10(abs(fft(y1).^2)/N);  
    f=(0:length(P)-1);  
subplot(2,1,2);
plot(f,P);  xlabel('Ƶ��(Hz)');ylabel('����(dB)');title('�ݲ���������');grid;     

    %     ************�˲����ź�ͼ**********
figure  
subplot (4,1,1); 
plot(M(1:10000,1)*1000,m1); axis([0 6000 -5.5*10^-4 5.5*10^-4 ]);  xlabel('t(s)');ylabel('��ֵ');title('ԭʼ�ź�');grid;  
subplot(4,1,2);
plot(M(1:10000,1)*1000,m2); axis([0 6000 -5.5*10^-4 5.5*10^-4 ]);  xlabel('t(s)');ylabel('��ֵ');title('�����˲����ź�');grid;    
subplot (4,1,3); 
plot(M(1:10000,1)*1000,m1);axis([450 2300 -5.5*10^-4 5.5*10^-4 ]);  xlabel('t(s)');ylabel('��ֵ');title('ԭʼ�ź�');grid;  
subplot(4,1,4);
plot(M(1:10000,1)*1000,m2);axis([450 2300 -5.5*10^-4 5.5*10^-4 ]);  xlabel('t(s)');ylabel('��ֵ');title('�����˲����ź�');grid;  
%% EMD�ֽ�
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
title(sprintf('IMF%d��Ƶ��', k))
xlabel('f/Hz')
ylabel('|IMF(f)|');

subplot(323)
plot(t, yenvelope)
title(sprintf('IMF%d�İ���', k))
xlabel('Time/s')
ylabel('envelope');

subplot(324)
plot(t(1:end-1), yfreq)
title(sprintf('IMF%d��˲ʱƵ��', k))
xlabel('Time/s')
ylabel('Frequency/Hz');

subplot(325)
plot(t, yModulate)
title(sprintf('IMF%d�ĵ����ź�', k))
xlabel('Time/s')
ylabel('modulation');

subplot(326)
plot(f, YMf)
title(sprintf('IMF%d�����źŵ�Ƶ��', k))
xlabel('f/Hz')
ylabel('|YMf(f)|');

%%   EMD�˲�








%%  ������ȡ
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
xlabel('ȡ������/1');
ylabel('�����ѹֵ/uV');







