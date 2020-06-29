% Ƶ�׷���
function [Y, f] = FFTAnalysis(y, Ts)
Fs = 1/Ts;
L = length(y);
NFFT = 2^nextpow2(L);

y = y - mean(y);
Y = fft(y, NFFT)/L;
Y = 2*abs(Y(1:NFFT/2+1));
f = Fs/2*linspace(0, 1, NFFT/2+1);
end
% ��������������������������������
% ��Ȩ����������ΪCSDN������summer_upc����ԭ�����£���ѭCC 4.0 BY-SA��ȨЭ�飬ת���븽��ԭ�ĳ������Ӽ���������
% ԭ�����ӣ�https://blog.csdn.net/summer_upc/article/details/50955724