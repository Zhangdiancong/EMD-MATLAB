% Hilbert����
function [yenvelope, yf, yh, yangle] = HilbertAnalysis(y, Ts)
yh = hilbert(y);
yenvelope = abs(yh);% ����
yangle = unwrap(angle(yh));% ��λ
yf = diff(yangle)/2/pi/Ts;% ˲ʱƵ��
end
% ��������������������������������
% ��Ȩ����������ΪCSDN������summer_upc����ԭ�����£���ѭCC 4.0 BY-SA��ȨЭ�飬ת���븽��ԭ�ĳ������Ӽ���������
% ԭ�����ӣ�https://blog.csdn.net/summer_upc/article/details/50955724