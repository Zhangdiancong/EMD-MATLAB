function n = findpeaks(x)
% Find peaks. �Ҽ���ֵ�㣬���ض�Ӧ����ֵ�������
n = find(diff(diff(x) > 0) < 0); % �൱���Ҷ��׵�С��0�ĵ�
u = find(x(n+1) > x(n));
n(u) = n(u)+1;        % ��1��������Ӧ����ֵ��

% ͼ�ν�����������
% figure
% subplot(611)
% x = x(1:100);
% plot(x, '-o')
% grid on
%
% subplot(612)
% plot(1.5:length(x), diff(x) > 0, '-o')
% grid on
% axis([1,length(x),-0.5,1.5])
%
% subplot(613)
% plot(2:length(x)-1, diff(diff(x) > 0), '-o')
% grid on
% axis([1,length(x),-1.5,1.5])
%
% subplot(614)
% plot(2:length(x)-1, diff(diff(x) > 0)<0, '-o')
% grid on
% axis([1,length(x),-1.5,1.5])
%
% n????= find(diff(diff(x) > 0) < 0);
% subplot(615)
% plot(n, ones(size(n)), 'o')
% grid on
% axis([1,length(x),0,2])
%
% u????= find(x(n+1) > x(n));
% n(u) = n(u)+1;
% subplot(616)
% plot(n, ones(size(n)), 'o')
% grid on
% axis([1,length(x),0,2])
% ��������������������������������
% ��Ȩ����������ΪCSDN������summer_upc����ԭ�����£���ѭCC 4.0 BY-SA��ȨЭ�飬ת���븽��ԭ�ĳ������Ӽ���������
% ԭ�����ӣ�https://blog.csdn.net/summer_upc/article/details/50955724