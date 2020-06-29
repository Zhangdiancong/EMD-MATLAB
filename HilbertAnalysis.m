% Hilbert分析
function [yenvelope, yf, yh, yangle] = HilbertAnalysis(y, Ts)
yh = hilbert(y);
yenvelope = abs(yh);% 包络
yangle = unwrap(angle(yh));% 相位
yf = diff(yangle)/2/pi/Ts;% 瞬时频率
end
% ————————————————
% 版权声明：本文为CSDN博主「summer_upc」的原创文章，遵循CC 4.0 BY-SA版权协议，转载请附上原文出处链接及本声明。
% 原文链接：https://blog.csdn.net/summer_upc/article/details/50955724