function y = FbgFunc( x )
%FbgFunc 计算抗体-抗原亲合力（适应度函数，目标函数）
%   根据目标函数形式编制本文件，为了保证亲合力不能为负数，需要对目标函数进行适当修改
%by dxb 20140808

assignin('base','ZIO',x(1));
assignin('base','ZDO',x(2));
assignin('base','ZP',x(3));
assignin('base','ZI',x(4));
assignin('base','ZD',x(5));
% assignin('base','YIO',x(6));
% assignin('base','YDO',x(7));
% assignin('base','YP',x(8));
% assignin('base','YI',x(9));
% assignin('base','YD',x(10));
% assignin('base','RIO',x(11));
% assignin('base','RDO',x(12));
% assignin('base','RP',x(13));
% assignin('base','RI',x(14));
% assignin('base','RD',x(15));
try
[t,xx,yy]=sim('SiFenZhiYiWeiYiSiFuKongZhi_opt',10);
catch
    yy = 1e+14;
end
y=yy(end);
% y=sin(x(1))^2+16*exp(sin(x(2)))-x(2)^5+x(5)*x(4)^2*cos(x(3));

end

