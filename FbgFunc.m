function y = FbgFunc( x )
%FbgFunc ���㿹��-��ԭ�׺�������Ӧ�Ⱥ�����Ŀ�꺯����
%   ����Ŀ�꺯����ʽ���Ʊ��ļ���Ϊ�˱�֤�׺�������Ϊ��������Ҫ��Ŀ�꺯�������ʵ��޸�
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

