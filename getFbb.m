function [ Fbbf ] = getFbb( pop, PopSize, rng )
%getFbb ����һ�鿹���п���-�����׺Ͷ� 
%   pop - ���鿹��Ļ�����Ϊĳ����Ļ���Σ���Ϊ��ͬ����
%   popSize - ���鿹��ĸ���
%   Fbbf - �������鿹���Ӧ����ÿ���-�����׺Ͷ�
%by dxb 20140807

    Fbbf = zeros(PopSize,PopSize-1);   %�����׺ͶȾ���ռ�
    parfor i= 1:PopSize
        Pop = pop;
        tempRow = zeros(1,PopSize);
        for j=1:1:PopSize
%             tempRow(j)=exp(sqrt(sum((pop(j,:)-Pop(i,:)).^2))/sum(rng));  %����ŷ����þ��� 
            tempRow(j)=exp(sqrt(sum(((pop(j,:)-Pop(i,:))./rng).^2)));  %����ŷ����þ��� 
        end
        tempRow(i) = [];
        Fbbf(i,:)=tempRow;
    end
end
