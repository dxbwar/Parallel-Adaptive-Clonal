function [ Fbbff ] = getFbb2( pop, PopSize, rng, Var)
%getFbb 计算一组抗体中抗体-抗体亲和度 
%   pop - 这组抗体的基因，行为某抗体的基因段，列为不同抗体
%   popSize - 这组抗体的个数
%   Fbbf - 返回这组抗体对应的最好抗体-抗体亲和度
%by dxb 20140807

    Fbbf = zeros(PopSize,PopSize);   %分配亲和度矩阵空间
    Fbbff = zeros(PopSize,PopSize-1);   %分配亲和度矩阵空间
    parfor i = 1:PopSize
        Pop = pop;
        tempRow = zeros(1,PopSize);
        for j=i+1:1:PopSize
%             tempRow(j)=exp(sqrt(sum((pop(j,:)-Pop(i,:)).^2))/sum(rng));  %计算欧几里得距离 
            tempRow(j)=sqrt(sum(((pop(j,:)-Pop(i,:))./rng).^2)/Var);  %计算欧几里得距离 
        end
%         tempRow(i) = [];
        Fbbf(i,:)=tempRow;
    end
    Fbbf = Fbbf + Fbbf';
    for i = 1:PopSize               %消掉0值项
        tempRow=Fbbf(i,:);
        tempRow(i) = [];
        Fbbff(i,:) = tempRow;
    end
end
