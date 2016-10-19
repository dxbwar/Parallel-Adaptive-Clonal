function [b, TG] = geneOps3( pop, q, numvar, rng, bound, pm, TG)
%geneOps 基因操作
%   将抗体组克隆后按概率进行变异
%   pop - 这组抗体的基因，行为某抗体的基因段，列为不同抗体
%   popSize - 这组抗体的个数
%   q - 原抗体组对应的克隆个数
%   numvar - 基因段数，变量个数
%   rng - 变量取值范围
%   bound - 变量取值边界
%   pm - 每段基因对应的变异概率
%   b - 返回变异后的抗体-抗原亲合力最优抗体组
%by dxb 20140807
    PopSize = size(pop, 1); 
    b = [];                                                    %克隆后的最优抗体集合
    mutateQ = zeros(PopSize,1);                                 %每个抗体克隆后发生变异的个数
%     b = [b,zeros(size(b,1)];                                    %加入变异标识
    Fbg = pop(:,numvar+1);                                      %提取亲合度
    avrFbg = mean(Fbg);
    pK = pm + ((1-pm)^-1-1+exp(avrFbg./Fbg)).^-1;                         %抗体变异概率计算
    
%     for i=1:numvar
%         if rand < 0.5
%             TG(i) = TG(i) / 1.3;
%         else
%             TG(i) = TG(i) * 1.3;
%         end
%         if TG(i) > (rng(i)/3)
%             TG(i) = rng(i)/3;
%         end
%         if TG(i) < (rng(i)/300)
%             TG(i) = rng(i)/300;
%         end
%     end
    
%     parfor i = 1:PopSize                                        %不同父代基因
    for i = 1:PopSize  
        %为并行运算做的修改
        Pop = pop;
        Rng =rng;
        Bound = bound;    
  
%         temp = [];                                            %
%         temp = ones(q,1)*Pop(i,:);                          %克隆抗体
        for j =1:q                                          %变异循环
%             ifmutate = 0;                                   %变异标识
            tempmutate = pop(i,1:numvar);
            for k=1:numvar
                                       
                mutate = Pop(i,k) + randn * TG(k);          %抗体变异 
%                 length = abs(mutate-Pop(i,k))/rng(k);       %变异基因与原基因相对欧几里得距离
%                 pg=pK(i)^length;                            %变异概率
%                 if rand<pg
% %                 mutate = random('norm',Pop(i,k),TG(i));          %抗体变异
                if (mutate>Bound(k,1) && mutate < Bound(k,1) + Rng(k))
                    if rand<pK(i)
%                         temp = [temp, mutate];
                        tempmutate(k) = mutate;
%                         ifmutate = 1;                           %发生变异后改变变异标识                    
                    end
                end
            end
            if tempmutate == Pop(i,1:numvar)
            
            else
                b = [b; tempmutate]; 
%                 temp(j,numvar+1) = FbgFunc(temp(j,:));
            end
        end
        mutateQ(i) = size(b,1);
%         b = [b;temp];
    end
    [b, mutateQ] = getFbgb(b, mutateQ(PopSize), mutateQ);
    b = [b, ones(size(b,1),1)*TG];
    if mutateQ(1)==0
        b(1:mutateQ(1),numvar+2:numvar+6) = b(1:mutateQ(1),numvar+2:numvar+6).*1.3;
    else
        if mean(b(1:mutateQ(1),numvar+1)) < pop(1,numvar+1)
            b(1:mutateQ(1),numvar+2:numvar+6) = b(1:mutateQ(1),numvar+2:numvar+6)./1.3;
        end
        if mean(b(1:mutateQ(1),numvar+1)) > pop(1,numvar+1)
            b(1:mutateQ(1),numvar+2:numvar+6) = b(1:mutateQ(1),numvar+2:numvar+6).*1.3;
        end
    end
    for i = 2:PopSize
        if mutateQ(i)==0
            b(mutateQ(i-1)+1:mutateQ(i),numvar+2:numvar+6) = b(mutateQ(i-1)+1:mutateQ(i),numvar+2:numvar+6).*1.3;
        else
            if mean(b(mutateQ(i-1)+1:mutateQ(i),numvar+1)) < pop(i,numvar+1)
                b(mutateQ(i-1)+1:mutateQ(i),numvar+2:numvar+6) = b(mutateQ(i-1)+1:mutateQ(i),numvar+2:numvar+6)./1.3;
            end
            if mean(b(mutateQ(i-1)+1:mutateQ(i),numvar+1)) > pop(i,numvar+1)
                b(mutateQ(i-1)+1:mutateQ(i),numvar+2:numvar+6) = b(mutateQ(i-1)+1:mutateQ(i),numvar+2:numvar+6).*1.3;
            end
        end
    end
end
