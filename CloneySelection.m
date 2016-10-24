%PACS算法 并行自适应克隆选择算法
%By dxb 201501013

clear;
clc;
parpool(6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%初始化%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bound = [0 1;0 1;0 30; 0 50; 0 30];         %变量取值范围
MAXGEN = 200;                               %迭代次数
nc = 80;                                    %克隆规模 大于种群数量
pm = 0.5;                                   %抗体变异概率
popSize = 40;                               %种群规模
numVar = 5;                                 %变量个数
d = 4;                                      %替换抗体数量
newBorn = 20;                               %随机生成新子代
Fs = 20;                                    %选择池大小
[pop,rng] = initializeCS(popSize, bound);   %随机生成初始总群
TG = rng/30;                                %变异系数          3sigma
pop = getFbg(pop, popSize);                 %获取亲合度
pop = [pop, ones(size(pop,1),1)*TG];        %将变异系数加入到总群中
result = zeros(MAXGEN, numVar*2+3);          %设置结果存储空间
Count = 1;                                  %初始化迭代次数
popNew = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%初始化完成%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while ((Count < MAXGEN))                    %终止条件判断
    %选择
    pop = sortrows(pop,numVar+1);           %按亲合度排序
    if size(pop,1)<Fs
        popFs = pop;
    else
        popFs = pop(1:Fs,:);                    %选择亲合度好的抗体
    end
    popFs = [popFs; popNew];                %加入选择后的新抗体
    
    %计算克隆数量
    q=ceil((sum(popFs(:,numVar+1))-popFs(:,numVar+1))./((Fs+d-1)*sum(popFs(:,numVar+1))).*nc);

    %基因操作(克隆+变异)
    [b, TG] = geneOps3(popFs, q, numVar, rng, bound, pm, TG); 
    pop = [pop; b];
    
    %再选择                                             
    pop = sortrows(pop,numVar+1);                   %按亲合度排序
    pop = pop(1:popSize,:);                         %淘汰亲合度差的抗体
    
    %reborn
    [popNew,rngTmp] = initializeCS(newBorn, bound); %随机生成新抗体
    popNew = getFbg(popNew, newBorn);               %计算新抗体亲合度
    popNew = [popNew, ones(size(popNew,1),1)*TG];   %初始化新抗体变异系数
    popNew = sortrows(popNew,numVar+1);             %按亲合度排序
    if size(popNew,1)>d        
        popNew = popNew(1:d,:);                         %选择最好的新抗体  
    end
             
    %结果整理
    result(Count,:)=[pop(1,:), Fs, nc];  
    
    %数据显示
    pop(:,numVar+1)
    result(Count,1:numVar)
    result(Count, numVar+1)
    Count = Count + 1
end

plot(1:MAXGEN-1,result(1:MAXGEN-1, numVar+1))
delete(gcp)