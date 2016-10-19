function [ pop, mutateQ] = getFbgb( pop, popSize, mutateQ )
%getFbg 计算一组抗体的抗体-抗原亲合度
%   pop - 这组抗体的基因，行为某抗体的基因段，列为不同抗体
%   popSize - 这组抗体的个数
%   Fbgf - 返回这组抗体对应的抗体-抗原亲合度
%by dxb 201501013
    mutateQSize = size(mutateQ,1);
    numvar = size(pop, 2);
    Fbgf = zeros(popSize,1);            %设置Fbgf存储空间
    parfor i=1:popSize
%     for i=1:popSize
        Fbgf(i) = FbgFunc(pop(i,:));
    end
    pop = [pop, Fbgf];
    i = 1;
    k = 1;
    while i<=popSize
        if pop(i,numvar+1) == 1e+14
            pop(i,:) = [];              %将出错的参数剔除
            popSize = popSize - 1;
            mutateQ(k:mutateQSize) = mutateQ(k:mutateQSize) - 1;
            if mutateQ(k)==0;
                k = k+1;
            end
        else
            i = i+1;
            if i>mutateQ(k)
                k = k+1;
            end
        end
    end
end

