function [ pop, mutateQ] = getFbgb( pop, popSize, mutateQ )
%getFbg ����һ�鿹��Ŀ���-��ԭ�׺϶�
%   pop - ���鿹��Ļ�����Ϊĳ����Ļ���Σ���Ϊ��ͬ����
%   popSize - ���鿹��ĸ���
%   Fbgf - �������鿹���Ӧ�Ŀ���-��ԭ�׺϶�
%by dxb 201501013
    mutateQSize = size(mutateQ,1);
    numvar = size(pop, 2);
    Fbgf = zeros(popSize,1);            %����Fbgf�洢�ռ�
    parfor i=1:popSize
%     for i=1:popSize
        Fbgf(i) = FbgFunc(pop(i,:));
    end
    pop = [pop, Fbgf];
    i = 1;
    k = 1;
    while i<=popSize
        if pop(i,numvar+1) == 1e+14
            pop(i,:) = [];              %������Ĳ����޳�
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

