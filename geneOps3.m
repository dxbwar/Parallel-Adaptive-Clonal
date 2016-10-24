function [b, TG] = geneOps3( pop, q, numvar, rng, bound, pm, TG)
%geneOps �������
%   ���������¡�󰴸��ʽ��б���
%   pop - ���鿹��Ļ�����Ϊĳ����Ļ���Σ���Ϊ��ͬ����
%   popSize - ���鿹��ĸ���
%   q - ԭ�������Ӧ�Ŀ�¡����
%   numvar - �����������������
%   rng - ����ȡֵ��Χ
%   bound - ����ȡֵ�߽�
%   pm - ÿ�λ����Ӧ�ı������
%   b - ���ر����Ŀ���-��ԭ�׺������ſ�����
%by dxb 20140807
    PopSize = size(pop, 1); 
    b = [];                                                    %��¡������ſ��弯��
    mutateQ = zeros(PopSize,1);                                 %ÿ�������¡��������ĸ���
%     b = [b,zeros(size(b,1)];                                    %��������ʶ
    Fbg = pop(:,numvar+1);                                      %��ȡ�׺϶�
    avrFbg = mean(Fbg);
    pK = pm + ((1-pm)^-1-1+exp(avrFbg./Fbg)).^-1;                         %���������ʼ���
    
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
    
%     parfor i = 1:PopSize                                        %��ͬ��������
    for i = 1:PopSize  
        %Ϊ�������������޸�
        Pop = pop;
        Rng =rng;
        Bound = bound;    
  
%         temp = [];                                            %
%         temp = ones(q,1)*Pop(i,:);                          %��¡����
        for j =1:q                                          %����ѭ��
%             ifmutate = 0;                                   %�����ʶ
            tempmutate = pop(i,1:numvar);
            for k=1:numvar
                                       
                mutate = Pop(i,k) + randn * TG(k);          %������� 
%                 length = abs(mutate-Pop(i,k))/rng(k);       %���������ԭ�������ŷ����þ���
%                 pg=pK(i)^length;                            %�������
%                 if rand<pg
% %                 mutate = random('norm',Pop(i,k),TG(i));          %�������
                if (mutate>Bound(k,1) && mutate < Bound(k,1) + Rng(k))
                    if rand<pK(i)
%                         temp = [temp, mutate];
                        tempmutate(k) = mutate;
%                         ifmutate = 1;                           %���������ı�����ʶ                    
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
