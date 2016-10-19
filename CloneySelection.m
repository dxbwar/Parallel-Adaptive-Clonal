%PACS�㷨 ��������Ӧ��¡ѡ���㷨
%By dxb 201501013

clear;
clc;
parpool(6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%��ʼ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bound = [0 1;0 1;0 30; 0 50; 0 30];         %����ȡֵ��Χ
MAXGEN = 200;                               %��������
nc = 80;                                    %��¡��ģ ������Ⱥ����
pm = 0.5;                                   %����������
popSize = 40;                               %��Ⱥ��ģ
numVar = 5;                                 %��������
d = 4;                                      %�滻��������
newBorn = 20;                               %����������Ӵ�
Fs = 20;                                    %ѡ��ش�С
[pop,rng] = initializeCS(popSize, bound);   %������ɳ�ʼ��Ⱥ
TG = rng/30;                                %����ϵ��          3sigma
pop = getFbg(pop, popSize);                 %��ȡ�׺϶�
pop = [pop, ones(size(pop,1),1)*TG];        %������ϵ�����뵽��Ⱥ��
result = zeros(MAXGEN, numVar*2+3);          %���ý���洢�ռ�
Count = 1;                                  %��ʼ����������
popNew = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%��ʼ�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while ((Count < MAXGEN))                    %��ֹ�����ж�
    %ѡ��
    pop = sortrows(pop,numVar+1);           %���׺϶�����
    if size(pop,1)<Fs
        popFs = pop;
    else
        popFs = pop(1:Fs,:);                    %ѡ���׺϶ȺõĿ���
    end
    popFs = [popFs; popNew];                %����ѡ�����¿���
    
    %�����¡����
    q=ceil((sum(popFs(:,numVar+1))-popFs(:,numVar+1))./((Fs+d-1)*sum(popFs(:,numVar+1))).*nc);

    %�������(��¡+����)
    [b, TG] = geneOps3(popFs, q, numVar, rng, bound, pm, TG); 
    pop = [pop; b];
    
    %��ѡ��                                             
    pop = sortrows(pop,numVar+1);                   %���׺϶�����
    pop = pop(1:popSize,:);                         %��̭�׺϶Ȳ�Ŀ���
    
    %reborn
    [popNew,rngTmp] = initializeCS(newBorn, bound); %��������¿���
    popNew = getFbg(popNew, newBorn);               %�����¿����׺϶�
    popNew = [popNew, ones(size(popNew,1),1)*TG];   %��ʼ���¿������ϵ��
    popNew = sortrows(popNew,numVar+1);             %���׺϶�����
    if size(popNew,1)>d        
        popNew = popNew(1:d,:);                         %ѡ����õ��¿���  
    end
             
    %�������
    result(Count,:)=[pop(1,:), Fs, nc];  
    
    %������ʾ
    pop(:,numVar+1)
    result(Count,1:numVar)
    result(Count, numVar+1)
    Count = Count + 1
end

plot(1:MAXGEN-1,result(1:MAXGEN-1, numVar+1))
delete(gcp)