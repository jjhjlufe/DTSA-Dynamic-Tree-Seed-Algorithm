clear all;
 %mex cec14_func.cpp -DWINDOWS
func_num=1;
D=100; %维度开关
Xmin=-100; %最小值开关
Xmax=100; % 最大值开关
pop_size=30; %总体个数
iter_max=500; %最大迭代次数  
runs=2; %一般为30
fhd=str2func('cec14_func');

convergenceData = cell(30, 1);


for i=1:30%测试哪个函数开关
    func_num=i; 
    for j=1:runs %每个函数测试多少趟开关，一般为30趟

        [TSA(i,j).gbest,TSA_gbestVal(i,j),TSA(i,j).convergence,TSA(i,j).FES]= TSA_func(fhd,D,pop_size,iter_max,Xmin,Xmax,func_num);=
        [DE_TSA2(i,j).gbest,DE_TSA_gbestVal2(i,j),DE_TSA2(i,j).convergence, DE_TSAtime2(i,j)]  = DTSA_func(fhd,D,pop_size,iter_max,Xmin,Xmax,func_num);

        convergenceData{i}{j, 1} = TSA(i,j).convergence;
        convergenceData{i}{j, 6} = DE_TSA2(i,j).convergence;
    end

 
     f_mean_TSA(i)=mean(TSA_gbestVal(i,:)); 
      TSA(i).convergMean=mean(cat(1,TSA(i,:).convergence)); %记录TSA收敛曲线
       f_std_TSA(i)=std(TSA_gbestVal(i,:));
      
     f_mean_DE_TSA2(i)=mean(DE_TSA_gbestVal2(i,:)); 
      DE_TSA2(i).convergMean=mean(cat(1,DE_TSA2(i,:).convergence)); %记录TSA收敛曲线
     f_std_DE_TSA2(i)=std(DE_TSA_gbestVal2(i,:)); 
     f_mean_DE_TSAtime2(i)=mean(DE_TSAtime2(i,:));

  figure(i)


    % 绘制数据，使用红色 ('r') 进行对数坐标绘图
    semilogy(TSA(i).convergMean, 'b');
    hold on
    semilogy(DE_TSA2(i).convergMean,'Color', 'r');
    hold on
    indices = 1:50:500;
    plot(indices, DE_TSA2(i).convergMean(indices), 'r*'); % 标记红色圆点

    legend('TSA','DTSA');
    title('Convergence curve');
    xlabel('Iteration');
    ylabel('Best score obtained so far'); 

end

     




% for i=1:30
% eval(['load input_data/shift_data_' num2str(i) '.txt']);
% eval(['O=shift_data_' num2str(i) '(1:10);']);
% f(i)=cec14_func(O',i);
% end