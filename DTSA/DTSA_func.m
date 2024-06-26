function  [gbest,gbestval,convergence,fitcount]= DE_TSA_func2(fhd,Dimension,Particle_Number,Max_Gen,VRmin,VRmax,varargin)
      
    % 参数设置
    N =Particle_Number;             % 初始树的数量
    low = ceil(0.1 * N);   % 种子数量的下限
    high = ceil(0.25 * N);  % 种子数量的上限
    D = Dimension;              % 问题维度
    ST = 0.1;          % 搜索倾向参数，用于选择方程
    dmin = VRmin;        % 搜索空间下界
    dmax =VRmax;         % 搜索空间上界

    fes=N;              % is the fes counter.
    max_fes=500000;    % is the termination condition, 与PSO的函数评估次数相等
  
    trees = zeros(N, D);    % 树种群
    obj = zeros(1, N);      % 每棵树的目标函数值
    vmax =2;
  
    count = 0;

    t = 1;
    convergence = ones(1,Max_Gen);
  
for i = 1:N
    % 使用Bine混沌映射生成随机数（bine)
    chebyshev_x = chebyshev_map(D, 1);
     for j = 1:D
        v(i,j)=vmax*chebyshev_x(j);
        trees(i, j) = dmin + (dmax - dmin) *chebyshev_x(j);   % 树在搜索空间的位置
    end
end
    obj = feval(fhd, trees',varargin{:});

    % 寻找当前最佳树的位置
    [~, indis] = min(obj);
    indis = indis(end);
    best = obj(indis); % 记录最佳位置的函数值
    best_params = trees(indis, :);  % 将最佳位置的属性进行赋值

        while(t<= Max_Gen)
            
                %权重自适应
                w = 0.8 * exp((-2.5 * t) / Max_Gen);

                k = 2 -2 * (t /Max_Gen);

                %  每一颗树的种子生成
                for i = 1:N   
                    ns = low + randi(high - low + 1) - 1;
                    ns = ns(1);
                    seeds = zeros(ns, D);
                    objs = zeros(1, ns);
                    vs = ones(ns+1,D);
                    % 使用循环将 vs 的每一行设置为 v 的对应行的值
                    for u = 1:size(vs, 1) % 对 vs 的每一行进行迭代
                        vs(u, :) = v(i, :); % 将 vs 的第 i 行设置为 v 的第 i 行的值
                    end
                   

                    for j = 1:ns
                        % Determination of the neighbour for ith tree.    
                             r = fix(rand * N) + 1;
                            while(i == r)  
                                r = fix(rand * N) + 1;
                            end  
                  

                        %%% SEEDS ARE CREATING
                        for d = 1:D
                             
                                if(rand < ST)

                                     if count >=15
                                         vs(j+1,d) = - (w * vs(j,d)  +rand  *( best_params(d) -  trees(i, d) ) * k * cos(2*pi*rand));
                                          
                                     else
                                          vs(j+1,d) = w * vs(j,d)  +rand  *( best_params(d) -  trees(i, d) ) * k * cos(2*pi*rand);
                                     end
                                     
                                    if abs(vs(j,d))>vmax
                                        if vs(j,d)<0
                                            vs(j,d)=-vmax;
                                        end
                                        if vs(j,d)>0
                                            vs(j,d)=vmax;
                                        end
                                   end

                                    seeds(j,d) =  trees(i, d) +( (best_params(d) - trees(r, d)) * (rand - 0.5) * 2  )  +vs(j+1,d);%更新位置     
                        

                                else
                                    
                                                            

                                    if count >=15
                                          vs(j+1,d) =-( w * vs(j,d)  +rand  *( trees(r,d)   - trees(i, d) ) * k * cos(2*pi*rand));
                                          
                                    else
                                         vs(j+1,d) = w * vs(j,d)  +rand  *( trees(r,d)   - trees(i, d) ) * k * cos(2*pi*rand); 
                                    end
                                    if abs(vs(j,d))>vmax
                                        if vs(j,d)<0
                                            vs(j,d)=-vmax;
                                        end
                                        if vs(j,d)>0
                                            vs(j,d)=vmax;
                                        end
                                   end           
                                      seeds(j,d) = trees(i, d) +( (trees(r, d) -  trees(i, d) ) *(rand - 0.5) * 2)  + vs(j+1,d);


                                    
         
                                end
                                
                                
                                
                                

                                    if(seeds(j, d) > dmax || seeds(j, d) < dmin)
                                        seeds(j, d) = dmin + (dmax - dmin) * rand;
                                    end
                                
                        end
                       

                    end %种子循环结束

                    
                    
                    
                      objs=feval(fhd,seeds',varargin{:});

                    

                % 寻找当前最佳种子的位置
                    [~,seed_indis] = min(objs);
                   seed_indis = seed_indis(end);
                    best_seed = objs(seed_indis); % 记录最佳位置的函数值
                    best_seed_params = seeds(seed_indis,:);

                 % 更新树的位置
                if(best_seed < obj(i))
                    trees(i, :) = best_seed_params;
                    v(i,:)=vs(seed_indis,:);
                    obj(i) =best_seed;
                end
                  fes=fes+ns;       


                % 更新当前最佳树的位置
                [~, indis] = min(obj);
                indis = indis(end);
                temp_best = obj(indis);
               
                prev_best = best;
                prev_best_params  =  trees(indis, :);

                if(temp_best < best)
                    best = temp_best;
                    best_params = trees(indis, :);
                    count = 0;
                end

               
             
                end %每棵树迭代结束

             if (best == prev_best) || isequal(prev_best_params, best_params)
                    count = count +1;
                end
            
            
            % 将适应值划分为好与不好两个群体
        
            
           
           
            % 按照适应度值从小到大对树进行排序
            [~, sort_idx] = sort(obj);
            sort_trees = trees(sort_idx,:);
            sort_v = v(sort_idx,:);

            betterGroup1_idx = sort_idx(1:ceil(0.1*N));
            betterGroup2_idx = sort_idx(ceil(0.1*N)+1 : floor(0.5*N));
            badGroup1_idx = sort_idx(floor(0.5*N)+1 : floor(0.9*N));
            badGroup2_idx = sort_idx(floor(0.9*N)+1 : N);
            

            
            
            betterGroup1_v = sort_v(betterGroup1_idx,:);
            betterGroup2_v = sort_v(betterGroup2_idx,:);
            badGroup1_v = sort_v(badGroup1_idx,:);
            badGroup2_v = sort_v(badGroup2_idx,:);
        
            betterGroup1 = sort_trees(betterGroup1_idx,:);
            betterGroup2 = sort_trees(betterGroup2_idx,:);
            badGroup1 = sort_trees(badGroup1_idx,:);
            badGroup2 = sort_trees(badGroup2_idx,:);

                   
            % 交叉操作
             
            for i = 1:size(badGroup1,1)
                rand_num = rand();
                badGroup1(i,:) = rand_num .* betterGroup2(i,:) + (1-rand_num) .* badGroup1(i,:);
            end

            for i = 1:size(badGroup1,1)
                rand_num = rand();
                badGroup1_v(i,:) = rand_num .* betterGroup2_v(i,:) + (1-rand_num) .* badGroup1_v(i,:);
            end
            



            % 自然选择操作            
            badGroup2(:,1:D) = betterGroup1(:,1:D);
            badGroup2_v(:) = betterGroup1_v(:);
            

            % 更新坏的树位置及速度
            trees(badGroup1_idx,:) = badGroup1;
            trees(badGroup2_idx,:) = badGroup2;
            
            for i=1:size(badGroup1_v)
                num_v = badGroup1_idx(i);
                v(num_v,:) = badGroup1_v(i,:);
            end
            for i=1:size(badGroup2_v)
                  num_v = badGroup2_idx(i);
                v(num_v,:) = badGroup2_v(i,:);
            end
            


            % 记录convergence情况
            convergence(t)=best;
             t =t + 1;
        end
    gbest = best_params; %返回最佳值
    gbestval =best; %返回最佳位置
     fitcount=fes; %返回函数评估次数
 
end