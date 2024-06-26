function  [gbest,gbestval,convergence,fitcount]= DE_TSA_func2(fhd,Dimension,Particle_Number,Max_Gen,VRmin,VRmax,varargin)
      
    % ��������
    N =Particle_Number;             % ��ʼ��������
    low = ceil(0.1 * N);   % ��������������
    high = ceil(0.25 * N);  % ��������������
    D = Dimension;              % ����ά��
    ST = 0.1;          % �����������������ѡ�񷽳�
    dmin = VRmin;        % �����ռ��½�
    dmax =VRmax;         % �����ռ��Ͻ�

    fes=N;              % is the fes counter.
    max_fes=500000;    % is the termination condition, ��PSO�ĺ��������������
  
    trees = zeros(N, D);    % ����Ⱥ
    obj = zeros(1, N);      % ÿ������Ŀ�꺯��ֵ
    vmax =2;
  
    count = 0;

    t = 1;
    convergence = ones(1,Max_Gen);
  
for i = 1:N
    % ʹ��Bine����ӳ�������������bine)
    chebyshev_x = chebyshev_map(D, 1);
     for j = 1:D
        v(i,j)=vmax*chebyshev_x(j);
        trees(i, j) = dmin + (dmax - dmin) *chebyshev_x(j);   % ���������ռ��λ��
    end
end
    obj = feval(fhd, trees',varargin{:});

    % Ѱ�ҵ�ǰ�������λ��
    [~, indis] = min(obj);
    indis = indis(end);
    best = obj(indis); % ��¼���λ�õĺ���ֵ
    best_params = trees(indis, :);  % �����λ�õ����Խ��и�ֵ

        while(t<= Max_Gen)
            
                %Ȩ������Ӧ
                w = 0.8 * exp((-2.5 * t) / Max_Gen);

                k = 2 -2 * (t /Max_Gen);

                %  ÿһ��������������
                for i = 1:N   
                    ns = low + randi(high - low + 1) - 1;
                    ns = ns(1);
                    seeds = zeros(ns, D);
                    objs = zeros(1, ns);
                    vs = ones(ns+1,D);
                    % ʹ��ѭ���� vs ��ÿһ������Ϊ v �Ķ�Ӧ�е�ֵ
                    for u = 1:size(vs, 1) % �� vs ��ÿһ�н��е���
                        vs(u, :) = v(i, :); % �� vs �ĵ� i ������Ϊ v �ĵ� i �е�ֵ
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

                                    seeds(j,d) =  trees(i, d) +( (best_params(d) - trees(r, d)) * (rand - 0.5) * 2  )  +vs(j+1,d);%����λ��     
                        

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
                       

                    end %����ѭ������

                    
                    
                    
                      objs=feval(fhd,seeds',varargin{:});

                    

                % Ѱ�ҵ�ǰ������ӵ�λ��
                    [~,seed_indis] = min(objs);
                   seed_indis = seed_indis(end);
                    best_seed = objs(seed_indis); % ��¼���λ�õĺ���ֵ
                    best_seed_params = seeds(seed_indis,:);

                 % ��������λ��
                if(best_seed < obj(i))
                    trees(i, :) = best_seed_params;
                    v(i,:)=vs(seed_indis,:);
                    obj(i) =best_seed;
                end
                  fes=fes+ns;       


                % ���µ�ǰ�������λ��
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

               
             
                end %ÿ������������

             if (best == prev_best) || isequal(prev_best_params, best_params)
                    count = count +1;
                end
            
            
            % ����Ӧֵ����Ϊ���벻������Ⱥ��
        
            
           
           
            % ������Ӧ��ֵ��С���������������
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

                   
            % �������
             
            for i = 1:size(badGroup1,1)
                rand_num = rand();
                badGroup1(i,:) = rand_num .* betterGroup2(i,:) + (1-rand_num) .* badGroup1(i,:);
            end

            for i = 1:size(badGroup1,1)
                rand_num = rand();
                badGroup1_v(i,:) = rand_num .* betterGroup2_v(i,:) + (1-rand_num) .* badGroup1_v(i,:);
            end
            



            % ��Ȼѡ�����            
            badGroup2(:,1:D) = betterGroup1(:,1:D);
            badGroup2_v(:) = betterGroup1_v(:);
            

            % ���»�����λ�ü��ٶ�
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
            


            % ��¼convergence���
            convergence(t)=best;
             t =t + 1;
        end
    gbest = best_params; %�������ֵ
    gbestval =best; %�������λ��
     fitcount=fes; %���غ�����������
 
end