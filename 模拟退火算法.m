clc,clear            %清空环境中的变量
tic
iter = 1;                                                                                  % 迭代次数初值
a=0.99;                                                                                    %温度衰减系数
t0=120;                                                                                    %初始温度
tf=1;                                                                                      %最后温度
t=t0;
Markov=10000;                                                                              %Markov链长度
load cities.txt                                                                            %读入城市的坐标
city=cities;
n = size(city,1);                                                                          %城市距离初始化
D = zeros(n,n);                                                    
for i = 1:n
    for j = 1:n
            D(i,j) = sqrt(sum((city(i,:) - city(j,:)).^2));
    end    
end                                                                                
route=1:n;   
route_new=route;
best_length=Inf;
Length=Inf;
best_route=route;
%%
while t>=tf
            for j=1:Markov
                    %进行扰动,长生新的序列route_new;
                    if (rand<0.7)
                        %交换两个数的顺序
                           ind1=0;ind2=0;
                           while(ind1==ind2&&ind1>=ind2)
                                    ind1=ceil(rand*n);
                                    ind2=ceil(rand*n);
                           end                      
                                      temp=route_new(ind1);
                                      route_new(ind1)=route_new(ind2);
                                      route_new(ind2)=temp;
                    else
                          ind=zeros(3,1);
                          L_ind=length(unique(ind));
                          while (L_ind<3)
                                    ind=ceil([rand*n rand*n rand*n]);
                                    L_ind=length(unique(ind));
                          end
                          ind0=sort(ind);
                          a1=ind0(1);b1=ind0(2);c1=ind0(3);
                         route0=route_new;
                         route0(a1:a1+c1-b1-1)=route_new(b1+1:c1);
                         route0(a1+c1-b1:c1)=route_new(a1:b1);
                         route_new=route0;    
                    end 
                     %计算路径的距离,Length_new 
                          length_new = 0;
                          Route=[route_new route_new(1)];
                              for j = 1:n
                                  length_new = length_new+ D(Route(j),Route(j + 1));
                              end
                     if length_new<Length      
                              Length=length_new;
                              route=route_new;
                           %对最优路线和距离更新
                           if       length_new<best_length
                                    iter = iter + 1;    
                                     best_length=length_new;
                                     best_route=route_new;
                           end
                     else
                             if rand<exp(-(length_new-Length)/t)
                                  route=route_new;
                                  Length=length_new;
                              end
                     end
                       route_new=route; 
                end              
                        t=t*a;
end

%% 结果显示 
toc
Route=[best_route best_route(1)];
plot([city(Route ,1)], [city(Route ,2)],'o-');
    disp('最优解为：')
    disp(best_route)
    disp('最短距离：')
    disp(best_length)
    disp('最优解迭代次数：')
    disp(iter)
for i = 1:n
    %对每个城市进行标号
    text(city(i,1),city(i,2),['   ' num2str(i)]);
end
xlabel('城市位置横坐标')
ylabel('城市位置纵坐标')
title(['模拟退火算法(最短距离）:' num2str(best_length) ''])