function f = genetic_algorithm(D, mute, Pm)
% 遗传算法解TSP问题的简单实现
% 注释我就不写了 很简单能看懂。另外一些细节看MD。

if ~exist('mute', 'var')
    mute = 0; % 是否显示各种提示信息
end
if ~exist('Pm', 'var')
    Pm = 0.2; % 变异概率，越大收敛越慢但是解一般越好
end
if ~exist('D', 'var')
    D= [0, 1, 2, 3, 4, 5;
        1, 0, 6, 7, 8, 9;
        2, 6, 0, 8, 7, 6;
        3, 7, 8, 0, 5, 4;
        4, 8, 7, 5, 0, 3;
        5, 9, 6, 4, 3, 0]; % D为城市间的距离矩阵；可以不对称
end
rng(1);
n = size(D, 1);
N = 100; % 群体规模
TOL = 20; % 最大容忍次数(连续TOL次rate不上升，或找不到更优解，则停止迭代)

solutions = zeros(N, n);
fs = zeros(N, 1);
for i = 1:N
    solutions(i, :) = [1, randperm(n-1) + 1]; % 生成N个解，假定从1开始
    fs(i) = TSP_distance(D, solutions(i, :));
end
Pu = max(fs) - fs + 1;
P = Pu/sum(Pu);
cumP = cumsum(P);
best = min(fs);
avg = mean(fs);
rate = best/avg;
if ~mute
disp('初始解的群体中最短的路径长度为：');
disp(best);
disp('初始解的群体中平均路径长度为：');
disp(avg);
end

tol = 0;
count = 0;
while 1
    count = count + 1;
    if ~mute
        fprintf('当前第%d次迭代\n', count);
    end
    parents = zeros(size(solutions));
    for i = 1:N % 使用轮盘赌的方式选出父代的染色体
        index = sum(cumP <= rand) + 1;
        parents(i, :) = solutions(index, :);
    end
    new_solutions = zeros(size(solutions));
    assert(mod(N, 2) == 0);
    for i = 1:N/2 % 交配操作；这里默认N为偶数，每两个父代一起产生两个子代
        % 产生的子代1取父代1的前一半染色体，后一半则由父代2提供；同理于子代2
        p1 = parents(2*i-1, :);
        p2 = parents(2*i, :);
        middle = ceil(n/2);
        s1 = p1(1:middle);
        res1 = setdiff(p2, s1, 'stable');
        s1 = [s1, res1];
        s2 = p2(1:middle);
        res2 = setdiff(p1, s2, 'stable');
        s2 = [s2, res2];
        new_solutions(2*i-1, :) = s1;
        new_solutions(2*i, :) = s2;
    end
    for i = 1:N % 变异操作；变异的方式为随机取两个城市的顺序交换
        if rand < Pm
            temp = randperm(n-1) + 1;
            k = temp(1);
            new_solutions(i, [1, k]) = new_solutions(i, [k, 1]);
        end
    end
    % 至此，新的种群已经生成完毕，代替旧种群后开始新一轮的计算
    solutions = new_solutions;
    for i = 1:N
        fs(i) = TSP_distance(D, solutions(i, :));
    end
    Pu = max(fs) - fs + 1;
    P = Pu/sum(Pu);
    cumP = cumsum(P);
    best_new = min(fs);
    avg = mean(fs);
    rate_new = best_new/avg;
    if ~mute
        disp('最短的路径长度为：');
        disp(best_new);
        disp('平均的路径长度为：');
        disp(avg);
    end
    tol = tol + 1;
    if best_new < best || rate_new > rate
        best = best_new;
        rate = rate_new;
        tol = 0;
    end
    if tol >= TOL
        break
    end
    if count > 5000
        break
    end
end
[f, index] = min(fs);
solution = solutions(index, :);
if ~mute
fprintf('最后搜索得到的最优路径为：\n');
disp(solution);
disp('路径长度为：');
disp(f);
end

function f = TSP_distance(D, solution)
% 本函数计算给定solution的距离，其中城市之间的距离由D给出。
n = numel(solution);
sum = 0;
for i = 1:n-1
    sum = sum + D(solution(i), solution(i+1));
end
sum = sum + D(solution(n), solution(1));
f = sum;