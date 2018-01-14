function f = simulated_annealing(D, mute, MAX_ITER)
% 模拟退火解TSP问题的简单实现
% 注释我就不写了 很简单能看懂。另外一些细节看MD。
if ~exist('mute', 'var')
    mute = 0; % 是否显示各种提示信息
end
if ~exist('MAX_ITER', 'var')
    MAX_ITER = 50; % 每个温度内的最大状态交换尝试次数
                   % 这个参数很玄乎没什么固定规律的样子
end
if ~exist('D', 'var')
    D= [0, 1, 2, 3, 4, 5;
        1, 0, 6, 7, 8, 9;
        2, 6, 0, 8, 7, 6;
        3, 7, 8, 0, 5, 4;
        4, 8, 7, 5, 0, 3;
        5, 9, 6, 4, 3, 0]; % 城市间的距离矩阵；可以不对称
end
rng(0);
n = size(D, 1);
T_range_factor = exp(0:-0.1:-5); % 温度的范围系数

solution = [1, randperm(n-1) + 1]; % 生成一个解，假定从1开始
P = generate_neighbors(solution); % 生成初始解的邻域P
if ~mute
fprintf('邻域内总共有%d个解。\n', numel(P));
end
f = TSP_distance(D, solution); % 求出当前解的总距离
if ~mute
disp('初始的路径为：');   
disp(solution);
disp('路径长度为：');
disp(f);
end
TMAX = f;
T_range = T_range_factor * TMAX; % 温度范围动态地随初始解的好坏而变化
                                 % 尽量保证开始的温度足够高 结束的足够低
first = zeros(1, MAX_ITER);
final = zeros(1, MAX_ITER);
for t = T_range
    for i = 1:MAX_ITER
        index = randi(numel(P), 1); % 产生一个1~|P|之间的随机整数
        neighbor = P{index}; % 在P中随机取一个解
        f_neighbor = TSP_distance(D, neighbor); % 计算这个解的距离
        Pt = exp(-(f_neighbor - f)/t); % 转移概率
        if Pt > 1, Pt = 1; end
        if Pt > rand
            f = f_neighbor;
            solution = neighbor;
            P = generate_neighbors(solution);
        end
        if ~mute && t == T_range(1)
            first(i) = Pt;
        end
        if ~mute && t == T_range(end)
            final(i) = Pt;
        end
    end
    if ~mute && t == T_range(1)
        fprintf('初始温度下的转移概率中位数为%f\n', median(first));
    end
    if ~mute && t == T_range(end)
        fprintf('最后温度下的转移概率中位数为%f\n', median(final));
    end
end
if ~mute
fprintf('最后搜索得到的最优路径为：\n');
disp(solution);
disp('路径长度为：');
disp(f);
end

function P = generate_neighbors(solution)
% 本函数根据交换任意两个城市的原则，产生solution的邻域
n = numel(solution);
P = cell(1, n*(n-1)/2);
t = 1;
for i = 1:n
    for j = i+1:n
        solution_new = solution;
        solution_new([i, j]) = solution_new([j, i]); % 交换两个元素
        P{t} = solution_new;
        t = t + 1;
    end
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