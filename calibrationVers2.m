global input
input = csvread('dataform2018.csv');
alpha = [0 1];
prop = [0.8 0.2];

group = randsrc(256,90,[alpha;prop]);  
patience_gen = 3;    %break the loop

while true
    tmp = zeros(256,91);
    tmp_1 = tmp;
    for k = 1:128
        a=normrnd(0,100,1,2);
        a1=1+mod(abs(fix((a(1)))),90);
        a2=1+mod(abs(fix((a(2)))),90);
        tmp(2*k-1:2*k,1:90) = exchange(group(a1,:),group(a2,:));%crossing over and produce a new generation
    end
    %each individual mutate randomly according to their level of fitness
    for j = 1:256
        tmp(j,1:90)=mutate(tmp(j,1:90),4);
    end
    for k = 1:256
        tmp(k,91) = fitness(tmp(k,1:90));
    end
    tmp = choose(tmp);
    tmp = sortrows(tmp,+91);%sort the group according to their fitness descending
    group = tmp(1:256,1:90);
    disp(tmp(1,91));
end

out(tmp,gen);
function new = choose(old)
dim = size(old);
new = zeros(dim(1), dim(2));
for k = 1:dim(1)
    rad = randi([1 256],1,2);
    val_1 = old(rad(1),91);
    val_2 = old(rad(2),91);
    if val_1<val_2
        m = rad(1);
    else m = rad(2);
    end
    new(k,:) = old(m,:);
end
end
    
 
function [] = out(new,g)
disp(new(1,91));
disp(g);
end

%exchange a random section from 2 individuals
function children = exchange(m1, m2)
x = randi(7);
cut = randi(89,1,x);
temp = m1(1:cut(1));
m1(1:cut(1)) = m2(1:cut(1));
m2(1:cut(1)) = temp;
for k = 1 : x - 1
    temp = m1(cut(k):cut(k+1));
    m1(cut(k):cut(k+1)) = m2(cut(k):cut(k+1));
    m2(cut(k):cut(k+1)) = temp;
end
temp = m1(cut(x):90);
m1(cut(x):90) = m2(cut(x):90);
m2(cut(x):90) = temp;
children=[m1;m2];
end

%individuals mutant randomly
function output = mutate(m,p)
temp = [rand(1,p);randi(90,1,p)];
for i =1:p
    if temp(1,i) > 0.95
        m(temp(2,i)) = 1 - m(temp(2,i));        
    end
end
output=m;
end

function s = fitness(v)    
global input
measure = ones(1,90)*50;  
cost_measure = dot(v,measure);
cost_error = 0;
for k = 1:500
    x = v.*input(2*k-1,:);
    y = v.*input(2*k,:);
        y(y==0) = [];
    if v(21)==0
        x(x==0) = [];
    else
        x1 = x(1,1:20);
        x2 = x(1,22:90);
        x1(x1==0) = [];
        x2(x2==0) = [];
        x = [x1 0 x2];
    end
    estimate = spline(x, y, input(2*k-1,:));
    error = abs(estimate-input(2*k,:));
    for m = 1:90
        a = judge(error(m));
        cost_error = cost_error+a;
    end
end
s = cost_measure+cost_error/500;
end

function s = judge(k)   
if k<=0.5
    s = 0;
elseif k<=1
    s = 1;
elseif k<=1.5
    s = 4;
elseif k<=2
    s = 10;
elseif k>2
    s = 10000;
end
end 