%% 清空缓存
clc;
clear all;

%% 初始化参数

% DH参数
d2 = 160.5;
a3 = 612.0;
a4 = 590.0;
d5 = 156.5;
d6 = 125.0;
d7 = 108.5;
DH = [d2, a3, a4, d5, d6, d7];

degree = pi/180;

% 当前关节角度
a = -0 * degree;
b = -0 * degree;
c = -0 * degree;
d = 0 * degree;
e = 0 * degree;
f = 0 * degree;
CurQ = [a, b, c, d, e, f];

% 连杆角度
aa = 0 * degree;
ba = 0 * degree;
ca = 0 * degree;
da = 0 * degree;
ea = 0 * degree;
fa = 0 * degree;
Alpha = [aa,ba,ca,da,ea,fa]; 

% 计算当前末端姿态
CurT = ForwardKine(CurQ, Alpha, DH);
CurR = [CurT(1,1), CurT(1,2), CurT(1,3);
          CurT(2,1), CurT(2,2), CurT(2,3);
          CurT(3,1), CurT(3,2), CurT(3,3)];
Cureuler = rotm2eul(CurR, "ZYX");
CurP = [CurT(1,4),CurT(2,4),CurT(3,4)];
z = Cureuler(1);
y = Cureuler(2);
x = Cureuler(3);
Te = [cos(y)*cos(z),-sin(z), 0;
      cos(y)*sin(z), cos(z), 0;
      -sin(y), 0,     1];

% 计算目标末端姿态
a = -90 * degree;
b = 90 * degree;
c = 90 * degree;
d = 0 * degree;
e = 0 * degree;
f = 0 * degree;
[-1.57079632679490	-3.10565188673572	-1.62970441239366	-1.23280319506608	-5.68434188608080e-14	-0.360960450473385]
TarQ = [a, b, c, d, e, f]; 
TarT = ForwardKine(TarQ, Alpha, DH)
TarR = [TarT(1,1), TarT(1,2), TarT(1,3);
          TarT(2,1), TarT(2,2), TarT(2,3);
          TarT(3,1), TarT(3,2), TarT(3,3)];
Tareuler = rotm2eul(TarR, "ZYX");
TarP = [TarT(1,4), TarT(2,4), TarT(3,4)];
% [
%    1.0e+03 *
% 
%     0.0001   -0.0002   -0.0010   -0.4972
%    -0.0004    0.0009   -0.0003    0.8060
%     0.0009    0.0005   -0.0000    1.0387
%          0         0         0    0.0010]
%% 进入牛顿差值迭代，计算关节角度
for i=1:1000
deltaP = [TarP(1)-CurP(1),TarP(2)-CurP(2),TarP(3)-CurP(3)];
deltaEuler = [Tareuler(3)-Cureuler(3);Tareuler(2)-Cureuler(2);Tareuler(1)-Cureuler(1)];

rotational = Te * deltaEuler;
Fq = [deltaP(1), deltaP(2), deltaP(3), rotational(1), rotational(2), rotational(3)];

J = Jacobian_Forward(CurQ, Alpha, DH);

Jcross = pinv(J'*J)*J';

% 起点或终点过奇异点时
deltaQ = GetDeltaQ(Fq, J);

% 起点和终点不过奇异点时
% deltaQ = Jcross * Fq'
CurQ = CurQ + deltaQ';

%更新当前位置
CurT = ForwardKine(CurQ, Alpha, DH)
CurR = [CurT(1,1), CurT(1,2), CurT(1,3);
          CurT(2,1), CurT(2,2), CurT(2,3);
          CurT(3,1), CurT(3,2), CurT(3,3)];
Cureuler = rotm2eul(CurR, "ZYX");
CurP = [CurT(1,4),CurT(2,4),CurT(3,4)];

%更新关系矩阵
z = Cureuler(1);
y = Cureuler(2);
x = Cureuler(3);
Te = [cos(y)*cos(z),-sin(z), 0;
      cos(y)*sin(z), cos(z), 0;
      -sin(y), 0,     1];
end

for i=1:6
    if abs(CurQ(i)) > pi
            CurQ(i) = CurQ(i) - round(CurQ(i) / (2 * pi)) * 2 * pi;
    end
end