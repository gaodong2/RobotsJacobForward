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
a = -27.931 * degree;
b = -20.406 * degree;
c = -91.235 * degree;
d = 62.928 * degree;
e = 88.816 * degree;
f = 73.872 * degree;
CurQ = [a, b, c, d, e, f];

% 连杆角度
aa = 0 * degree;
ba = 0.01 * degree;
ca = 0 * degree;
da = 0.04 * degree;
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
a = -15.422 * degree;
b = -26.804 * degree;
c = -35.530 * degree;
d = 53.589 * degree;
e = 42.961 * degree;
f = 62.884 * degree;
TarQ = [a, b, c, d, e, f]; 
TarT = ForwardKine(TarQ, Alpha, DH);
TarR = [TarT(1,1), TarT(1,2), TarT(1,3);
          TarT(2,1), TarT(2,2), TarT(2,3);
          TarT(3,1), TarT(3,2), TarT(3,3)];
Tareuler = rotm2eul(TarR, "ZYX");
TarP = [TarT(1,4), TarT(2,4), TarT(3,4)];

%% 进入牛顿差值迭代，计算关节角度
for i=1:10
deltaP = [TarP(1)-CurP(1),TarP(2)-CurP(2),TarP(3)-CurP(3)];
deltaEuler = [Tareuler(3)-Cureuler(3);Tareuler(2)-Cureuler(2);Tareuler(1)-Cureuler(1)];

rotational = Te * deltaEuler;
Fq = [deltaP(1), deltaP(2), deltaP(3), rotational(1), rotational(2), rotational(3)]

J = Jacobian_Forward(CurQ, Alpha, DH);

Jcross = pinv(J'*J)*J';
deltaQ = Jcross * Fq';
CurQ = CurQ + deltaQ';

%更新当前位置
CurT = ForwardKine(CurQ, Alpha, DH);
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

