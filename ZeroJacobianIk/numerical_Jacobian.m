%清空缓存
clc;
clear all;
%程序说明：这是6公斤弧焊机器人的标定求解雅克比矩阵的程序，带有减速比和耦合系统的求解
%参数说明：
%a~f分别代表了6个关节角
%aa,ab,ba,bb,ca,da,ea分别代表了各杆杆长
%ak、bk、ck dk ek fk 代表6个关节的减速比，jk代表5轴和6轴的耦合系数
%fa,fb,fc代表了工具坐标系的偏差
%经测试发现a,和aa还有fa都没有实际的意义，加入fa会使雅克比矩阵不满秩，而aa和a本身没有意义；

syms a b c d e f;%关节角
syms aa ba ca da ea fa;%减速比
syms d2 a3 a4 d5 d6 d7;%杆长

% theta
Rz1=[1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1];

Rz2=[cos(a+pi/2) -sin(a+pi/2) 0 0;
     sin(a+pi/2)  cos(a+pi/2) 0 0;
     0               0              1 0;
     0               0              0 1];

Rz3=[cos(b-pi/2) -sin(b-pi/2) 0 0;
     sin(b-pi/2)  cos(b-pi/2) 0 0;
     0               0              1 0;
     0               0              0 1];

Rz4=[cos(c) -sin(c) 0 0;
     sin(c)  cos(c) 0 0;
     0          0         1 0;
     0          0         0 1];

Rz5=[cos(d-pi/2) -sin(d-pi/2) 0 0;
     sin(d-pi/2)  cos(d-pi/2) 0 0;
     0               0              1 0;
     0               0              0 1];

Rz6=[cos(e) -sin(e) 0 0;
     sin(e)  cos(e) 0 0;
     0          0         1 0;
     0          0         0 1];

Rz7=[cos(f) -sin(f) 0 0;
     sin(f)  cos(f) 0 0;
     0          0         1 0;
     0          0         0 1];

%a
Tx1 = [1 0 0 0;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];

Tx2 = [1 0 0 0;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];

Tx3 = [1 0 0 a3;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];

Tx4 = [1 0 0 a4;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];

Tx5 = [1 0 0 0;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];

Tx6 = [1 0 0 0;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];

Tx7 = [1 0 0 0;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];

%d
Tz1 = [1 0 0 0;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];

Tz2 = [1 0 0 0;
       0 1 0 0;
       0 0 1 -d2;
       0 0 0 1];

Tz3 = [1 0 0 0;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];

Tz4 = [1 0 0 0;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];

Tz5 = [1 0 0 0;
       0 1 0 0;
       0 0 1 d5;
       0 0 0 1];

Tz6 = [1 0 0 0;
       0 1 0 0;
       0 0 1 d6;
       0 0 0 1];

Tz7 = [1 0 0 0;
       0 1 0 0;
       0 0 1 d7;
       0 0 0 1];

%alpha
Rx1=[1  0  0 0;
     0 -1  0 0;
     0  0 -1 0;
     0  0  0 1];

Rx2=[1 0             0            0;
     0 cos(aa+pi/2) -sin(aa+pi/2) 0;
     0 sin(aa+pi/2)  cos(aa+pi/2) 0;
     0 0             0            1];

Rx3=[1 0        0       0;
     0 cos(ba) -sin(ba) 0;
     0 sin(ba)  cos(ba) 0;
     0 0        0       1];

Rx4=[1 0           0          0;
     0 cos(ca+pi) -sin(ca+pi) 0;
     0 sin(ca+pi)  cos(ca+pi) 0;
     0 0           0          1];

Rx5=[1 0             0            0;
     0 cos(da-pi/2) -sin(da-pi/2) 0;
     0 sin(da-pi/2)  cos(da-pi/2) 0;
     0 0             0            1];

Rx6=[1 0  0 0;
     0 cos(ea+pi/2) -sin(ea+pi/2) 0;
     0 sin(ea+pi/2)  cos(ea+pi/2) 0;
     0 0  0 1];

Rx7=[1 0        0       0;
     0 cos(fa) -sin(fa) 0;
     0 sin(fa)  cos(fa) 0;
     0 0        0       1];

%正运动学
T00 = Rz1*Tx1*Tz1*Rx1;
T01 = T00*Rz2*Tx2*Tz2*Rx2;
T02 = T01*Rz3*Tx3*Tz3*Rx3;
T03 = T02*Rz4*Tx4*Tz4*Rx4;
T04 = T03*Rz5*Tx5*Tz5*Rx5;
T05 = T04*Rz6*Tx6*Tz6*Rx6;
T06 = T05*Rz7*Tx7*Tz7*Rx7;
G=Rz1*Tx1*Tz1*Rx1*Rz2*Tx2*Tz2*Rx2*Rz3*Tx3*Tz3*Rx3*Rz4*Tx4*Tz4*Rx4*Rz5*Tx5*Tz5*Rx5*Rz6*Tx6*Tz6*Rx6*Rz7*Tx7*Tz7*Rx7;

%只关心x,y,z，因为不考虑姿态，所以没有加入姿态的雅克比矩阵
%如需要求解姿态的求解，姿态求解不能简单的求解，还需要对正运动学矩阵做一次求Ra,Rb,Rc的代表式

x=G(1,4);
y=G(2,4);
z=G(3,4);

%对角速度求解时加入
% A = G(3,2)/G(3,3);
% B = -G(3,1);
% C = G(2,1) /G(1,1);

%求雅克比矩阵
%J=jacobian([x;y;z],[a b c d e f g d2 a3 a4 d5 d6 d7 fa fb fc ak bk ck dk ek fk gk])
J=jacobian([x;y;z],[a b c d e f]);
%J=jacobian([x;y;z],[ a b c d e f]);
% S = svd(J)

%化简：用simplify(simple)对函数进行自动化简！太棒了
JA = [J(1,1), J(1,2), J(1,3), J(1,4), J(1,5), J(1,6);
      J(2,1), J(2,2), J(2,3), J(2,4), J(2,5), J(2,6);
      J(3,1), J(3,2), J(3,3), J(3,4), J(3,5), J(3,6);
    T00(1,3), T01(1,3), T02(1,3), T03(1,3), T04(1,3), T05(1,3);
    T00(2,3), T01(2,3), T02(2,3), T03(2,3), T04(2,3), T05(2,3);
    T00(3,3), T01(3,3), T02(3,3), T03(3,3), T04(3,3), T05(3,3)];

%化简的另一种公式
% B = simple(J)
[m,n] = size(JA);

for i= 1:m
    for j =1:n
        JA(i,j)
    end
end
%求雅克比矩阵的秩
 %det(A)


