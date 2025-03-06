function T = ForwardKine(Q, Alpha, DH)

% theta
Rz1=[1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1];

Rz2=[cos(Q(1)+pi/2) -sin(Q(1)+pi/2) 0 0;
     sin(Q(1)+pi/2)  cos(Q(1)+pi/2) 0 0;
     0               0              1 0;
     0               0              0 1];

Rz3=[cos(Q(2)-pi/2) -sin(Q(2)-pi/2) 0 0;
     sin(Q(2)-pi/2)  cos(Q(2)-pi/2) 0 0;
     0               0              1 0;
     0               0              0 1];

Rz4=[cos(Q(3)) -sin(Q(3)) 0 0;
     sin(Q(3))  cos(Q(3)) 0 0;
     0          0         1 0;
     0          0         0 1];

Rz5=[cos(Q(4)-pi/2) -sin(Q(4)-pi/2) 0 0;
     sin(Q(4)-pi/2)  cos(Q(4)-pi/2) 0 0;
     0               0              1 0;
     0               0              0 1];

Rz6=[cos(Q(5)) -sin(Q(5)) 0 0;
     sin(Q(5))  cos(Q(5)) 0 0;
     0          0         1 0;
     0          0         0 1];

Rz7=[cos(Q(6)) -sin(Q(6)) 0 0;
     sin(Q(6))  cos(Q(6)) 0 0;
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

Tx3 = [1 0 0 DH(2);
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];

Tx4 = [1 0 0 DH(3);
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
       0 0 1 -DH(1);
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
       0 0 1 DH(4);
       0 0 0 1];

Tz6 = [1 0 0 0;
       0 1 0 0;
       0 0 1 DH(5);
       0 0 0 1];

Tz7 = [1 0 0 0;
       0 1 0 0;
       0 0 1 DH(6);
       0 0 0 1];

%alpha
Rx1=[1  0  0 0;
     0 -1  0 0;
     0  0 -1 0;
     0  0  0 1];

Rx2=[1 0             0            0;
     0 cos(Alpha(1)+pi/2) -sin(Alpha(1)+pi/2) 0;
     0 sin(Alpha(1)+pi/2)  cos(Alpha(1)+pi/2) 0;
     0 0             0            1];

Rx3=[1 0        0       0;
     0 cos(Alpha(2)) -sin(Alpha(2)) 0;
     0 sin(Alpha(2))  cos(Alpha(2)) 0;
     0 0        0       1];

Rx4=[1 0           0          0;
     0 cos(Alpha(3)+pi) -sin(Alpha(3)+pi) 0;
     0 sin(Alpha(3)+pi)  cos(Alpha(3)+pi) 0;
     0 0           0          1];

Rx5=[1 0             0            0;
     0 cos(Alpha(4)-pi/2) -sin(Alpha(4)-pi/2) 0;
     0 sin(Alpha(4)-pi/2)  cos(Alpha(4)-pi/2) 0;
     0 0             0            1];

Rx6=[1 0  0 0;
     0 cos(Alpha(5)+pi/2) -sin(Alpha(5)+pi/2) 0;
     0 sin(Alpha(5)+pi/2)  cos(Alpha(5)+pi/2) 0;
     0 0  0 1];

Rx7=[1 0        0       0;
     0 cos(Alpha(6)) -sin(Alpha(6)) 0;
     0 sin(Alpha(6))  cos(Alpha(6)) 0;
     0 0        0       1];

%正运动学
T00 = Rz1*Tx1*Tz1*Rx1;
T01 = T00*Rz2*Tx2*Tz2*Rx2;
T02 = T01*Rz3*Tx3*Tz3*Rx3;
T03 = T02*Rz4*Tx4*Tz4*Rx4;
T04 = T03*Rz5*Tx5*Tz5*Rx5;
T05 = T04*Rz6*Tx6*Tz6*Rx6;
T = T05*Rz7*Tx7*Tz7*Rx7;

end

