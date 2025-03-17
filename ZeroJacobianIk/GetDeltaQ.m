function  deltaQ = GetDeltaQ(Fq, jacob)

% calculate W_n
constant_Wn = eye(6) * 1;
Eq = [Fq(1) 0 0 0 0 0
    0 Fq(2) 0 0 0 0
    0 0 Fq(3) 0 0 0
    0 0 0 Fq(4) 0 0
    0 0 0 0 Fq(5) 0
    0 0 0 0 0 Fq(6)];
Wn = Eq + constant_Wn;

% calculate H_k
constant_We = eye(6) * 1;
Hk = jacob' * constant_We * jacob + Wn;

% calculate Gk
Gk = jacob' * constant_We * Fq';

%calculate deltaQ
deltaQ = inv(Hk) * Gk;
end

