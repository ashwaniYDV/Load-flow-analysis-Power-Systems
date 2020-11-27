% Load Flow Solution using G-S Method 

clc
clear all
alpha = 1.6; 
%% Data 

% Line Data 
%            From   To    R     X      Y/2
linedata = [ 1      2     0   0.1     0.01
             1      3     0   0.1     0.01
             2      3     0   0.1     0.01];

% Bus Data
%         Bus No   Bus Type    Pg        Qg    Pd       Qd      |V|    delta    Qmin    Qmax      
busdata = [ 1       1          0          0     0        0       1       0         0       0
            2       2          0.6661     0     0        0       1.05    0         0.2     1.2
            3       3          0          0     2.8653   1.2244  1       0         0       0  ];
        
%  Bus Type  Slack =1, PV = 2, PQ = 3

%% Y Bus Formation
R = linedata(:,3); X=linedata(:,4); B=1i*linedata(:,5);
Z = R+1i*X; Y = 1./Z;
nline= length(linedata(:,1));
nbus = max(max(linedata(:,1),linedata(:,2)));

Ybus = zeros(nbus,nbus);
for k=1:nline
    
    % Off Diagaonal Elements
    Ybus(linedata(k,1),linedata(k,2)) = Ybus(linedata(k,1),linedata(k,2)) - Y(k);
    Ybus(linedata(k,2),linedata(k,1)) = Ybus(linedata(k,1),linedata(k,2));
    
    % Diagonal Elements
    Ybus(linedata(k,1),linedata(k,1)) = Ybus(linedata(k,1),linedata(k,1))+ Y(k) + B(k); 
    Ybus(linedata(k,2),linedata(k,2)) = Ybus(linedata(k,2),linedata(k,2))+ Y(k) + B(k);
end

%% Bus Data collection
type = busdata(:,2);
Pg = busdata(:,3); Qg = busdata(:,4); Pd = busdata(:,5); Qd = busdata(:,6); Qmin=busdata(:,9); Qmax = busdata(:,10);
Vmag = busdata(:,7); delta = busdata(:,8); 
V = Vmag.*(cos(delta) + 1i*sin(delta));
P = Pg-Pd; Q = Qg-Qd;  accuracy = 1; 
%% Iteration
iter = 1;
while accuracy >=0.001 && iter < 2 
% Since Bus 1 is a slack bus, We need not calculate its voltage.
for n=2:nbus
    
    if type(n) == 2  % PV Bus
        
        % Finding Q
        % Finding Summation of  Q equation
        A=0;
        for k=1:nbus
            A = A + Ybus(n,k)*V(k);
        end
        Q(n) = -imag(conj(V(n))*A);
        
        % Qlimit check
        if Q(n) >= Qmin(n) && Q(n) <= Qmax(n)
            
             % Finding summation of Voltage equation
            A = 0;
            for k=1:nbus
                if k~=n
                    A = A + Ybus(n,k)*V(k);
                end 
            end
            Vnew(n) = 1/Ybus(n,n) * ((P(n) - 1i*Q(n))/conj(V(n)) - A );
            Vnew_corr(n) = Vmag(n) * Vnew(n)/abs(Vnew(n));
            DV(n) = Vnew_corr(n) - V(n);
            V(n) = Vnew_corr(n);
            
        else                        % PV bus will become PQ bus
            
            if Q(n) < Qmin(n)
                Q(n) = Qmin(n);
            else 
                Q(n) = Qmax(n);
            end
             % To Find summation in the voltage equation
            A = 0;
            for k=1:nbus
                if k~=n
                    A = A + Ybus(n,k)*V(k);
                end 
            end

            Vnew(n) = 1/Ybus(n,n) * ((P(n) - 1i*Q(n))/conj(V(n)) - A ); 
            Vnew_acc(n) = V(n) + alpha*(Vnew(n) - V(n));
            DV(n) = Vnew_acc(n) - V(n);
            V(n) = Vnew_acc(n);
        end   
    else          % PQ Bus
        
        % To Find summation in the voltage equation
        A = 0;
        for k=1:nbus
            if k~=n
                A = A + Ybus(n,k)*V(k);
            end 
        end
        
        Vnew(n) = 1/Ybus(n,n) * ((P(n) - 1i*Q(n))/conj(V(n)) - A );
        Vnew_acc(n) = V(n) + alpha*(Vnew(n) - V(n));
        DV(n) = Vnew_acc(n) - V(n);
        V(n) = Vnew_acc(n);
    end
 
end
accuracy = norm(DV);
iter = iter+1;
end

