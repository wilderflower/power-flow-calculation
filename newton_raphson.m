clear;clc
%重新编号，把原题中的节点1,2,3,4,5重新依次编号为5,1,2,3,4，其中1-4号为PQ节点，5号为平衡节点
y=zeros(5,5);%各支路导纳原始数据
%输入原始数据
z = [0.0000 + 0.0000i   0.06+0.18i         0.06+0.18i         0.04+0.12i         0.02+0.06i
     0.0000 + 0.0000i   0.0000 + 0.0000i   0.01+0.03i         0.0000 + 0.0000i   0.08+0.24i
     0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.08+0.24i         0.0000 + 0.0000i
     0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i];
for i=1:5
    for j=i+1:5
        if z(i,j) ~= 0
            y(i,j) = 1/z(i,j);
            y(j,i)=y(i,j);
        end
    end
end
Y=0; %求节点导纳矩阵
%求互导纳
for i=1:5
    for j=1:5
        if i~=j
            Y(i,j)=-y(i,j);
        end
    end
end
%求自导纳
for i=1:5
    Y(i,i)=sum(y(i,:));
end
Y
G=real(Y);
B=imag(Y);
%原始节点注入功率
S = [0.2+0.2i, -0.45-0.15i, -0.4-0.05i, -0.6-0.1i, 0] %平衡节点功率设为0,待求
P=real(S);
Q=imag(S);
%赋初值
U=ones(1,5);U(5)=1.06;%节点电压幅值 设平衡节点电压保持1.06
e=zeros(1,5); %节点电压角度
ox=ones(8,1);%修正方程变量delata x (delata_theta; delta_u')
fx=ones(8,1);%修正方程函数fx
count=0; %计算迭代次数
while max(fx)>1e-5
    for i=1:4 
        for j=1:4
            H(i,j)=0;N(i,j)=0;M(i,j)=0;L(i,j)=0;oP(i)=0;oQ(i)=0;%
        end
    end
    for i=1:4
        for j=1:5 %写成修正方程的形式,求deltaP deltaQ
            oP(i)=oP(i)-U(i)*U(j)*(G(i,j)*cos(e(i)-e(j))+B(i,j)*sin(e(i)-e(j)));
            oQ(i)=oQ(i)-U(i)*U(j)*(G(i,j)*sin(e(i)-e(j))-B(i,j)*cos(e(i)-e(j)));
        end
        oP(i)=oP(i)+P(i); oQ(i)=oQ(i)+Q(i);
    end
    fx=[oP,oQ]'; %jacobi矩阵等号左边
    %求雅克比矩阵
    %当i~=j时候求H,N,M,L 如下：
    for i=1:4
        for j=1:4
            if i~=j 
                H(i,j)=-U(i)*U(j)*(G(i,j)*sin(e(i)-e(j))-B(i,j)*cos(e(i)-e(j)));
                N(i,j)=-U(i)*U(j)*(G(i,j)*cos(e(i)-e(j))+B(i,j)*sin(e(i)-e(j)));
                L(i,j)=H(i,j);
                M(i,j)=-N(i,j);
            end
        end
    end
    %H,N,M,L
    %当i=j 时H,N,M,L如下：
    for i=1:4
        for j=1:5
            if i~=j
                H(i,i)=H(i,i)+U(i)*U(j)*(G(i,j)*sin(e(i)-e(j))-B(i, j)*cos (e(i)-e(j)));                     
                N(i,i)=N(i,i)-U(i)*U(j)*(G(i, j)*cos(e(i)-e(j))+B(i,j)*sin(e(i)-e(j)));
                M(i,i)=M(i,i)-U(i)*U(j)*(G(i,j)*cos(e(i)-e(j))+B(i,j)*sin(e(i)-e(j))); 
                L(i,i)=L(i,i)-U(i)*U(j)*(G(i,j)*sin(e(i)-e(j))-B(i,j)*cos(e(i)-e(j)));
            end
        end
        N(i,i)=N(i,i)-2*(U(i))^2*G(i,i);
        L(i,i)=L(i,i)+2*(U(i))^2*B(i,i);
    end
    J=[H,N;M,L]; %J 为雅克比矩阵
    ox=-((inv(J))*fx); %解修正方程,得到修正量delta_theta,delta_U'
    for i=1:4
        oe(i)=ox(i); 
        oU(i)=ox(i+4)*U(i);
    end
    for i=1:4
        e(i)=e(i)+oe(i);
        U(i)=U(i)+oU(i);
    end
    count=count+1;
end

%ox,
U,e,count %各节点电压幅值，相位角，总迭代次数
%求节点注入的净功率
i=5;
for j=1:5 %计算平衡节点的发电机发出的有功和无功功率
    P(i)=U(i)*U(j)*(G(i,j)*cos(e(i)-e(j))+B(i,j)*sin(e(i)-e(j)))+P(i);
    Q(i)=U(i)*U(j)*(G(i,j)*sin(e(i)-e(j))-B(i,j)*cos(e(i)-e(j)))+Q(i);
end
S(5)=P(5)+Q(5)*sqrt(-1);
S %各节点注入功率
I=Y*U'%求各节点注入电流
 
 for i=1:5 %各支路流过的功率
    for j=1:5 
        L_S(i,j)=conj(U(i))*(conj(U(i))-conj(U(j)))*conj(y(i,j));
    end 
 end 
 L_S

for i=1:5 %各支路电流
    for j=1:5 
        L_I(i,j)=conj(L_S(i,j))/conj(U(i));
    end 
end
L_I
%%%%%%%%%%%%%下面这段代码支路损耗不确定是否计算正确%%%%%%%%%%%%%%%%%
S_loss = zeros(5,5); %
for i=1:4 %各支路损耗
    for j=i+1:5 
        if y(i,j) ~= 0
            S_loss(i,j)=conj(L_I(i,j))*L_I(i,j)*z(i,j); 
        end
    end 
end
for i=1:5
    for j=i:5
        S_loss(j,i)=S_loss(i,j);
    end
end
S_loss
