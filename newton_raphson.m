clear;clc
%���±�ţ���ԭ���еĽڵ�1,2,3,4,5�������α��Ϊ5,1,2,3,4������1-4��ΪPQ�ڵ㣬5��Ϊƽ��ڵ�
y=zeros(5,5);%��֧·����ԭʼ����
%����ԭʼ����
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
Y=0; %��ڵ㵼�ɾ���
%�󻥵���
for i=1:5
    for j=1:5
        if i~=j
            Y(i,j)=-y(i,j);
        end
    end
end
%���Ե���
for i=1:5
    Y(i,i)=sum(y(i,:));
end
Y
G=real(Y);
B=imag(Y);
%ԭʼ�ڵ�ע�빦��
S = [0.2+0.2i, -0.45-0.15i, -0.4-0.05i, -0.6-0.1i, 0] %ƽ��ڵ㹦����Ϊ0,����
P=real(S);
Q=imag(S);
%����ֵ
U=ones(1,5);U(5)=1.06;%�ڵ��ѹ��ֵ ��ƽ��ڵ��ѹ����1.06
e=zeros(1,5); %�ڵ��ѹ�Ƕ�
ox=ones(8,1);%�������̱���delata x (delata_theta; delta_u')
fx=ones(8,1);%�������̺���fx
count=0; %�����������
while max(fx)>1e-5
    for i=1:4 
        for j=1:4
            H(i,j)=0;N(i,j)=0;M(i,j)=0;L(i,j)=0;oP(i)=0;oQ(i)=0;%
        end
    end
    for i=1:4
        for j=1:5 %д���������̵���ʽ,��deltaP deltaQ
            oP(i)=oP(i)-U(i)*U(j)*(G(i,j)*cos(e(i)-e(j))+B(i,j)*sin(e(i)-e(j)));
            oQ(i)=oQ(i)-U(i)*U(j)*(G(i,j)*sin(e(i)-e(j))-B(i,j)*cos(e(i)-e(j)));
        end
        oP(i)=oP(i)+P(i); oQ(i)=oQ(i)+Q(i);
    end
    fx=[oP,oQ]'; %jacobi����Ⱥ����
    %���ſ˱Ⱦ���
    %��i~=jʱ����H,N,M,L ���£�
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
    %��i=j ʱH,N,M,L���£�
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
    J=[H,N;M,L]; %J Ϊ�ſ˱Ⱦ���
    ox=-((inv(J))*fx); %����������,�õ�������delta_theta,delta_U'
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
U,e,count %���ڵ��ѹ��ֵ����λ�ǣ��ܵ�������
%��ڵ�ע��ľ�����
i=5;
for j=1:5 %����ƽ��ڵ�ķ�����������й����޹�����
    P(i)=U(i)*U(j)*(G(i,j)*cos(e(i)-e(j))+B(i,j)*sin(e(i)-e(j)))+P(i);
    Q(i)=U(i)*U(j)*(G(i,j)*sin(e(i)-e(j))-B(i,j)*cos(e(i)-e(j)))+Q(i);
end
S(5)=P(5)+Q(5)*sqrt(-1);
S %���ڵ�ע�빦��
I=Y*U'%����ڵ�ע�����
 
 for i=1:5 %��֧·�����Ĺ���
    for j=1:5 
        L_S(i,j)=conj(U(i))*(conj(U(i))-conj(U(j)))*conj(y(i,j));
    end 
 end 
 L_S

for i=1:5 %��֧·����
    for j=1:5 
        L_I(i,j)=conj(L_S(i,j))/conj(U(i));
    end 
end
L_I
%%%%%%%%%%%%%������δ���֧·��Ĳ�ȷ���Ƿ������ȷ%%%%%%%%%%%%%%%%%
S_loss = zeros(5,5); %
for i=1:4 %��֧·���
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
