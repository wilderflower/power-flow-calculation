clear
%������һ����������ڵ���ڶ���,��ƽ��ڵ���Ҫ�����
NN=input('����������ڵ���NN:')
r=input('����������PV�ڵ���r:');
N1=NN-1;          % N1����Ҫ����Ľڵ���

A = input('�������������:');
y = input('������֧·���ɾ���:');
Y=A*y*A';
'���ɾ������'
Y

G=real(Y);
B=imag(Y);

%��ֵ����
for chuzhi=1:N1
    delt(chuzhi)=0;
end
 
for chuzhi=1:N1-r
    u(chuzhi)=1.0;
end 
   
 
%����PQ�ڵ㹦������
'������ڵ�����빦��'
for jiediangl=1:N1
    sp(jiediangl)=input('�����������ڵ�����빦��:');
end
p=real(sp);
q=imag(sp);
for jiediangl=1:r
    u(jiediangl+N1-r)=input('�����������ڵ�ĵ�ѹ:');
end

k=0;precision=1;

geidingdy=input('������ƽ��ڵ������ѹ:');

while precision>0.00001;
    delt(NN)=0;u(NN)=geidingdy;
   for m=1:N1
       for n=1:N1+1
           pt(n)=u(m)*u(n)*(G(m,n)*cos(delt(m)-delt(n))+B(m,n)*sin(delt(m)-delt(n)));
       end
       pp(m)=p(m)-sum(pt);
       p
       sum(pt)
   end
 
   for m=1:N1-r
       for n=1:N1+1
           qt(n)=u(m)*u(n)*(G(m,n)*sin(delt(m)-delt(n))-B(m,n)*cos(delt(m)-delt(n)));
       end
       qq(m)=q(m)-sum(qt);
   end

   for m=1:N1
       for n=1:N1+1
           h0(n)=u(m)*u(n)*(G(m,n)*sin(delt(m)-delt(n))-B(m,n)*cos(delt(m)-delt(n)));
           n0(n)=-u(m)*u(n)*(G(m,n)*cos(delt(m)-delt(n))+B(m,n)*sin(delt(m)-delt(n)));
           j0(n)=-u(m)*u(n)*(G(m,n)*cos(delt(m)-delt(n))+B(m,n)*sin(delt(m)-delt(n)));
           L0(n)=-u(m)*u(n)*(G(m,n)*sin(delt(m)-delt(n))-B(m,n)*cos(delt(m)-delt(n)));
       end
       H(m,m)=sum(h0)-u(m)^2*(G(m,m)*sin(delt(m)-delt(m))-B(m,m)*cos(delt(m)-delt(m)));
       N(m,m)=sum(n0)-2*u(m)^2*G(m,m)+u(m)^2*(G(m,m)*cos(delt(m)-delt(m))+B(m,m)*sin(delt(m)-delt(m)));
       J(m,m)=sum(j0)+u(m)^2*(G(m,m)*cos(delt(m)-delt(m))+B(m,m)*sin(delt(m)-delt(m)));
       L(m,m)=sum(L0)+2*u(m)^2*B(m,m)+u(m)^2*(G(m,m)*sin(delt(m)-delt(m))-B(m,m)*cos(delt(m)-delt(m)));
   end
 
   for m=1:N1
       for n=1:N1
           if m==n
           else
               H(m,n)=-u(m)*u(n)*(G(m,n)*sin(delt(m)-delt(n))-B(m,n)*cos(delt(m)-delt(n)));
               J(m,n)=u(m)*u(n)*(G(m,n)*cos(delt(m)-delt(n))+B(m,n)*sin(delt(m)-delt(n)));
               N(m,n)=-J(m,n);
               L(m,n)=H(m,n);
           end
       end
   end
   for m=1:N1
       PP(m)=pp(m);
   end
   for m=1:N1-r
       PP(m+N1)=qq(m);
   end
   
   JJ=[H N;J L];
   Jt=JJ(:,1:2*N1-r);
   Jn=Jt(1:2*N1-r,:);
   uu=-inv(Jn)*PP';
   precision=max(abs(uu));
   for n=1:N1
       delt(n)=delt(n)+uu(n);
   end
   for n=1:N1-r
       u(n)=u(n)+uu(n+N1);
   end

   k=k+1;
end
%������Ǽ��㹦��
U=u.*exp(i*delt)
S=conj(conj(diag(U))*Y*conj(U'))
for m=1:N1+1
    for n=1:N1+1
        SS(m,n)=U(m)*(conj(U(m))-conj(U(n)))*conj(-Y(m,n));
    end
end

'��������'
k-1
'�ڵ��ѹ'
u
'�ڵ�Ƕ�'   
delt
'���빦��'
S
'��·����'
SS