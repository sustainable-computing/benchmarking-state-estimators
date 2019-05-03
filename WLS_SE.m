function [V,k]=WLS_SE(Y,z,zType,R,iter_max,threshold,Vtrue)
% SE
n=length(Y);%number of nodes
G=real(Y);%real part of Y
B=imag(Y);%imaginary part of Y
%initial guess
% x=[.01*randn(1,n-1),ones(1,n)+.03*randn(1,n)];
% x=[zeros(1,n-1),ones(1,n)];
% x=[[-2*pi/3,-4*pi/3,repmat([0,-2*pi/3,-4*pi/3],1,n/3-1)],ones(1,n)];
x=[angle(Vtrue(2:end)),abs(Vtrue)]; %for testing the accuracy with correct solution as starting point
for k=1:iter_max
    V=x(n:end);%voltage magnitudes
    Th=[0,x(1:n-1)];%voltage angles (Theta(1)=0 reference)
    %calculating measurement function given the states h(x)
    h=zeros(length(z),1);
    for m=1:length(z)
        if zType(m,1)==1 %Pij
            i=zType(m,2);
            j=zType(m,3);
            h(m)=V(i)^2*(-G(i,j))-V(i)*V(j)*(-G(i,j)*cos(Th(i)-Th(j))-B(i,j)*sin(Th(i)-Th(j)));
        elseif zType(m,1)==2 %Pi
            i=zType(m,2);
            h(m)=0;
            for j=1:length(Y)
                h(m)=h(m)+V(i)*V(j)*(G(i,j)*cos(Th(i)-Th(j))+B(i,j)*sin(Th(i)-Th(j)));
            end
        elseif zType(m,1)==3 %Qij
            i=zType(m,2);
            j=zType(m,3);
            h(m)=-V(i)^2*(-B(i,j))-V(i)*V(j)*(-G(i,j)*sin(Th(i)-Th(j))+B(i,j)*cos(Th(i)-Th(j)));
        elseif zType(m,1)==4 %Qi
            i=zType(m,2);
            h(m)=0;
            for j=1:length(Y)
                h(m)=h(m)+V(i)*V(j)*(G(i,j)*sin(Th(i)-Th(j))-B(i,j)*cos(Th(i)-Th(j)));
            end
        elseif zType(m,1)==5 %|Vi|
            i=zType(m,2);
            h(m)=V(i);
        elseif zType(m,1)==6 %Theta Vi
            i=zType(m,2);
            h(m)=Th(i);
        elseif zType(m,1)==7 %real Iij
            i=zType(m,2);%sending end
            j=zType(m,3);%receiving end
            ph=zType(m,4);%phase
            a1=3*(i-1)+1;
            b1=3*(i-1)+2;
            c1=3*(i-1)+3;
            a2=3*(j-1)+1;
            b2=3*(j-1)+2;
            c2=3*(j-1)+3;
            y=-Y([a1,b1,c1],[a2,b2,c2]);
            g=real(y);
            b=imag(y);
            h(m)=g(ph,1)*(V(a1)*cos(Th(a1))-V(a2)*cos(Th(a2)))-b(ph,1)*(V(a1)*sin(Th(a1))-V(a2)*sin(Th(a2)))+...
                 g(ph,2)*(V(b1)*cos(Th(b1))-V(b2)*cos(Th(b2)))-b(ph,2)*(V(b1)*sin(Th(b1))-V(b2)*sin(Th(b2)))+...
                 g(ph,3)*(V(c1)*cos(Th(c1))-V(c2)*cos(Th(c2)))-b(ph,3)*(V(c1)*sin(Th(c1))-V(c2)*sin(Th(c2))); 
        elseif zType(m,1)==8 %imag Iij
            i=zType(m,2);%sending end
            j=zType(m,3);%receiving end
            ph=zType(m,4);%phase
            a1=3*(i-1)+1;
            b1=3*(i-1)+2;
            c1=3*(i-1)+3;
            a2=3*(j-1)+1;
            b2=3*(j-1)+2;
            c2=3*(j-1)+3;
            y=-Y([a1,b1,c1],[a2,b2,c2]);
            g=real(y);
            b=imag(y);
            h(m)=g(ph,1)*(V(a1)*sin(Th(a1))-V(a2)*sin(Th(a2)))+b(ph,1)*(V(a1)*cos(Th(a1))-V(a2)*cos(Th(a2)))+...
                 g(ph,2)*(V(b1)*sin(Th(b1))-V(b2)*sin(Th(b2)))+b(ph,2)*(V(b1)*cos(Th(b1))-V(b2)*cos(Th(b2)))+...
                 g(ph,3)*(V(c1)*sin(Th(c1))-V(c2)*sin(Th(c2)))+b(ph,3)*(V(c1)*cos(Th(c1))-V(c2)*cos(Th(c2)));
        end
    end

    %building the jacobian H
    H=zeros(length(z),length(x));

    for m=1:length(z)
        if zType(m,1)==1 %Pij
            i=zType(m,2);
            j=zType(m,3);
            if i~=1
                H(m,i-1)=V(i)*V(j)*(-G(i,j)*sin(Th(i)-Th(j))+B(i,j)*cos(Th(i)-Th(j)));
            end
            if j~=1
                H(m,j-1)=-V(i)*V(j)*(-G(i,j)*sin(Th(i)-Th(j))+B(i,j)*cos(Th(i)-Th(j)));
            end
            H(m,i+n-1)=-V(j)*(-G(i,j)*cos(Th(i)-Th(j))-B(i,j)*sin(Th(i)-Th(j)))+2*(-G(i,j))*V(i);
            H(m,j+n-1)=-V(i)*(-G(i,j)*cos(Th(i)-Th(j))-B(i,j)*sin(Th(i)-Th(j)));
        elseif zType(m,1)==2 %Pi
            i=zType(m,2);
            for j=1:n
                if i~=j
                    if j>1
                        H(m,j-1)=V(i)*V(j)*(G(i,j)*sin(Th(i)-Th(j))-B(i,j)*cos(Th(i)-Th(j)));
                    end
                    H(m,j+n-1)=V(i)*(G(i,j)*cos(Th(i)-Th(j))+B(i,j)*sin(Th(i)-Th(j)));
                end
            end
            if i>1
                H(m,i-1)=-V(i)^2*B(i,i);
                for j=1:n
                    H(m,i-1)=H(m,i-1)+V(i)*V(j)*(-G(i,j)*sin(Th(i)-Th(j))+B(i,j)*cos(Th(i)-Th(j)));
                end
            end
            H(m,i+n-1)=+V(i)*G(i,i);
            for j=1:n
                H(m,i+n-1)=H(m,i+n-1)+V(j)*(G(i,j)*cos(Th(i)-Th(j))+B(i,j)*sin(Th(i)-Th(j)));
            end
        elseif zType(m,1)==3 %Qij
            i=zType(m,2);
            j=zType(m,3);
            if i~=1
                H(m,i-1)=-V(i)*V(j)*(-G(i,j)*cos(Th(i)-Th(j))-B(i,j)*sin(Th(i)-Th(j)));
            end
            if j~=1
                H(m,j-1)=V(i)*V(j)*(-G(i,j)*cos(Th(i)-Th(j))-B(i,j)*sin(Th(i)-Th(j)));
            end
            H(m,i+n-1)=-V(j)*(-G(i,j)*sin(Th(i)-Th(j))+B(i,j)*cos(Th(i)-Th(j)))-2*(-B(i,j))*V(i);
            H(m,j+n-1)=-V(i)*(-G(i,j)*sin(Th(i)-Th(j))+B(i,j)*cos(Th(i)-Th(j)));
        elseif zType(m,1)==4 %Qi
            i=zType(m,2);
            for j=1:n
                if i~=j
                    if j>1
                        H(m,j-1)=V(i)*V(j)*(-G(i,j)*cos(Th(i)-Th(j))-B(i,j)*sin(Th(i)-Th(j)));
                    end
                    H(m,j+n-1)=V(i)*(G(i,j)*sin(Th(i)-Th(j))-B(i,j)*cos(Th(i)-Th(j)));
                end
            end
            if i>1
                H(m,i-1)=-V(i)^2*G(i,i);
                for j=1:n
                    H(m,i-1)=H(m,i-1)+V(i)*V(j)*(G(i,j)*cos(Th(i)-Th(j))+B(i,j)*sin(Th(i)-Th(j)));
                end
            end
            H(m,i+n-1)=-V(i)*B(i,i);
            for j=1:n
                H(m,i+n-1)=H(m,i+n-1)+V(j)*(G(i,j)*sin(Th(i)-Th(j))-B(i,j)*cos(Th(i)-Th(j)));
            end
        elseif zType(m,1)==5 %|Vi|
            i=zType(m,2);
            H(m,i+n-1)=1;
        elseif zType(m,1)==6 %Theta Vi
            i=zType(m,2);
            H(m,i-1)=1;
        elseif zType(m,1)==7 %real Iij
            i=zType(m,2);
            j=zType(m,3);
            ph=zType(m,4);%phase
            a1=3*(i-1)+1;
            b1=3*(i-1)+2;
            c1=3*(i-1)+3;
            a2=3*(j-1)+1;
            b2=3*(j-1)+2;
            c2=3*(j-1)+3;
            y=-Y([a1,b1,c1],[a2,b2,c2]);
            g=real(y);
            b=imag(y);
            %angles
            if a1>1
                H(m,a1-1)=-g(ph,1)*V(a1)*sin(Th(a1))-b(ph,1)*V(a1)*cos(Th(a1));
            end
            H(m,b1-1)=-g(ph,2)*V(b1)*sin(Th(b1))-b(ph,2)*V(b1)*cos(Th(b1));
            H(m,c1-1)=-g(ph,3)*V(c1)*sin(Th(c1))-b(ph,3)*V(c1)*cos(Th(c1));
            H(m,a2-1)= g(ph,1)*V(a2)*sin(Th(a2))+b(ph,1)*V(a2)*cos(Th(a2));
            H(m,b2-1)= g(ph,2)*V(b2)*sin(Th(b2))+b(ph,2)*V(b2)*cos(Th(b2));
            H(m,c2-1)= g(ph,3)*V(c2)*sin(Th(c2))+b(ph,3)*V(c2)*cos(Th(c2));
            %magnitudes
            H(m,a1+n-1)= g(ph,1)*cos(Th(a1))-b(ph,1)*sin(Th(a1));
            H(m,b1+n-1)= g(ph,2)*cos(Th(b1))-b(ph,2)*sin(Th(b1));
            H(m,c1+n-1)= g(ph,3)*cos(Th(c1))-b(ph,3)*sin(Th(c1));
            H(m,a2+n-1)=-g(ph,1)*cos(Th(a2))+b(ph,1)*sin(Th(a2));
            H(m,b2+n-1)=-g(ph,2)*cos(Th(b2))+b(ph,2)*sin(Th(b2));
            H(m,c2+n-1)=-g(ph,3)*cos(Th(c2))+b(ph,3)*sin(Th(c2));
            
        elseif zType(m,1)==8 %imag Iij
            i=zType(m,2);
            j=zType(m,3);
            ph=zType(m,4);%phase
            a1=3*(i-1)+1;
            b1=3*(i-1)+2;
            c1=3*(i-1)+3;
            a2=3*(j-1)+1;
            b2=3*(j-1)+2;
            c2=3*(j-1)+3;
            y=-Y([a1,b1,c1],[a2,b2,c2]);
            g=real(y);
            b=imag(y);
            %angles
            if a1>1
                H(m,a1-1)= g(ph,1)*V(a1)*cos(Th(a1))-b(ph,1)*V(a1)*sin(Th(a1));
            end
            H(m,b1-1)= g(ph,2)*V(b1)*cos(Th(b1))-b(ph,2)*V(b1)*sin(Th(b1));
            H(m,c1-1)= g(ph,3)*V(c1)*cos(Th(c1))-b(ph,3)*V(c1)*sin(Th(c1));
            H(m,a2-1)=-g(ph,1)*V(a2)*cos(Th(a2))+b(ph,1)*V(a2)*sin(Th(a2));
            H(m,b2-1)=-g(ph,2)*V(b2)*cos(Th(b2))+b(ph,2)*V(b2)*sin(Th(b2));
            H(m,c2-1)=-g(ph,3)*V(c2)*cos(Th(c2))+b(ph,3)*V(c2)*sin(Th(c2));
            %magnitudes
            H(m,a1+n-1)= g(ph,1)*sin(Th(a1))+b(ph,1)*cos(Th(a1));
            H(m,b1+n-1)= g(ph,2)*sin(Th(b1))+b(ph,2)*cos(Th(b1));
            H(m,c1+n-1)= g(ph,3)*sin(Th(c1))+b(ph,3)*cos(Th(c1));
            H(m,a2+n-1)=-g(ph,1)*sin(Th(a2))-b(ph,1)*cos(Th(a2));
            H(m,b2+n-1)=-g(ph,2)*sin(Th(b2))-b(ph,2)*cos(Th(b2));
            H(m,c2+n-1)=-g(ph,3)*sin(Th(c2))-b(ph,3)*cos(Th(c2));
        end
    end

    %The right hand side of the equation
    RHS=H.'*R^-1*(z-h);

    %The gain matrix
    GM=H.'*R^-1*H;

    %Cholesky decomposition
    % L = chol(G,'lower');

    delta_x=GM\RHS;
    if abs(delta_x)<threshold
        break
    end
    x=x+(delta_x).';
end
Vmag=x(length(Y):end);%voltage magnitudes
Theta=[0,x(1:length(Y)-1)];%voltage angles (Theta(1)=0 reference)
V=Vmag.*(cos(Theta)+1i*sin(Theta));
end
