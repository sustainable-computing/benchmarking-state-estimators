function [V,k]=WLS_SE(Y,z,zType,R,iter_max,threshold,Vtrue)
% SE
n=length(Y);%number of nodes
G=real(Y);%real part of Y
B=imag(Y);%imaginary part of Y
%initial guess
% x=[.01*randn(1,n-1),ones(1,n)+.03*randn(1,n)];
% x=[zeros(1,n-1),ones(1,n)];
x=[[-2*pi/3,-4*pi/3,repmat([0,-2*pi/3,-4*pi/3],1,n/3-1)],ones(1,n)];
% x=[angle(Vtrue(2:end)),abs(Vtrue)]; %for testing the accuracy with correct solution as starting point
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
        elseif zType(m,1)==7 %|Iij|
            i=zType(m,2);
            j=zType(m,3);
            h(m)=sqrt((G(i,j)^2+B(i,j)^2)*(V(i)^2+V(j)^2-2*V(i)*V(j)*cos(Th(i)-Th(j))));
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
        elseif zType(m,1)==7 %|Iij|
            i=zType(m,2);
            j=zType(m,3);
            if i>1
                H(m,i-1)=(G(i,j)^2+B(i,j)^2)*V(i)*V(j)*sin(Th(i)-Th(j))/h(m);
            end
            if j>1
                H(m,j-1)=-(G(i,j)^2+B(i,j)^2)*V(i)*V(j)*sin(Th(i)-Th(j))/h(m);
            end
            H(m,i+n-1)=(G(i,j)^2+B(i,j)^2)*(V(i)-V(j)*cos(Th(i)-Th(j)))/h(m);
            H(m,i+n-1)=(G(i,j)^2+B(i,j)^2)*(V(j)-V(i)*cos(Th(i)-Th(j)))/h(m);
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