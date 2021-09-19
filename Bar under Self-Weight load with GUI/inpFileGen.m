function [nodeData,elementData] = inpFileGen(LT,A1,A2,t,E,rho,P,x0,n,BC)

coordinates = zeros(n+1,2);
L(1:n)  = LT/n;
As      = zeros(1,n+1);
As(1)   = A1;
As(end) = A2;
A       = zeros(1,n);
f       = zeros(1,n+1);

for i = 1 : n
    x = L(i) * i;
    coordinates(i+1,2) = -x;
    
    As(i+1) = A2 + (A1 - A2)*(1-x/LT);
    A(i) = (As(i) + As(i+1))/2;
    
    F = (A(i)*L(i)*t  *  rho)/2*[1,1];
    f(i:i+1) = f(i:i+1) + F;
end

% calculating exact values with integration
syms xx real
Ax = A2 + (A1 - A2)*(1-xx/LT);
Vx = Ax * t;
f_exact = zeros(1,n+1);
for i = 1 :n
    x1 = (i-1) * L(i);
    x2 = i     * L(i);
    N = [1 - (xx-x1)/L(i) ; (xx-x1)/L(i)];
    Fdx = N * Vx * rho ;
    F_exact = double(int(Fdx,x1,x2))';
    f_exact(i:i+1) = f_exact(i:i+1) + F_exact;
end
f;
f_exact;
f = f_exact;

% Generating the Corresponding Input File

fid = fopen('inputFile.txt','w');
fprintf(fid,'%d   %d\n',n+1,n);
for i = 1 : n+1
    if i == 1 
        bc = BC(1);
    elseif i == n+1
        bc = BC(2);
    else
        bc = 1;
    end
    if abs(coordinates(i,2)) == x0
        f(i) = f(i) + P;
    end
    fprintf(fid, '%d   %f  %f  %d  %d  %f  %f \n',i,...
            coordinates(i,1), coordinates(i,2), 0, bc, 0.0, -f(i));
end
for i = 1 : n
    fprintf(fid,'%d   %d  %d    %f  %f \n',i, i, i+1, E, A(i));
end
fclose(fid);


elementData=[A',L'];
nodeData = [coordinates(:,2),As',f'];

end