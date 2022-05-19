clear all
close all
clc

%% 
img_ref = imread("ref_1.png");
[~,~,P] = readBarcode(img_ref,"QR-CODE");

img_des = imread('des.png');
[n,m,~] = size(img_des);

%%
X1 = P(1,:)';
X2 = P(2,:)';
X3 = P(3,:)';

x1 = [1 n]';
x2 = [1 1]';
x3 = [m 1]';

% Parametros para optimizar
dx = 0;
dy = 0;
theta = 0;
s = 0.3;

q = [dx dy theta s]';

xp1 = Transformacion_Similitud(q,x1);
xp2 = Transformacion_Similitud(q,x2);
xp3 = Transformacion_Similitud(q,x3);

e1 = Distancia_Euclidiana(X1,xp1);
e2 = Distancia_Euclidiana(X2,xp2);
e3 = Distancia_Euclidiana(X3,xp3);

f = (1/6)*(e1+e2+e3);


Imprimir_Imagenes(q,img_des,img_ref)


%% Funciones
function xp = Transformacion_Similitud (qi,xi)
    dx = qi(1);
    dy = qi(2);
    theta = qi(3);
    s = qi(4);
    
    xp = [s*cos(theta) -s*sin(theta); s*sin(theta) s*cos(theta)]*xi + [dx dy]';
end

function e = Distancia_Euclidiana (X,x)
    e = sqrt((X(1)-x(1))^2+(X(2)-x(2))^2);
end

function Imprimir_Imagenes (q,img_des,img_ref)
    dx = q(1);
    dy = q(2);
    theta = q(3);
    s = q(4);
    
    T = [s*cos(theta) -s*sin(theta) dx; s*sin(theta) s*cos(theta) dy; 0 0 1];
    Tp = projective2d(T');
    
    [N,M,~] = size(img_ref);
    [n,m,~] = size(img_des);

    panoramaView = imref2d([N M]);
    Iwarp = imwarp(img_des,Tp,'OutputView',panoramaView);
    Imask = imwarp(true(n,m),Tp,'OutputView',panoramaView);
    
    blender = vision.AlphaBlender('Operation','Binary mask','MaskSource','Input port');
    panorama = step(blender,img_ref,Iwarp,Imask);
    
    imshow(panorama)
end
