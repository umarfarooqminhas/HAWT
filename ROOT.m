%% DATA
clear all
clc
%% DATA EXTRACTION

%ROOT
s=0.05;

load A_sp4721_200K.txt;
load A_sp4721_500K.txt;
load A_sp4721_1000K.txt;
load A_sp4721_2000K.txt;
load A_sp4721_5000K.txt;

[A_sp4721_200K] = palas(A_sp4721_200K,s);
[A_sp4721_500K] = palas(A_sp4721_500K,s);
[A_sp4721_1000K] = palas(A_sp4721_1000K,s);
[A_sp4721_2000K] = palas(A_sp4721_2000K,s);
[A_sp4721_5000K] = palas(A_sp4721_5000K,s);

%PRIMARY
s=0.05;

load B_sp4621_200K.txt;
load B_sp4621_500K.txt;
load B_sp4621_1000K.txt;
load B_sp4621_2000K.txt;
load B_sp4621_5000K.txt;

[B_sp4621_200K] = palas(B_sp4621_200K,s);
[B_sp4621_500K] = palas(B_sp4621_500K,s);
[B_sp4621_1000K] = palas(B_sp4621_1000K,s);
[B_sp4621_2000K] = palas(B_sp4621_2000K,s);
[B_sp4621_5000K] = palas(B_sp4621_5000K,s);

%TIP
s=0.05;

load C_S828_50K.txt;
load C_S828_100K.txt;
load C_S828_200K.txt;
load C_S828_500K.txt;
load C_S828_1000K.txt;

[C_S828_50K] = palas(C_S828_50K,s);
[C_S828_100K] = palas(C_S828_100K,s);
[C_S828_200K] = palas(C_S828_200K,s);
[C_S828_500K] = palas(C_S828_500K,s);
[C_S828_1000K] = palas(C_S828_1000K,s);

%%



