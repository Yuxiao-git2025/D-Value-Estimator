# D-Value-Estimator
An optimized D-value estimator implemented in MATLAB

```matlab
clc;clear;
load("PackageSAF_Fractal.mat");
time=timen;
time=time-time(1);
mag=magn;
lon=lonn;
lat=latn;
dep=depn;
clear lonn latn magn depn timen
%%
% From 500m to 80km
% window size is designed by yourself
r1=0.5*1e3;
r2=80*1e3;
Num=35;
win=80;
dwin=10;
dM=0.1;
loc=[lon,lat];
loctest=loc(100:150,:);
% Do some check
[D,Dstd,Methods,R2]=cal_Dvalue2(loctest,r1,r2,Num,1);
```

<div align=center><img width="700" height="500" alt="NCEDC" src="https://github.com/user-attachments/assets/1251256a-cb90-493b-ba93-449a29664c7c" /></div>

<div align=center><img width="700" height="500" alt="Sup12Select" src="https://github.com/user-attachments/assets/594880e1-62bd-4722-8e4c-186eb89a49b0" /></div>
