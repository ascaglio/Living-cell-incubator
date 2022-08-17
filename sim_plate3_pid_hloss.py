#Simulation using a PID control system

import numpy as np
import matplotlib.pyplot as plt 
import math
import pandas as pd
import random
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
import matplotlib.patches as patches
from sklearn.metrics import mean_squared_error



alpha=21 # (mm^2/s) thermal diffusivity, moist air
alpha2=18.8 # (mm^2/s) thermal diffusivity, stainless steel
#length of sides  of heated square plate
Lx=30 #(mm)
Ly=30 #(mm)
#number of points
N=30 #number of points in every dimension 

#Discretize my space
Xvec=np.linspace(0,Lx,N)
Yvec=np.linspace(0,Ly,N)
dx=Xvec[2]-Xvec[1] #dimensional step in x direction
dy=Yvec[2]-Yvec[1] ##dimensional step in y direction

#Discretize time
#dt=0.5*(dx**2)/(2*alpha) 
dt=0.01
time=3601 #how long i want to run the simulations (sec)
tvec=np.linspace(0,1,100) #this is how long i run the numerical approximation (sec)

#Inital Boundary Conditions

T=np.full((Xvec.size,Yvec.size),25.0) #entire space is at 25 C initially

T[0,10:20]=50.0 #50 degrees C applied to the plate
T[29,10:20]=50.0


#2nd order polynomial fit to predict deltaT in the grid
train = pd.read_excel('train_pid_loss.xlsx')
X, y = train[["T1_train", "T2_train", "T3_train"]], train["y_train"]

poly_reg = PolynomialFeatures(degree = 2)
X_poly = poly_reg.fit_transform(X)
 
pol_reg = LinearRegression()
pol_reg.fit(X_poly, y)


#PID function
'''
def PID(error,prev_error):
    kp = 4
    ki = 20
    kd = 2
    cum_error = error*1  #integral
    rate_error = error - prev_error   #derivative
    output = kp*error + ki*cum_error + kd*rate_error
    return(output)
'''

#NUMERICAL SOLVER FOR THE HEAT EQUATION

Tout=[]
diff_actual = []
diff_pred = []
Tgrid = []
temp1 = []
temp2 = []
temp3 = []
rmse = []
Tresistor = []
target = np.full((10,10),37)
actual = np.full((10,10),25)
pid_output = []
Tm = []
random = [random.uniform(-0.2,0.2) for x in range(time)]
Kp = 1
Ki = 0
Kd = 0
dt_pid = 1
cum_error = -12

for i in range (0,time):

     locals()["T" + str(i)]=np.full((Xvec.size,Yvec.size),25.0)
     locals()["T" + str(i)][0,10:20]=50.0
     locals()["T" + str(i)][29,10:20]=50.0

for i in range (0,time):
    
    if (i==0): 
        Told = T
    else:
        Told = locals()["T" + str(i-1)]
    
    for tx in range(1,Xvec.size-1):
        
        for ty in range(1,Yvec.size-1):

            if(i==0):
                du1=(dt*(alpha2*(Told[0,ty+1]-2*Told[0,ty]+Told[0,ty-1])/dy**2)+Told[0,ty])-T[0,ty]
                T[0,ty] = T[0,ty] + du1
                du2=(dt*(alpha2*(Told[29,ty+1]-2*Told[29,ty]+Told[29,ty-1])/dy**2)+Told[29,ty])-T[29,ty]
                T[29,ty] = T[29,ty] + du2
                
                du=(dt*(alpha*(Told[tx+1,ty]-2*Told[tx,ty]+Told[tx-1,ty])/dx**2+alpha*(Told[tx,ty+1]-2*Told[tx,ty]+Told[tx,ty-1])/dy**2)+Told[tx,ty])-T[tx,ty]
                T[tx,ty] = T[tx,ty] + du
                locals()["T" + str(i)] = T
            
            else:
                du3=(dt*(alpha2*(Told[0,ty+1]-2*Told[0,ty]+Told[0,ty-1])/dy**2)+Told[0,ty])-locals()["T" + str(i-1)][0,ty]
                locals()["T" + str(i)][0,ty] = locals()["T" + str(i-1)][0,ty] + du3
                du4=(dt*(alpha2*(Told[29,ty+1]-2*Told[29,ty]+Told[29,ty-1])/dy**2)+Told[29,ty])-locals()["T" + str(i-1)][29,ty]
                locals()["T" + str(i)][29,ty] = locals()["T" + str(i-1)][29,ty] + du4    
                
                du=(dt*(alpha*(Told[tx+1,ty]-2*Told[tx,ty]+Told[tx-1,ty])/dx**2+alpha*(Told[tx,ty+1]-2*Told[tx,ty]+Told[tx,ty-1])/dy**2)+Told[tx,ty])-locals()["T" + str(i-1)][tx,ty]
                locals()["T" + str(i)][tx,ty] = locals()["T" + str(i-1)][tx,ty] + du
                
    if(i>0):
        if(locals()["T" + str(i)][15,1]>25):
            locals()["T" + str(i)][1:28,1] = locals()["T" + str(i)][1:28,1]-(locals()["T" + str(i)][1:28,1])*0.001 
            locals()["T" + str(i)][1:28,29] = locals()["T" + str(i)][1:28,29]-(locals()["T" + str(i)][1:28,29])*0.001
        if(locals()["T" + str(i)][0,5]>25):
            locals()["T" + str(i)][0,0:9] = locals()["T" + str(i)][0,0:9]-(locals()["T" + str(i)][0,0:9])*0.001
            locals()["T" + str(i)][0,21:29] = locals()["T" + str(i)][0,21:29]-(locals()["T" + str(i)][0,21:29])*0.001
            locals()["T" + str(i)][29,0:9] = locals()["T" + str(i)][29,0:9]-(locals()["T" + str(i)][29,0:9])*0.001
            locals()["T" + str(i)][29,21:29] = locals()["T" + str(i)][29,21:29]-(locals()["T" + str(i)][29,21:29])*0.001

    locals()["T" + str(i)][:,0] = locals()["T" + str(i)][:,1]
    locals()["T" + str(i)][:,29] = locals()["T" + str(i)][:,28]
    
    if(i==0):
    
        temp1.append(locals()["T" + str(i)][5,15] + random[i])
        temp2.append(locals()["T" + str(i)][5,5] + random[i])
        temp3.append(locals()["T" + str(i)][15,5] + random[i])
    else:
        temp1.append(locals()["T" + str(i-1)][5,15] + random[i])
        temp2.append(locals()["T" + str(i-1)][5,5] + random[i])
        temp3.append(locals()["T" + str(i-1)][15,5] + random[i])
    
    temp = np.array([[temp1[i],temp2[i],temp3[i]]])
    actual = locals()["T" + str(i)][10:20,10:20] 
    diff = np.subtract(actual,target)
    mse = mean_squared_error(actual,target)
    root_mse = np.sqrt(mse)
    rmse.append(root_mse)
    deltaT_actual = diff.mean() 
    deltaT_pred = pol_reg.predict(poly_reg.fit_transform(temp))
    diff_pred.append(deltaT_pred)
    diff_actual.append(deltaT_actual)
    Ta = (locals()["T" + str(i)][28,15]*10 + locals()["T" + str(i)][29,9] + locals()["T" + str(i)][29,21] + 10*25)/22 
        
    if (i<dt_pid):
        locals()["T" + str(i)][0,10:20]=50.0
        locals()["T" + str(i)][29,10:20]=50.0
    else:
        #out_control = PID(diff_pred[i],diff_pred[i-1])
        #cum_error += diff_pred[i]*dt_pid  #integral
        #rate_error = (diff_pred[i] - diff_pred[i-dt_pid])/dt_pid   #derivative
        #out_control = Kp*diff_pred[i] + Ki*cum_error + Kd*rate_error
        cum_error += diff_actual[i]*dt_pid  #integral
        rate_error = (diff_actual[i] - diff_actual[i-dt_pid])/dt_pid   #derivative
        out_control = Kp*diff_actual[i] + Ki*cum_error + Kd*rate_error
        resistor_control = np.full(10,out_control)
        if (out_control<0):    
            #locals()["T" + str(i)][0,10:20]=50.0
            #locals()["T" + str(i)][29,10:20]=50.0
            locals()["T" + str(i)][0,10:20] = locals()["T" + str(i-1)][0,10:20] - resistor_control
            locals()["T" + str(i)][29,10:20] = locals()["T" + str(i-1)][29,10:20] - resistor_control
            if (locals()["T" + str(i)][0,11]>50):
                locals()["T" + str(i)][0,10:20] = 50
                locals()["T" + str(i)][29,10:20] = 50
        else:
            if (locals()["T" + str(i)][0,11]>25):
                locals()["T" + str(i)][0,10:20]= Ta + (locals()["T" + str(i-1)][0,15] - Ta)*math.exp(-0.052) 
                locals()["T" + str(i)][29,10:20]= Ta + (locals()["T" + str(i-1)][29,15] - Ta)*math.exp(-0.052)
        pid_output.append(out_control)    
    Tout.append(locals()["T" + str(i)]) 
    Tresistor.append(locals()["T" + str(i)][0,11])
    Tm.append(Tout[i][15,15])

actual_diffT = []

for i in range(0,len(diff_actual)):
    e = abs(diff_actual[i])
    actual_diffT.append(e)



'''
#GUARDAR DATOS

df =pd.DataFrame()
df['T1'] = temp1
df['T2'] = temp2
df['T3'] = temp3
df['deltaT_actual'] = diff_actual

df['actual_dt'] = diff_actual
#df['predicted_dt'] = diff_pred
df.to_excel('diff_actual_hloss.xlsx', index=False)


'''

#GRAPHICS       

'''
    
X=np.linspace(0,31,31)
Y=np.linspace(0,31,31)

for i in range(0,500):
    fig, ax = plt.subplots()
    im = ax.pcolormesh(X, Y, Tout[i], vmin=0, vmax=60, cmap='jet')
    plt.hlines(y=30, xmin=0, xmax=31)
    plt.hlines(y=1, xmin=0, xmax=31)
    plt.plot(15,5, marker="o", color="black")
    ax.text(16,4,'T1=',fontsize=12)
    ax.text(19,4,round(Tout[i][5,15],1),fontsize=12)
    ax.text(22,4,'°C',fontsize=12)
    plt.plot(5,5, marker="o", color="black")
    ax.text(2,2,'T2=',fontsize=12)
    ax.text(5,2,round(Tout[i][5,5],1),fontsize=12)
    ax.text(8,2,'°C',fontsize=12)
    plt.plot(5,15, marker="o", color="black")
    ax.text(1,12,'T3=',fontsize=12)
    ax.text(4,12,round(Tout[i][15,5],1),fontsize=12)
    ax.text(7,12,'°C',fontsize=12)
    ax.add_patch(
     patches.Rectangle(
        (10, 10),
        10,
        10,
        edgecolor = 'black',
        fill=False
     ) )
    for j in range(10,20):
        plt.hlines(y=j, xmin=10, xmax=20)
        plt.vlines(x=j, ymin=10, ymax=20)
    fig.colorbar(im)
    ax.text(25, 27,i,fontsize=20)
    ax.text(25, 24.5,'(sec)',fontsize=20)
    ax.text(1,27,'RMSE(°C):',fontsize=15)
    ax.text(1,24.5,round(rmse[i],2),fontsize=15)
    plt.xlabel('X(cm)',fontsize=20)
    plt.ylabel('Y(cm)',fontsize=20)
    plt.savefig('fig'+str(i)+'.jpg')
    plt.close()

'''


#Surface of error vs time

'''
t_pos=1
diff_cum_pos = 0
t_neg = 1
diff_cum_neg = 0
areas = []

for i in range(0,time):
    if(diff_actual[i]<0):
       areas.append(diff_cum_pos/t_pos)
       diff_cum_pos = 0
       t_pos = 1
       t_neg+=1
       diff_cum_neg += diff_actual[i]
    else:
        areas.append(diff_cum_neg/t_neg)
        diff_cum_neg = 0
        t_neg = 1
        t_pos+=1
        diff_cum_pos += diff_actual[i]
        
'''

#Amplitud error

t_pos=1
diff_cum_pos = [0]
t_neg = 1
diff_cum_neg = [0]
areas = []

for i in range(0,time):
    if(diff_actual[i]<0):
       areas.append(max(diff_cum_pos))
       diff_cum_pos = [0]
       t_pos = 1
       t_neg+=1
       diff_cum_neg.append(diff_actual[i])
    else:
        areas.append(min(diff_cum_neg))
        diff_cum_neg = [0]
        t_neg = 1
        t_pos+=1
        diff_cum_pos.append(diff_actual[i])

  
error_max = round(max(diff_actual[3300:3600]),2)
error_min = round(min(diff_actual[3300:3600]),2)
  
Tsp = np.full((time),37.0)
cero = np.full((time),0)
timeVec = [] 
t=0
for i in range(0,time):
    t=t+1
    timeVec.append(t)

for i in range(0,dt_pid):
    pid_output.insert(i-1,0)

fig, axs = plt.subplots(4,1)
axs[0].plot(timeVec, diff_actual,label='$\Delta T_a$')
axs[0].plot(timeVec,diff_pred,alpha=0.5,label='$\Delta T_p$')
axs[0].plot(timeVec,cero,'--')
axs[0].set(ylabel='T($\degree$C)')
axs[0].set(ylim=(-3,3))
axs[0].set(xticks=[])
axs[0].legend(loc='upper left')
axs[1].stem(timeVec,areas,markerfmt='_',label='max e/cycle')
axs[1].set(ylabel='T($\degree$C)')
axs[1].legend(loc='upper left')
axs[1].set(ylim=(-1.5,1.5))
axs[1].set(xticks=[])
axs[1].text(0,-0.6,' '+str(error_max)+'$\degree$C',fontsize=8)
axs[1].text(0,-1.2,str(error_min)+'$\degree$C',fontsize=8)
axs[2].plot(timeVec,pid_output,color='green',label='PID')
axs[2].set(ylim=(-5,5))
axs[2].text(700,-2,'Kp='+str(Kp)+'  Ki='+str(Ki)+'  Kd='+str(Kd)+'  dt='+str(dt_pid)+'sec')
axs[2].set(ylabel='T($\degree$C)')
axs[2].plot(timeVec,cero,'--')
axs[2].set(xticks=[])
axs[2].legend(loc='upper left')
axs[3].plot(timeVec,Tresistor,color='black',label='Tr')
axs[3].set(xlabel='time(sec)', ylabel='T($\degree$C)')
axs[3].legend(loc="lower left")
#axs[3].set(ylim=(24,51))


'''

plt.plot(timeVec,Tm,label='T in center')
plt.plot(timeVec,Tsp,'--')
plt.ylabel('T($\degree$C)')
plt.xlabel('Time(sec)')
plt.legend()

'''
