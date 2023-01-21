import numpy as np
import matplotlib.pyplot as plt

rho=800
rho0=1000
v0=-4.43
eta=10e-3
r=0.02
g=9.81
V=4/3*np.pi*r**3

a=6*np.pi*eta*r/(rho*V)
b=g*(rho0-rho)/rho

def acc(v,z):
    return -a*v+b


N=50000
dt=0.001
z=np.zeros(N)
z[0]=0
vz=np.zeros(N)
vz[0]=v0
t=np.zeros(N)

for i in range (0,N-1):
    k1=dt*vz[i]
    w1=dt*acc(vz[i],z[i])

    k2=dt*(vz[i]+w1/2)
    w2=dt*acc(vz[i]+w1/2,z[i]+k1/2)

    z[i+1]=z[i]+k2
    vz[i+1]=vz[i]+w2
    t[i+1]=t[i]+dt


def vzt(t):
    return (v0-b/a)*np.exp(-a*t)+b/a

def zt(t):
    return b/a*t+(1/a)*(b/a-v0)*(np.exp(-a*t)-1)

t0=-1/a*np.log(b/(b-a*v0))
print(zt(t0))

fig,(ax1,ax2)=plt.subplots(2)
ax1.plot(t,z)
ax2.plot(t,vz)
ax2.axhline(y=b/a,xmin=0,xmax=t[N-1],color='r')
plt.show()

