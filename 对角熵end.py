import sympy as sp
from sympy.physics.quantum.spin import Jx,Jy,Jz,JzKet,JzBra
# from sympy.physics.quantum import Dagger
from sympy.physics.quantum.qapply import qapply
from sympy import I,sin,cos
from sympy.abc import x
import numpy as np
from scipy.special import jv
import matplotlib.pyplot as plt

# 提前设置符号
t = sp.symbols('t',real=True)
N = sp.symbols('N',real=True,positive=True)
nu = sp.symbols('nu',real=True)
mu = sp.symbols('mu',real=True)
B2,B3,B8,B9,B11,B12 = sp.symbols('B_2,B_3,B_8,B_9,B_{11},B_{12}')
state0 = JzKet(N/2,-N/2) # 初态
x33 = B2*Jy+B3*Jz+B8*Jy*Jz+B9*Jz*Jy+B11*Jy**2+B12*Jz**2
ed = sp.exp(I*x)
rd = 1+I*x-x**2/2-I*x**3/6+x**4/24
# 体系的态
print("开始计算体系的态")
state = qapply((cos(mu)-I*Jx*sin(mu))*(cos(nu)-I*Jy*sin(nu))*rd.subs({x:x33})*state0)
# 展开系数

c0 = qapply(JzBra(N/2,-N/2)*state)
print("计算展开系数：c0")
c1 = qapply(JzBra(N/2,1-N/2)*state)
print("计算展开系数：c1")
c2 = qapply(JzBra(N/2,2-N/2)*state)
print("计算展开系数：c2")
c3 = qapply(JzBra(N/2,3-N/2)*state)
print("计算展开系数：c3")
c4 = qapply(JzBra(N/2,4-N/2)*state)
print("计算展开系数：c4")
c5 = qapply(JzBra(N/2,5-N/2)*state)
print("计算展开系数：c5")
c6 = qapply(JzBra(N/2,6-N/2)*state)
print("计算展开系数：c6")
c7 = qapply(JzBra(N/2,7-N/2)*state)
print("计算展开系数：c7")
c8 = qapply(JzBra(N/2,8-N/2)*state)
print("计算展开系数：c8")
c9 = qapply(JzBra(N/2,9-N/2)*state)
print("计算展开系数：c9")
c10 = qapply(JzBra(N/2,10-N/2)*state)
print("计算展开系数：c10")
# 开始计算对角熵
print("开始计算对角熵")
cs = [c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10]
S = 0
num = 0
for c in cs:
    p =c*sp.conjugate(c) 
    S+=p*sp.log(p)/sp.log(2)

    print("num={}".format(num))
    num+=1
# sp.init_printing()
# print(S)
# 尝试开始替换 # 这里我已经将第三方库里面的hbar替换成了1了
# ================================================= 变量替换
Delta = 1
N=10
l = 2 # lambda=2
v0=1
T=0.1*np.pi
w0=1
t = np.linspace(0,T,100)
v = np.pi
u = np.pi-v0*(t/2-(T*np.sin(2*np.pi*t/T))/(4*np.pi))

parameter = -v0*T/4/np.pi
j_0_p = jv(0,parameter)
j_1_p = jv(1,parameter)
j_0_2p = jv(0,2*parameter)
j_1_2p = jv(1,2*parameter)
B2 = 4*Delta*(j_0_p/v0+8*np.pi*j_1_p/(16*np.pi**2-v0**2*T**2))*np.sin(v0*T/4)**2
B3 = -2*Delta*(j_0_p/v0+8*np.pi*j_1_p/(16*np.pi**2-v0**2*T**2))*np.sin(v0*T/2)
B8 = l/N*(j_0_2p/v0+4*np.pi*T*j_1_2p/(4*np.pi**2-v0**2*T**2))



