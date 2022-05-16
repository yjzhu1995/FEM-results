import numpy as np
import math

'''
Diffusion creep flow law from Burov, 2011 table 1 for dry olivine
'''
def burov(tep,str_rate):
    A=7.7e-8
    d=1e-3 #mm
    R=8.314
    H=536e3
    m=3
    strength = str_rate*math.pow(d,m)/A*np.exp(H/(R*tep))
    return strength

def burgmann(tep,str_rate):
    A=math.pow(10,9.2)
    d=1 #micro meter
    H=375e3
    m=3
    strength = str_rate*math.pow(d,m)/A*np.exp(H/(R*tep))    
    return strength

def dislocation(tem,str_rate):
    tem=tem + 273
    A=7
    n=3
    H=523e3
    R=8.314
    strength = (str_rate/A*1e14*math.exp(H/(R*tem)))**(1/n)
    return strength

def burgmann_dislocation(tem,str_rate):
    tem=tem + 273
    A=4.e-15
    n=3
    H=523e3
    R=8.314
    strength = (str_rate/A*math.exp(H/(R*tem)))**(1/n)
    return strength

def dorn(tem,str_rate):
    tem=tem + 273
    A=5.7e11
#   q=2
    H=549e3
    R=8.314    
    p_stress = 8.5e9
    strength = p_stress*(1-math.sqrt(R*tem/H*math.log(A/str_rate)))
    return strength

def envelope(b1,b2,d1):
    if np.size(b1) == np.size(b2) and np.size(b1) == np.size(d1):
        print("Format correct, start enveloping!")
    else:
        print("Format error, please check the input file!")
        return None
    for i in range(len(b1[:,0])):
        if b1[i,1]> d1[i,1]:
            b1[i,1] = d1[i,1]
        if b2[i,1] > d1[i,1]:
            b2[i,1] = d1[i,1]
    return b1,b2




#####parameters of flow for diffusion creep: A is constant (s^-1), mu is rigidity (Pa), d is grain size (m), H is activation athalpy (kj/mol)



'''
Parameters from Karato & Wu, 1993
A = 8.7e15 
mu = 8e10
d = 1e-6
R = 8.314
H = 300e3
m = 2.5
str_rate = 1e-15
b = 0.5e-9
'''

'''
Parameters from Hirth & Kohlstedt, 2003
A = 1.5e9 
mu = 1
d = 1
R = 8.314
H = 375e3
m = 3
str_rate = 1e-14
b = 1
'''
A = 8.7e15 
mu = 8e10
d = 1e-6
R = 8.314
H = 300e3
m = 2.5
str_rate = 1e-15
b = 0.5e-9
## input temperature profile name
in_tp_name = "t.1"
in_tp = np.loadtxt(in_tp_name)


dep, tem = in_tp[:,0], in_tp[:,1]
#press = 
#strength = str_rate * mu * math.pow(d/b,m) / A *np.exp(H/(R*tep))
strength = burgmann(tem+273,str_rate)
out = np.hstack((np.mat(dep).T,np.mat(strength).T))
#out_name = 'diffusion_size' + str(d) + '_rate' + str(str_rate)
out_name = 'test.dat'
np.savetxt(out_name,out)

count = len(dep)
strength = np.zeros(count)
strength_new = strength.copy()
strength_dorn = strength.copy()
for i in range(count):
    strength[i] = dislocation(tem[i], str_rate)
    strength_new[i] = burgmann_dislocation(tem[i], str_rate)
    strength_dorn[i] = dorn(tem[i],str_rate)

out = np.hstack((np.mat(dep).T,np.mat(strength).T))
out_name = 'dislocation.dat'
np.savetxt(out_name,out)
out = np.hstack((np.mat(dep).T,np.mat(strength_new).T))
out_name = 'dislocation_new.dat'
np.savetxt(out_name,out)
print(math.pow(7e-14,1/3))
out = np.hstack((np.mat(dep).T,np.mat(strength_dorn).T))
out_name = 'dorn.dat'
np.savetxt(out_name,out)


in1 = np.loadtxt('b.1')
in2 = np.loadtxt('b.2')
in3 = np.loadtxt('d.1')
bd1, bd2 = envelope(in1,in2,in3)
np.savetxt('bd1',bd1)
np.savetxt('bd2',bd2)