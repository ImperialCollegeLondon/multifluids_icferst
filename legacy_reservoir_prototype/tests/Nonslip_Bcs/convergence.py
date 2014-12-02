from fluidity_tools import stat_parser as stat
import pylab as p
import matplotlib as mpl

mpl.rcParams['font.size']=14



def stats(fname):
    error=stat(fname+".stat")["Fluid"]["Error%magnitude"]["l2norm"][:]
    t=stat(fname+".stat")["ElapsedTime"]["value"][:]
    return t,error

def statsMP(fname):
    error=stat(fname+".stat")["phase1"]["Error%magnitude"]["l2norm"][:]
    t=stat(fname+".stat")["ElapsedTime"]["value"][:]
    return t,error


h=(0.2,0.1,0.05,0.025,0.0125)


f=('A','B','C','D','E')

error=[]
errorMp=[]

for H,F in zip(h,f):
    t,e= stats('newtonian_chris_'+F)
    
    error.append(e[-1])
    print t[-1],e[-1]

h0=[error[0]/(2**i) for i in range(len(h))]
h1=[error[0]/(4**i) for i in range(len(h))]

p.figure()
p.xscale('log',basex=2)
p.yscale('log',basey=2)
p.scatter(h,error,s=60)
p.plot(h,h0,'k--',lw=2)
#p.plot(h,h1,'k:',lw=2)

p.xlabel("Mesh length")
p.ylabel(r"$L_2$ norm of velocity error")

for H,F in zip(h[:-1],f[:-1]):
    t,e= statsMP('newtonian_multiphase_'+F)
    
    errorMp.append(e[-1])
    print t[-1],e[-1]

h0=[errorMp[0]/(2**i) for i in range(len(errorMp))]
h1=[errorMp[0]/(4**i) for i in range(len(errorMp))]

p.figure()
p.xscale('log',basex=2)
p.yscale('log',basey=2)
p.scatter(h[:-1],errorMp,s=60)
p.plot(h[:-1],h0,'k--',lw=2)
#p.plot(h,h1,'k:',lw=2)

t,e= statsMP('newtonian_multiphase_C_adaptive')
ha=[0.085]

p.scatter(ha,[e[-1]],c='r',s=60)
p.xlabel("Mesh length")
p.ylabel(r"$L_2$ norm of velocity error")

