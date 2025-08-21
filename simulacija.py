import sudar2 as s 
import matplotlib.pyplot as plt 
import math
import numpy as np 
import random 
from matplotlib.animation import PillowWriter

v_c = 299792458
v = 0.99*v_c
r1 = 0.87*10**(-15)
r2 = 0.87*10**(-15)
m = 1.6726*10**(-27)

def start(r1,r2, n):
    x1 = -10**(-15)
    y1 = 0
    x2 = 10**(-15)
    y2 = 0.2*r1
    X = [x1, x2]
    Y = [y1,y2]

    '''gaussova raspodjela'''
    r = abs(r1+0*0.1*10**(-15))
    l = -abs(r1+0*0.1*10**(-15))
    sr = (l+r)/2
    x = np.arange(-r1, r2)
    sm = 0
    for i in range (len(x)): 
        sm += (x[i]-sr)**2
    
    st_dev = 2*r1/6
    y = np.random.normal(sr, st_dev)

    '''stvaranje čestica'''
    a = s.Object()  

    '''                      v,  x,  y, kut, r,  m'''
    a.set_initial_conditions(v, x1, y1, 0, r1, m, 'blue')

   
    b = s.Object() 
    b.set_initial_conditions(-v, x2, y, 0, r2, m, 'red')
    
    c = s.Collision(10*10**(-24), 10**(-26))  
    c.add_object(a, 1)
    c.add_object(b, 2)

    '''određivanje položaja i stvaranje "kvarkova" '''
    xx = []
    yy = []
    r = r1/(1+(1/np.sin(np.pi/n)))

    kut_promjena = 2 * np.pi / n
    
    for i in range(n):
        angle = i * kut_promjena
        x_ = (r - r1) * np.cos(angle)
        y_ = (r - r1) * np.sin(angle)
        xx.append(x_)
        yy.append(y_)
    
    if n == 2:
        l = (abs(xx[1]-xx[0]))/2
        xx[0] = 0
        xx[1] = 0
        yy[0] = l
        yy[1] = -l 

    for i in range (n):
        d = s.Object()
        d.set_initial_conditions(v, xx[i]+X[0], yy[i]+Y[0], 0, r, m/(n), 'blue')
        c.add_object(d, 1)
        r = r1/(1+(1/np.sin(np.pi/n)))
        e = s.Object()
        e.set_initial_conditions(-v, xx[i]+X[1], yy[i]+y, 0, r, m/(n), 'red')
        c.add_object(e, 2)

    kut = c.impact()
    #c.plot()
    #c.animate()

    return kut, y

clr = ['blue', 'orange', 'green']
l_clr = ['lightblue', 'navajowhite', 'lightgreen']

for i in range(3):
    
    N = 14000/2**i 
    korak = 0
    kut1 = []
    kut2 = []
    kut_uk = []
    Y = []

    while korak < N:
        n = random.randint(2,8)
        a, b = start(r1, r2, n)
        kut1 += a[0]  #kutevi raspršenja čestica s lijeve strane
        kut2 += a[1]  #-||- s desne strane
        Y.append(b)
        korak += 1

    kut_uk = kut1 + kut2

    psdr = []
    theta= []
    for k in kut_uk:
        theta_rad = math.radians(k)
        theta.append(theta_rad)
        tan_half = math.tan(theta_rad / 2)
        if tan_half > 0:
            psdr.append(-math.log(tan_half))

    '''histogram'''
    bins = np.linspace(-5, 5, 51)  
    hist, edges = np.histogram(psdr, bins=bins)
    bin_centers = 0.5 * (edges[1:] + edges[:-1])
    bin_width = edges[1] - edges[0]


    dN_deta = hist / (bin_width)


    plt.plot(bin_centers, dN_deta, linestyle='--', color=l_clr[int(i)], zorder=1)  
    plt.plot(bin_centers, dN_deta, linestyle='None', marker='o', markersize=4, color=clr[int(i)], zorder=2) 
plt.xlabel(r"$\eta$")
plt.ylabel('N')
plt.grid()
plt.show()

# Optional: look at event-by-event multiplicity around η≈0 (experiments often quote this)
# print("⟨N_ch(|η|<0.5)⟩ per event =", np.mean(mult_eta0))
# plt.figure()
# plt.hist(mult_eta0, bins=40)
# plt.xlabel(r'$N_{\rm ch}$ in $| \eta | < 0.5$ per event'); plt.ylabel('Events')
# plt.grid(True); plt.show()
