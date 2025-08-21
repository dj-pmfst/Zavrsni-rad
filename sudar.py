import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import matplotlib.animation as animation

class Object:
    def __init__(self):
        self.t = []
        self.x = []
        self.y = []
        self.vx = []
        self.vy = []
        self.a = []
        self.kut = []
        self.c = []
        self.m = []

    def set_initial_conditions(self, v0, x0, y0, kut, r, m, b):
        self.c.append(b)
        self.v0 = v0
        self.kut.append(kut)
        self.m.append(m)
        self.r = [r]
        self.t.append(0)
        self.x.append(x0)
        self.y.append(y0)
        self.vx.append(v0 * np.cos(kut * np.pi/180))
        self.vy.append(v0 * np.sin(kut * np.pi/180))
        self.a.append(0)
        return self


class Collision:
    def __init__(self,t, dt):
        self.t  = [0]
        self.t_uk = t
        self.dt = dt
        self.xi = []
        self.yi = []
        self.xj = []
        self.yj = []
        self.objects1 = []  #lijevi
        self.objects2 = []  #desni
        self.tsudar = [0]
        self.sudar_test = 0 

    def add_object(self, objekt, x):
       if x == 1:
           self.objects1.append(objekt)

       elif x == 2:
           self.objects2.append(objekt)

    def __move(self):
        i = self.objects1[0]
        j = self.objects2[0]

        if np.sqrt((i.x[-1] - j.x[-1])**2 + (i.y[-1] - j.y[-1])**2) <= i.r[0] + j.r[0]:   #ne-centralni sudar
            if self.sudar_test == 0:  
                self.__angle()
                i.r.append(0*i.r[0]/2)
                j.r.append(0*j.r[0]/2)
                self.sudar_test += 1 
                for o in self.objects1:
                    v2 = (o.v0*(o.m[0]-j.m[0])+2*j.m[0]*j.v0)/(o.m[0]+j.m[0])
                    o.vx.append(abs(o.v0) * np.cos(o.kut[-1] )) #Kut rapršenja određuje smjer zbog čega je brzina u apsolutnoj vrijednosti
                    o.vy.append(abs(o.v0) * np.sin(o.kut[-1] ))
                for o in self.objects2:
                    v2 = (o.v0*(o.m[0]-i.m[0])+2*i.m[0]*i.v0)/(o.m[0]+i.m[0])
                    o.vx.append(abs(o.v0) * np.cos(o.kut[-1] ))
                    o.vy.append(abs(o.v0) * np.sin(o.kut[-1] ))
            else:
                pass

        if self.sudar_test != 0:
            self.tsudar.append(self.t[-1])
            for o in self.objects1:
                o.vx.append(o.vx[-1])
                o.vy.append(o.vy[-1])

            for o in self.objects2:
                o.vx.append(o.vx[-1])
                o.vy.append(o.vy[-1])
        
        if self.sudar_test == 0:
            self.tsudar.append(0)

        for o in self.objects1:
            o.t.append(o.t[-1] + self.dt)
            o.x.append(o.x[-1] + o.vx[-1]*self.dt)
            o.y.append(o.y[-1] + o.vy[-1]*self.dt)
            o.a.append(0)

        for o in self.objects2:
            o.t.append(o.t[-1] + self.dt)
            o.x.append(o.x[-1] + o.vx[-1]*self.dt)
            o.y.append(o.y[-1] + o.vy[-1]*self.dt)
            o.a.append(0)

        self.t.append(self.t[-1] - self.dt)


    def __angle(self):
        i = self.objects1[0]  ##lijevo
        j = self.objects2[0]  ##desno

        '''spremanje koordinata u trenutku sudara'''
        self.xi.append(i.x[-1])
        self.yi.append(i.y[-1])
        self.xj.append(j.x[-1])
        self.yj.append(j.y[-1])

        '''x koordinate u slucaju raspada'''
        if len(self.objects1) >= 2:
            for o in self.objects1:
                o.x[0] += self.xi[0]
            for o in self.objects2:
                o.x[0] += self.xj[0]
        
        for e in self.objects1:

            del self.xj[1:len(self.xj)]
            del self.yj[1:len(self.yj)]

            sudar = 0
            while sudar == 0:
                if e.x[0] + e.r[0] <= self.xj[-1] - j.r[0] and e.r[0] + j.r[0] < abs(e.y[0] - self.yj[0]):
                    sudar = 2 
                elif np.sqrt((e.x[0] - self.xj[-1])**2 + (e.y[0] - self.yj[-1])**2) <= e.r[0] + j.r[0]:
                    sudar = 1
                else:
                    self.xj.append(self.xj[-1] + j.vx[-1]*self.dt)
                    self.yj.append(self.yj[-1] + j.vy[-1]*self.dt)               

            if sudar == 1:
                a = (self.yj[-1] - e.y[0])/(self.xj[-1] - e.x[0]) #nagib pravca koji prolazi kroz srediste obe kruznice
                theta = np.arctan(a) #vraca kut uvijek izmedu 0 i 90° 
                kut_i = np.pi - 2*theta  #iznos kuta
                if e.y[-1] < self.yj[-1]: #provjera kvadranta
                    if abs(kut_i) < np.pi:  #npr. ako je lijeva kuglica ispod desne kut je u rasponu [pi,2pi] ili [o,-pi]
                        kut_i = np.negative(kut_i)
                    else:
                        kut_i = np.positive(kut_i)
                else:
                    if abs(kut_i) < np.pi:
                        kut_i = np.positive(kut_i)
                    else:
                        kut_i = np.negative(kut_i)
                if kut_i < 0:
                    kut_j = np.pi + kut_i
                else:
                    kut_j = np.pi - kut_i
                e.kut.append(kut_i)

        for e in self.objects2:

            del self.xi[1:len(self.xi)]
            del self.yi[1:len(self.yi)]

            sudar = 0
            while sudar == 0:
                if e.x[0] - e.r[0] <= self.xi[-1] + i.r[0] and e.r[0] + i.r[0] < abs(e.y[0] - self.yi[0]):
                    sudar = 2
                elif np.sqrt((self.xi[-1] - e.x[0])**2 + (self.yi[-1] - e.y[0])**2) <= i.r[0] + e.r[0]:
                    sudar = 1
                else:
                    self.xi.append(self.xi[-1] + i.vx[-1]*self.dt)
                    self.yi.append(self.yi[-1] + i.vy[-1]*self.dt)

            if sudar == 1:
                a = (self.yi[-1] - e.y[0])/(self.xi[-1]-e.x[0])
                theta = np.arctan(a) 
                kut_i = np.pi - 2*theta 
                if self.yi[-1] < e.y[-1]:
                    if abs(kut_i) < np.pi:
                        kut_i = np.negative(kut_i)
                    else:
                        kut_i = np.positive(kut_i)
                else:
                    if abs(kut_i) < np.pi:
                        kut_i = np.positive(kut_i)
                    else:
                        kut_i = np.negative(kut_i)
                if kut_i < 0:
                    kut_j = np.pi + kut_i
                else:
                    kut_j = np.pi - kut_i
                e.kut.append(kut_j)

    def impact(self):
        while self.t[-1] > -self.t_uk:
            self.__move()   
        kut1 = []
        kut2 = []
        for o in self.objects1:
            kut1.append(o.kut[-1]*57.2958)
        for o in self.objects2:
            kut2.append(o.kut[-1]*57.2958)
        
        del kut1[0]
        del kut2[0]

        ##------------------------------------------------------------------
        i = self.objects1[0]
        j = self.objects2[0]
        c = 299792458
        ev = 6.242e18

        prvo =  i.m[0]*i.v0 + j.m[0]*abs(j.v0)
        drugo = 0

        prvoe = ((i.m[0]*i.v0**2)/2 + (j.m[0]*j.v0**2)/2)
        drugoe = 0

        objects_1 = []
        objects_2 = []
        for o in self.objects1:
            if abs(o.kut[-1]) >= 1:
                objects_1.append(o)
        for o in self.objects2:
            if abs(o.kut[-1]) >= 1:
                objects_2.append(o)
        
        return kut1, kut2    

    def plot(self):
        fig, axs = plt.subplots()
        r = 0.87*10**(-15)

        ii = self.objects1[0]
        jj = self.objects2[0]

        while self.t[-1] > -self.t_uk:
            self.__move()

        objects = []
        for o in self.objects1:
            objects.append(o)
        for o in self.objects2:
            objects.append(o)

        T = len(self.t)
        rt = []

        for o in objects:
            for i in range(T):
                plt.grid()
                if self.tsudar[i] == 0:
                    plt.plot(ii.x[0:i], ii.y[0:i], color=ii.c[0], linewidth=2)
                    plt.plot(jj.x[0:i], jj.y[0:i], color=jj.c[0], linewidth=2)
                if self.tsudar[i] != 0:
                    rt.append(i)
                    if o == self.objects1[0] or o == self.objects2[0]:
                        pass
                    else:
                        plt.plot(o.x[i:-1],o.y[i:-1], color=o.c[0])

            circle = plt.Circle((ii.x[rt[0]], ii.y[rt[0]]), radius=ii.r[0], color = ii.c[0], fill = False)
            circle2 = plt.Circle((jj.x[rt[0]], jj.y[rt[0]]), radius=jj.r[0], color=jj.c[0], fill=False)
            if o == self.objects1[0] or o == self.objects2[0]:
                pass 
            else:
                circle3 = plt.Circle((o.x[rt[-1]], o.y[rt[-1]]), radius=o.r[0], color=o.c[0], fill=True)
                axs.add_patch(circle3)
            axs.add_patch(circle)
            axs.add_patch(circle2)

        plt.grid()
        plt.show()

    def animate(self):
        fig, axs = plt.subplots()
        fig = plt.figure()
        writer = animation.PillowWriter(50)
        r = 0.87*10**(-15)  

        ii = self.objects1[0]
        jj = self.objects2[0]

        while self.t[-1] > -self.t_uk:
            self.__move()

        objects = []
        for o in self.objects1:
            objects.append(o)
        for o in self.objects2:
            objects.append(o)
        n = []
        for o in objects:
            n.append(len(o.x))

        rt = 0

        with writer.saving(fig, "sudar86.gif", 100):
            plt.xlim(-3*r,3*r)
            plt.ylim(-3*r,3*r)
            plt.grid()
            for i in range(n[1]):
                if i%1 == 0:
                    plt.clf()
                    for o in objects:
                        axs = plt.subplot()
                        plt.xlim(-3*r,3*r)
                        plt.ylim(-3*r,3*r)
                        x = o.x
                        y = o.y
                        plt.grid()
                        plt.plot(x[:i], y[:i], linewidth = 0.6)
                        if self.tsudar[i] != 0:
                            rt = -1
                        circle = plt.Circle((ii.x[i], ii.y[i]), radius=ii.r[rt], color = ii.c[0], fill = True)
                        axs.add_patch(circle)
                        circle = plt.Circle((jj.x[i], jj.y[i]), radius=jj.r[rt], color = jj.c[0], fill = True)
                        axs.add_patch(circle)
                        if rt == -1:
                            circle = plt.Circle((o.x[i], o.y[i]), radius=o.r[rt], color = o.c[0], fill = True)
                            axs.add_patch(circle)
                        #plt.scatter(o.x[i], o.y[i], color=o.c[0])
                        # plt.scatter(x[i], y[i], s = 20000/8)
                    writer.grab_frame()
