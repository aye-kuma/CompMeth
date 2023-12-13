import numpy as np 
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython import display
import argparse

# Gravitational constant
G = 1

# Sets up argparse, makes a group where one of the options is mandatory 
parser = argparse.ArgumentParser(prog='ThreeBody')
group = parser.add_argument_group('Three Body Problem: Trisolaris', 'Animates solutions to the Three Body problem')
exclusive_group = group.add_mutually_exclusive_group(required=True)
exclusive_group.add_argument('--initial', action='store_true',
                             help='Select your own initial conditions for Three Body problem and see the solution animated')

exclusive_group.add_argument('--known',type=str, choices=['Figure8', 'Butterfly', 'Butterfly III', 'Butterfly IV',
                                                          'Bumblebee', 'Moth', 'Moth II', 'Moth III', 'Goggles',
                                                          'Dragonfly', 'Yin Yang'],
                                                          help='Predefined conditions periodic solutions. Use quotation marks for ones that contain more than one phrase')

args = parser.parse_args()


def ThreeBody(t,r):
    """
    Function that contains the ODEs for the three body problem
    """
    x1,y1,vx1,vy1,x2,y2,vx2,vy2,x3,y3,vx3,vy3 = r
    
    r12 = np.sqrt((x2-x1)**2+(y2-y1)**2)
    r13 = np.sqrt((x3-x1)**2+(y3-y1)**2)
    r23 = np.sqrt((x2-x3)**2+(y2-y3)**2)
    
    dvx1 = -G*mass[1]*(x1-x2)/r12**3 - G*mass[2]*(x1-x3)/r13**3
    dvx2 = -G*mass[2]*(x2-x3)/r23**3 - G*mass[0]*(x2-x1)/r12**3
    dvx3 = -G*mass[0]*(x3-x1)/r13**3 - G*mass[1]*(x3-x2)/r23**3
    
    dvy1 = -G*mass[1]*(y1-y2)/r12**3 - G*mass[2]*(y1-y3)/r13**3
    dvy2 = -G*mass[2]*(y2-y3)/r23**3 - G*mass[0]*(y2-y1)/r12**3
    dvy3 = -G*mass[0]*(y3-y1)/r13**3 - G*mass[1]*(y3-y2)/r23**3
    
    return np.array([ vx1,vy1,dvx1,dvy1,
                      vx2,vy2,dvx2,dvy2,
                      vx3,vy3,dvx3,dvy3 ])

def animate(i):    
    """
    This animates the plot generated from a given x and y set of data
    """
    trail_length = 50
    plt.style.use('dark_background')
    
    # Animates the marker portion of the plot
    marker1.set_data([x1[i]],[y1[i]])
    marker2.set_data([x2[i]],[y2[i]])
    marker3.set_data([x3[i]],[y3[i]])
    
    
    # Makes a line that trails the marker, returns markers and trails
    trail1.set_data(x1[max(0, i - trail_length):i+1],y1[max(0, i - trail_length):i+1])
    trail2.set_data(x2[max(0, i - trail_length):i+1],y2[max(0, i - trail_length):i+1])
    trail3.set_data(x3[max(0, i - trail_length):i+1],y3[max(0, i - trail_length):i+1])
    return marker1, marker2, marker3, trail1, trail2, trail3,

def Check(word):
    if word == 'mass' or word == 'Time':
        while True:
            try:
                x = float(input('Please enter a ' +word+ ':'))
                if x < 0:
                    print('Nonphysical value, please try again')
                else:
                    break
            except ValueError:
                print('Invalid input, please enter a valid float')

    else:
        while True:
            try:
                x = float(input('Please enter a ' +word+':'))
                break
            except ValueError:
                print('Invalid input, please enter a valid float ')

    return x

def Figure8():
    """
    Generates a figure 8 via solution from the three body problem
    """
    # Contains the period, the velocities pertaining to the three bodies, mass, and a time frame
    T = 6*6.3490473929 
    v1 = 0.3489048974  
    v2 = 0.5306305100   
    mass = np.array([1,1,1])
    t = np.linspace(0,T,1000)
    
    # Initial Conditions for Body 1
    x01 = -1.0024277970 
    y01 = 0.0041695061
    vx01 = v1
    vy01 = v2
    
    # Initial Conditions for Body 2
    x02 = 1.0024277970
    y02 = -0.0041695061
    vx02 = v1
    vy02 = v2

    # Initial Conditions for Body 3
    x03 = 0
    y03 = 0
    vx03 = -2*v1/mass[2]
    vy03 = -2*v2/mass[2]
    
    # Solution to the Three Body ODE
    Figure8 = solve_ivp(ThreeBody, (0,T), y0=[x01,y01,vx01,vy01,x02,y02,vx02,vy02,x03,y03,vx03,vy03], method='DOP853', t_eval=t)

    x1 = Figure8.y[0]
    y1 = Figure8.y[1]

    x2 = Figure8.y[4]
    y2 = Figure8.y[5]

    x3 = Figure8.y[8]
    y3 = Figure8.y[9]
    
    return x1,y1,x2,y2,x3,y3

def Butterfly():
    """
    Generates a butterfly via solution from the three body problem 
    """
    # Contains the period, the velocities pertaining to the three bodies, mass, and a time frame
    T = 3*6.2399303613 
    v1 = 0.3064392516  
    v2 = 0.1263673939  
    mass = np.array([1,1,1])
    t = np.linspace(0,T,1000)
    
    # Initial Conditions for Body 1
    x01 = -1.0005576155 
    y01 = -0.0029240248 
    vx01 = v1
    vy01 = v2

    # Initial Conditions for Body2
    x02 = 1.0005576155 
    y02 = -0.0029240248 
    vx02 = v1
    vy02 = v2

    # Initial Conditions for Body 3
    x03 = 0
    y03 = 0
    vx03 = -2*v1/mass[2]
    vy03 = -2*v2/mass[2]

    # Solution to the Three Body ODE
    Butterfly = solve_ivp(ThreeBody, (0,T), y0=[x01,y01,vx01,vy01,x02,y02,vx02,vy02,x03,y03,vx03,vy03], method='DOP853', t_eval=t)

    x1 = Butterfly.y[0]
    y1 = Butterfly.y[1]

    x2 = Butterfly.y[4]
    y2 = Butterfly.y[5]

    x3 = Butterfly.y[8]
    y3 = Butterfly.y[9]
    
    return x1,y1,x2,y2,x3,y3

def Bumblebee():
    """
    Generates a bumblebee via solution to the three body problem 
    """
    # Contains the period, the velocities pertaining to the three bodies, mass, and a time frame
    T = 4*64.2532204831 
    v1 = 0.1883232887  
    v2 = 0.5834831526  
    mass = np.array([1,1,1])
    t = np.linspace(0,T,1000)

    # Initial Conditions for Body 1
    x01 = -1.0074958476
    y01 = 0.0081648176
    vx01 = v1
    vy01 = v2

    # Initial Conditions for Body 2
    x02 = 1.0074958476
    y02 = -0.0081648176
    vx02 = v1
    vy02 = v2

    # Initial Conditions for Body 3
    x03 = 0
    y03 = 0
    vx03 = -2*v1/mass[2]
    vy03 = -2*v2/mass[2]

    # Solution to te Three Body ODE
    Bumblebee = solve_ivp(ThreeBody, (0,T), y0=[x01,y01,vx01,vy01,x02,y02,vx02,vy02,x03,y03,vx03,vy03], method='DOP853', t_eval=t)

    x1 = Bumblebee.y[0]
    y1 = Bumblebee.y[1]

    x2 = Bumblebee.y[4]
    y2 = Bumblebee.y[5]

    x3 = Bumblebee.y[8]
    y3 = Bumblebee.y[9]
    
    return x1,y1,x2,y2,x3,y3

def Moth():
    """
    Generates a moth via solution the three body problem 
    """
    # Contains the period, the velocities pertaining to the three bodies, mass, and a time frame
    T = 6*14.8698954200 
    v1 = 0.4646402601  
    v2 = 0.3963456869  
    mass = np.array([1,1,1])
    t = np.linspace(0,T,1000)

    # Initial Conditions of Body 1
    x01 = -0.9989071137
    y01 = -0.0001484864 
    vx01 = v1
    vy01 = v2

    # Initial Conditions of Body 2
    x02 = 0.9989071137
    y02 = 0.0001484864 
    vx02 = v1
    vy02 = v2

    # Initial Conditions of Body 3
    x03 = 0
    y03 = 0
    vx03 = -2*v1/mass[2]
    vy03 = -2*v2/mass[2]

    # Solution to the Three Body ODE
    Moth = solve_ivp(ThreeBody, (0,T), y0=[x01,y01,vx01,vy01,x02,y02,vx02,vy02,x03,y03,vx03,vy03], method='DOP853', t_eval=t)

    x1 = Moth.y[0]
    y1 = Moth.y[1]

    x2 = Moth.y[4]
    y2 = Moth.y[5]

    x3 = Moth.y[8]
    y3 = Moth.y[9]
    
    return x1,y1,x2,y2,x3,y3

def Butterfly_III():
    """ 
    Generates a butterfly via solution of the three body problem
    """
    # Contains the period, the velocities pertaining to the three bodies, mass, and a time frame
    T = 4*20.2667904949 
    v1 = 0.2446132140  
    v2 = 0.3305126876  
    mass = np.array([1,1,1])
    t = np.linspace(0,T,1000)

    # Initial Conditions of Body 1
    x01 = -1.1770534081
    y01 = -0.5225957568
    vx01 = v1
    vy01 = v2

    # Initial Conditions of Body 2
    x02 = 1.1770534081
    y02 = 0.5225957568
    vx02 = v1
    vy02 = v2

    # Initial Conditions of Body 3
    x03 = 0
    y03 = 0
    vx03 = -2*v1/mass[2]
    vy03 = -2*v2/mass[2]

    # Solution to the Three Body ODE
    Butterfly_III = solve_ivp(ThreeBody, (0,T), y0=[x01,y01,vx01,vy01,x02,y02,vx02,vy02,x03,y03,vx03,vy03], method='DOP853', t_eval=t)

    x1 = Butterfly_III.y[0]
    y1 = Butterfly_III.y[1]

    x2 = Butterfly_III.y[4]
    y2 = Butterfly_III.y[5]

    x3 = Butterfly_III.y[8]
    y3 = Butterfly_III.y[9]
    
    return x1,y1,x2,y2,x3,y3

def Moth_II():
    """
    Generates a moth via solution of the three body problem
    """
    # Contains the period, the velocities pertaining to the three bodies, mass, and a time frame
    T = 3*28.6600597106
    v1 = 0.4390528231  
    v2 = 0.4531713698   
    mass = np.array([1,1,1])
    t = np.linspace(0,T,1000)

    # Initial Conditions of Body 1
    x01 = -0.9997857309
    y01 = -0.0003533584
    vx01 = v1
    vy01 = v2

    # Initial Conditions of Body 2
    x02 = 0.9997857309
    y02 = -0.0003533584
    vx02 = v1
    vy02 = v2

    # Initial Conditions of Body 3
    x03 = 0
    y03 = 0
    vx03 = -2*v1/mass[2]
    vy03 = -2*v2/mass[2]

    # Solution to the Three Body ODE
    Moth_II = solve_ivp(ThreeBody, (0,T), y0=[x01,y01,vx01,vy01,x02,y02,vx02,vy02,x03,y03,vx03,vy03], method='DOP853', t_eval=t)

    x1 = Moth_II.y[0]
    y1 = Moth_II.y[1]

    x2 = Moth_II.y[4]
    y2 = Moth_II.y[5]

    x3 = Moth_II.y[8]
    y3 = Moth_II.y[9]
    
    return x1,y1,x2,y2,x3,y3

def Moth_III():
    """
    Generates a moth via solution of the three body problem
    """
    # Contains the period, the velocities pertaining to the three bodies, mass, and a time frame
    T = 3*26.0089024096 
    v1 = 0.3857847594 
    v2 = 0.3732858410   
    mass = np.array([1,1,1])
    t = np.linspace(0,T,1000)

    # Initial Conditions of Body 1
    x01 = -1.0043366457
    y01 = 0.0085104316 
    vx01 = v1
    vy01 = v2

    # Initial Conditions of Body 2
    x02 = 1.0043366457
    y02 = -0.0085104316 
    vx02 = v1
    vy02 = v2

    # Initial Conditions of Body 3
    x03 = 0
    y03 = 0
    vx03 = -2*v1/mass[2]
    vy03 = -2*v2/mass[2]

    # Solution to the Three Body ODE
    Moth_III = solve_ivp(ThreeBody, (0,T), y0=[x01,y01,vx01,vy01,x02,y02,vx02,vy02,x03,y03,vx03,vy03], method='DOP853', t_eval=t)

    x1 = Moth_III.y[0]
    y1 = Moth_III.y[1]

    x2 = Moth_III.y[4]
    y2 = Moth_III.y[5]

    x3 = Moth_III.y[8]
    y3 = Moth_III.y[9]
    
    return x1,y1,x2,y2,x3,y3

def Goggles():
    """
    Generates goggles via solution to the three body problem
    """
    # Contains the period, the velocities pertaining to the three bodies, mass, and a time frame
    T = 2*10.4589089289 
    v1 = 0.0836823874 
    v2 = 0.1276739661   
    mass = np.array([1,1,1])
    t = np.linspace(0,T,1000)

    # Initial Conditions of Body 1
    x01 = -0.9996174046
    y01 = 0.0028671814
    vx01 = v1
    vy01 = v2

    # Initial Conditions of Body 2
    x02 = 0.9996174046
    y02 = -0.0028671814
    vx02 = v1
    vy02 = v2

    # Initial Conditions of Body 3
    x03 = 0
    y03 = 0
    vx03 = -2*v1/mass[2]
    vy03 = -2*v2/mass[2]

    # Solution to the Three Body ODE
    Goggles = solve_ivp(ThreeBody, (0,T), y0=[x01,y01,vx01,vy01,x02,y02,vx02,vy02,x03,y03,vx03,vy03], method='DOP853', t_eval=t)

    x1 = Goggles.y[0]
    y1 = Goggles.y[1]

    x2 = Goggles.y[4]
    y2 = Goggles.y[5]

    x3 = Goggles.y[8]
    y3 = Goggles.y[9]
    
    return x1,y1,x2,y2,x3,y3

def Butterfly_IV():
    """
    Generates a butterfly via solution of the three body problem
    """
    # Contains the period, the velocities pertaining to the three bodies, mass, and a time frame
    T = 79.6790864881 
    v1 = 0.3501405955 
    v2 = 0.0778141747  
    mass = np.array([1,1,1])
    t = np.linspace(0,T,1000)

    # Initial Conditions of Body 1
    x01 = -1.0017029160
    y01 = 0.0041736578 
    vx01 = v1
    vy01 = v2

    # Initial Conditions of Body 2
    x02 = 1.0017029160
    y02 = -0.0041736578 
    vx02 = v1
    vy02 = v2

    # Initial Conditions of Body 3
    x03 = 0
    y03 = 0
    vx03 = -2*v1/mass[2]
    vy03 = -2*v2/mass[2]

    # Solution to the Three Body ODE
    Butterfly_IV = solve_ivp(ThreeBody, (0,T), y0=[x01,y01,vx01,vy01,x02,y02,vx02,vy02,x03,y03,vx03,vy03], method='DOP853', t_eval=t)

    x1 = Butterfly_IV.y[0]
    y1 = Butterfly_IV.y[1]

    x2 = Butterfly_IV.y[4]
    y2 = Butterfly_IV.y[5]

    x3 = Butterfly_IV.y[8]
    y3 = Butterfly_IV.y[9]
    
    return x1,y1,x2,y2,x3,y3

def Dragonfly():
    """
    Generates a dragonfly via solution to the three body problem 
    """
    # Contains the period, the velocities pertaining to the three bodies, mass, and a time frame
    T = 6*20.8385742723
    v1 = 0.0637911125 
    v2 = 0.5950102821  
    mass = np.array([1,1,1])
    t = np.linspace(0,T,1000)

    # Initial Conditions of Body 1
    x01 = -0.9859387540 
    y01 = -0.0288038725
    vx01 = v1
    vy01 = v2

    # Initial Conditions of Body 2
    x02 = 0.9859387540 
    y02 = 0.0288038725
    vx02 = v1
    vy02 = v2

    # Initial Conditions of Body 3
    x03 = 0
    y03 = 0
    vx03 = -2*v1/mass[2]
    vy03 = -2*v2/mass[2]

    # Solution to the Three Body ODE
    Dragonfly = solve_ivp(ThreeBody, (0,T), y0=[x01,y01,vx01,vy01,x02,y02,vx02,vy02,x03,y03,vx03,vy03], method='DOP853', t_eval=t)

    x1 = Dragonfly.y[0]
    y1 = Dragonfly.y[1]

    x2 = Dragonfly.y[4]
    y2 = Dragonfly.y[5]

    x3 = Dragonfly.y[8]
    y3 = Dragonfly.y[9]
    
    return x1,y1,x2,y2,x3,y3

def YinYang():
    """
    Generates a Yin Yang via solution to the three body problem
    """
    # Contains the period, the velocities pertaining to the three bodies, mass, and a time frame
    T = 4*10.6927139709
    v1 = 0.2710001824 
    v2 = 0.3415940623  
    mass = np.array([1,1,1])
    t = np.linspace(0,T,1000)

    # Initial Conditions of Body 1
    x01 = -0.9826146484 
    y01 = -0.0411837391 
    vx01 = v1
    vy01 = v2

    # Initial Conditions of Body 2
    x02 = 0.9826146484 
    y02 = 0.0411837391 
    vx02 = v1
    vy02 = v2

    # Initial Conditions of Body 3
    x03 = 0
    y03 = 0
    vx03 = -2*v1/mass[2]
    vy03 = -2*v2/mass[2]

    # Solution to the Three Body ODE
    YinYang = solve_ivp(ThreeBody, (0,T), y0=[x01,y01,vx01,vy01,x02,y02,vx02,vy02,x03,y03,vx03,vy03], method='DOP853', t_eval=t)

    x1 = YinYang.y[0]
    y1 = YinYang.y[1]

    x2 = YinYang.y[4]
    y2 = YinYang.y[5]

    x3 = YinYang.y[8]
    y3 = YinYang.y[9]
    
    return x1,y1,x2,y2,x3,y3


if __name__ == '__main__':
    
    # checks whether initial or known is selected, then produces solution based on choice(s)
    if args.initial:
        # get initial condition values for first body
        print('You will now have to enter initial x,y,vx, and vy values for Body 1')
        x01 = Check(word='position')
        y01 = Check(word='position')
        vx01 = Check(word='velocity')
        vy01 = Check(word='velocity')
        
        # get initial condition values for second body
        print('You will now have to enter initial x,y,vx, and vy values for Body 2')
        x02 = Check(word='position')
        y02 = Check(word='position')
        vx02 = Check(word='velocity')
        vy02 = Check(word='velocity')
        
        # get initial condition values for third body 
        print('You will now have to enter initial x,y,vx, and vy values for Body 3')
        x03 = Check(word='position')
        y03 = Check(word='position')
        vx03 = Check(word='velocity')
        vy03 = Check(word='velocity')
        
        # get period(time), and the masses of each body
        print('You will now have to enter a period (time), and the mass for each body' )
        T = Check(word='Time')
        m1 = Check(word='mass')
        m2 = Check(word='mass')
        m3 = Check(word='mass')
        
        # turns the masses into array, generates t values based on T
        mass = np.array([m1,m2,m3])
        t = np.linspace(0,T,1000)
        
        # Solves the three body problem based on given initial conditions by user
        Sol = solve_ivp(ThreeBody, (0,T), y0=[x01,y01,vx01,vy01,x02,y02,vx02,vy02,x03,y03,vx03,vy03], method='DOP853', t_eval=t)
        
        x1 = Sol.y[0]
        y1 = Sol.y[1]
        
        x2 = Sol.y[4]
        y2 = Sol.y[5]
        
        x3 = Sol.y[8]
        y3 = Sol.y[9]
        
        # Take the maximum of each x array value
        # this could be a bad idea, test it and see how it works out
        L = 0
        L1 = np.max(x1)
        L2 = np.max(x2)
        L3 = np.max(x3)
        
        # figures out which of those values are largest, so the set the xlim and ylim for plot
        if L1 > L2 and L1 > L3:
            L = L1
        elif L2 > L1 and L2 > L3:
            L = L2
        else:
            L = L3
        
        # sets up the plot
        plt.style.use('dark_background')
        fig, ax = plt.subplots()
        ax.grid()
        
        # makes markers for the plot
        marker1, = plt.plot([],[], '*', lw=4, markersize=10, color='#1f77b4')
        marker2, = plt.plot([],[], '*', lw=4, markersize=10, color='#ff7f0e')
        marker3, = plt.plot([],[], '*', lw=4, markersize=10, color='#2ca02c')
        
        # makes trail behind marker for the plot
        trail1, = plt.plot([],[], '-', color='#aec7e8')
        trail2, = plt.plot([],[], '-', color='#ffbb78')
        trail3, = plt.plot([],[], '-', color='#98df8a')
        
        # set limits on x and y value of graph
        ax.set_xlim(-L,L)
        ax.set_ylim(-L,L)
        
        # animates the solution to the three body problem 
        anim = animation.FuncAnimation(fig,animate,frames=1000,interval=20, blit=True)
        plt.show()
        
    else:
        mass = np.array([1,1,1])

        if args.known == 'Figure8':
            x1,y1,x2,y2,x3,y3 = Figure8()

        elif args.known == 'Bumblebee':
            x1,y1,x2,y2,x3,y3 = Bumblebee()
            
        elif args.known == 'Butterfly':
            x1,y1,x2,y2,x3,y3 = Butterfly()
            
        elif args.known == 'Butterfly III':
            x1,y1,x2,y2,x3,y3 = Butterfly_III()
            
        elif args.known == 'Butterfly IV':
            x1,y1,x2,y2,x3,y3 = Butterfly_IV()
            
        elif args.known == 'Moth I':
            x1,y1,x2,y2,x3,y3 = Moth()
            
        elif args.known == 'Moth II':
            x1,y1,x2,y2,x3,y3 = Moth_II()
            
        elif args.known == 'Moth III':
            x1,y1,x2,y2,x3,y3 = Moth_III()
            
        elif args.known == 'Goggles':
            x1,y1,x2,y2,x3,y3 = Goggles()
            
        elif args.known == 'Dragonfly':
            x1,y1,x2,y2,x3,y3 = Dragonfly()
            
        elif args.known == 'Yin Yang':
            x1,y1,x2,y2,x3,y3 = YinYang()

        else:
            print("Sorry, you didn't enter a known periodic solution, please try again.")

        plt.style.use('dark_background')
        fig, ax = plt.subplots()
        ax.grid()

        # Make marker for the plot
        marker1, = plt.plot([],[], '*', lw=4, markersize=10, color='#1f77b4')
        marker2, = plt.plot([],[], '*', lw=4, markersize=10, color='#ff7f0e')
        marker3, = plt.plot([],[], '*', lw=4, markersize=10, color='#2ca02c')

        # Make trail for the plot
        trail1, = plt.plot([],[], '-', color='#aec7e8')
        trail2, = plt.plot([],[], '-', color='#ffbb78')
        trail3, = plt.plot([],[], '-', color='#98df8a')

        ax.set_xlim(-1.5,1.5)
        ax.set_ylim(-1.5,1.5)
        
        # Animate the function
        anim = animation.FuncAnimation(fig,animate,frames=1000,interval=20,blit=True)
        plt.show()
        