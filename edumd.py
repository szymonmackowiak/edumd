# last modified 01.10.2024.

import matplotlib.pyplot as plt
import numpy as np
from itertools import product

import time

t1 = time.time()

#------------------------------------------------------------------------------

def simple_plot_1(x, y1, title, xlabel, ylabel, legend_y1, xlim, ylim):
    
    """example: simple_plot_1(sim_time, tot_en_list, "U/N", "t", "U", "U1", [0, 2], [-3, -2])"""
    
    plt.plot(x, y1, '-b', linewidth = 1)
    plt.ylabel(ylabel, fontsize = '14')
    plt.xlabel(xlabel, fontsize = '14')
    plt.title(f"{title}", fontsize = '16')
    plt.ylim(ylim[0], ylim[1])
    plt.xlim(xlim[0], xlim[1])
    plt.grid(color='grey', linestyle='--', linewidth=0.2)
    plt.legend(legend_y1)
    plt.show()

#------------------------------------------------------------------------------

def simple_plot_2(x, y1, y2, title, xlabel, ylabel, legend_y1, legend_y2, xlim, ylim):
    
    """example: simple_plot_2(sim_time, tot_en_list, tot_en_list, "U/N", "t", "U", "U1", "U2", [0, 2], [-3, -2])"""
    
    plt.plot(x, y1, '-b', x, y2, '-r', linewidth = 1)
    plt.ylabel(ylabel, fontsize = '14')
    plt.xlabel(xlabel, fontsize = '14')
    plt.title(f"{title}", fontsize = '16')
    plt.ylim(ylim[0], ylim[1])
    plt.xlim(xlim[0], xlim[1])
    plt.grid(color='grey', linestyle='--', linewidth=0.2)
    plt.legend([legend_y1, legend_y2])
    plt.show()

#------------------------------------------------------------------------------

def display(x_list, y_list, z_list):
    
    NATOMS = len(x_list)
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for i in range(NATOMS):
        x = x_list[i]
        y = y_list[i]
        z = z_list[i]
        ax.scatter(x, y, z, color='b', s=100)

    plt.show()
    
#------------------------------------------------------------------------------
    
def center_of_mass(RX, RY, RZ):
    
    NATOMS = len(RX)
    
    CX=0 
    CY=0 
    CZ=0

    for i in range(NATOMS): 
    
        CX=CX+RX[i]
        CY=CY+RY[i]
        CZ=CZ+RZ[i]
        
    CX=CX/NATOMS 
    CY=CY/NATOMS 
    CZ=CZ/NATOMS
    
    return CX, CY, CZ

#------------------------------------------------------------------------------

def center_system(RX, RY, RZ):
    
    NATOMS = len(RX)
    CX, CY, CZ = center_of_mass(RX, RY, RZ)
    
    for i in range(NATOMS):
        RX[i]=RX[i]-CX
        RY[i]=RY[i]-CY
        RZ[i]=RZ[i]-CZ
    
    return RX, RY, RZ
    
#------------------------------------------------------------------------------

def total_momentum(VX, VY, VZ):
    
    NATOMS = len(VX)
    
    total_px = 0
    total_py = 0
    total_pz = 0
          
    for i in range(NATOMS):
        
        total_px = total_px + VX[i]
        total_py = total_py + VY[i]
        total_pz = total_pz + VZ[i]
        
    return total_px, total_py, total_pz

#------------------------------------------------------------------------------

def zero_momentum(VX, VY, VZ):
    
    NATOMS = len(VX)
    total_px, total_py, total_pz = total_momentum(VX, VY, VZ)
    
    for i in range(NATOMS):
        
        VX[i]=VX[i]-total_px/NATOMS
        VY[i]=VY[i]-total_py/NATOMS
        VZ[i]=VZ[i]-total_pz/NATOMS
        
    return VX, VY, VZ

#------------------------------------------------------------------------------
def apply_pcb(x, pcb, L):
        
    return x-pcb*L*round(x/L)

#------------------------------------------------------------------------------
def apply_pcb_to_particle(p_RX, p_RY, p_RZ, PBC, DBOX):
    
    XPERIOD = PBC[0]
    YPERIOD = PBC[1]
    ZPERIOD = PBC[2]
    
    p_RX=p_RX-XPERIOD*(DBOX[0]*round(p_RX/DBOX[0]));
    p_RY=p_RY-YPERIOD*(DBOX[1]*round(p_RY/DBOX[1]));
    p_RZ=p_RZ-ZPERIOD*(DBOX[2]*round(p_RZ/DBOX[2]));
        
    return p_RX, p_RY, p_RZ

#------------------------------------------------------------------------------

def apply_pcb_to_system(RX, RY, RZ, PBC, DBOX):
    
    XPERIOD = PBC[0]
    YPERIOD = PBC[1]
    ZPERIOD = PBC[2]
    
    NATOMS = len(RX)
    
    for i in range(NATOMS):
        RX[i]=RX[i]-XPERIOD*(DBOX[0]*round(RX[i]/DBOX[0]));
        RY[i]=RY[i]-YPERIOD*(DBOX[1]*round(RY[i]/DBOX[1]));
        RZ[i]=RZ[i]-ZPERIOD*(DBOX[2]*round(RZ[i]/DBOX[2]));
    return RX, RY, RZ

#------------------------------------------------------------------------------

def block_fcc(units, d_box, rho, disorder_magnitude):
    
    np.random.seed(0)
    
    n_unit_x, n_unit_y, n_unit_z = units[0], units[1], units[2]
    
    b_x, b_y, b_z = d_box[0], d_box[1], d_box[2] 

    n_atoms = 4 * n_unit_x * n_unit_y * n_unit_z

    dist_x, dist_y, dist_z = b_x/n_unit_x, b_y/n_unit_y, b_z/n_unit_z

    RX = np.zeros(n_atoms)
    RY = np.zeros(n_atoms)
    RZ = np.zeros(n_atoms)

    RX[0], RY[0], RZ[0] = 0, 0, 0
    RX[1], RY[1], RZ[1] = 0, dist_y/2, dist_z/2
    RX[2], RY[2], RZ[2] = dist_x/2, 0, dist_z/2
    RX[3], RY[3], RZ[3] = dist_x/2, dist_y/2, 0

    m = 0

    for i in range(n_unit_z):
        for j in range(n_unit_y): 
            for k in range(n_unit_x):
                for l in range(4):
                    random_number = np.random.normal(0, 1)
                    RX[l+m] = RX[l] + dist_x * k + disorder_magnitude * random_number
                    
                    random_number = np.random.normal(0, 1)
                    RY[l+m] = RY[l] + dist_y * j + disorder_magnitude * random_number
                    
                    random_number = np.random.normal(0, 1)
                    RZ[l+m] = RZ[l] + dist_z * i + disorder_magnitude * random_number
                m = m + 4

    RX, RY, RZ = center_system(RX, RY, RZ)

    return RX, RY, RZ

#------------------------------------------------------------------------------

def shuffle_initial_speeds(NATOMS, T0):
    
    np.random.seed(0)
    
    VX = np.zeros(NATOMS)
    VY = np.zeros(NATOMS)
    VZ = np.zeros(NATOMS)
    
    for i in range(NATOMS):
        
        random_number = np.random.normal(0, 1)
        VX[i]=random_number * T0
    
        random_number = np.random.normal(0, 1)
        VY[i]=random_number * T0
    
        random_number = np.random.normal(0, 1)
        VZ[i]=random_number * T0
        
    VX, VY, VZ = zero_momentum(VX, VY, VZ)
    
    return VX, VY, VZ

#------------------------------------------------------------------------------

def verlet_list(RX, RY, RZ, PBC, DBOX, RCUT, SKIN):
    
    XPERIOD = PBC[0]
    YPERIOD = PBC[1]
    ZPERIOD = PBC[2]
    
    BX = DBOX[0]
    BY = DBOX[1]
    BZ = DBOX[2] 
    
    NATOMS = len(RX)
    
    v_list = []
    
    marker = []
    

    for i in range(NATOMS):
    
        v_row = []
        
        for j in range(NATOMS):
        
            X=0 
            Y=0 
            Z=0 
            R2=0 
            X=RX[i]-RX[j] 
            Y=RY[i]-RY[j] 
            Z=RZ[i]-RZ[j]
            
            X=X-XPERIOD*(BX*round(X/BX))
            Y=Y-YPERIOD*(BY*round(Y/BY))
            Z=Z-ZPERIOD*(BZ*round(Z/BZ))
            
            R2=X*X+Y*Y+Z*Z
                            
            if ((R2<((RCUT+SKIN)*(RCUT+SKIN))) and (i!=j)):
                v_row.append(j)
                
        v_list.append(v_row)
        marker.append(len(v_row))
    
    return v_list, marker

#------------------------------------------------------------------------------

def lj_energy_full_interaction_no_pbc(RX, RY, RZ):
    
    NATOMS = len(RX)
    
    a = 1
    b = 1
    N1 = 12
    N2 = 6
      
    U = np.zeros(NATOMS)
              
    for i in range(0, NATOMS):
        
        EUP=0.0
    
        for j in range(i+1, NATOMS):
            
            X=0 
            Y=0 
            Z=0 
            R2=0 
            R=0 
            RI=0
            X=RX[i]-RX[j] 
            Y=RY[i]-RY[j] 
            Z=RZ[i]-RZ[j]
              
            R2=X*X+Y*Y+Z*Z
            R=R2**0.5
                
            if R>0:
                RI=1/R
            else:
                print(f"WARNING: R={RI}, i={i}, j={j}")
                                  
            EUP=4*a*((RI**N1)-b*(RI**N2));

            U[i]=U[i]+0.5*EUP
            U[j]=U[j]+0.5*EUP

    return U

#------------------------------------------------------------------------------

def forces_lj_verlet_list(RX, RY, RZ, v_list, marker, PBC, DBOX, RCUT):
    
    XPERIOD = PBC[0]
    YPERIOD = PBC[1]
    ZPERIOD = PBC[2]
    
    BX = DBOX[0]
    BY = DBOX[1]
    BZ = DBOX[2] 
    
    NATOMS = len(RX)
    
    a = 1
    b = 1

    N1=12
    N2=6 

    shift=-4*(RCUT**(-N1)-RCUT**(-N2))

    FX = np.zeros(NATOMS)
    FY = np.zeros(NATOMS)
    FZ = np.zeros(NATOMS)
    U = np.zeros(NATOMS)

    for i in range(NATOMS):
        licznikomp=0
        
        for k in range(marker[i]):
            licznikomp=licznikomp+1
            j=v_list[i][k]
            if i>j:
                EUP=0.0
                FXP=0.0
                FYP=0.0
                FZP=0.0
                X=0
                Y=0
                Z=0
                R2=0
                R=0
                RI=0
                
                X=RX[i]-RX[j]
                Y=RY[i]-RY[j]
                Z=RZ[i]-RZ[j]
                
                X=X-XPERIOD*(BX*round(X/BX))
                Y=Y-YPERIOD*(BY*round(Y/BY))
                Z=Z-ZPERIOD*(BZ*round(Z/BZ))
                
                R2=X*X+Y*Y+Z*Z
                R=R2**0.5
                
                if R>0:
                    RI=1/R
                else:
                    print(f"WARNING: R={RI}")
                
                if R2<RCUT*RCUT:
                    
                    EUP=4*a*((RI**N1)-b*(RI**N2))
                    EUP=EUP+shift
                        
                    FXP=4*a*(N1*X*RI**(N1+2)-b*N2*X*RI**(N2+2))
                    FYP=4*a*(N1*Y*RI**(N1+2)-b*N2*Y*RI**(N2+2))
                    FZP=4*a*(N1*Z*RI**(N1+2)-b*N2*Z*RI**(N2+2))

                elif R2>=RCUT*RCUT:
                    EUP=0.0
                    FXP=0.0 
                    FYP=0.0 
                    FZP=0.0
                
                FX[i]=FX[i]+FXP
                FY[i]=FY[i]+FYP
                FZ[i]=FZ[i]+FZP 
                
                U[i]= U[i] + 0.5*EUP
  
                FX[j]=FX[j]-FXP 
                FY[j]=FY[j]-FYP
                FZ[j]=FZ[j]-FZP 
                
                U[j]= U[j] + 0.5*EUP
    
    return FX, FY, FZ, U


#------------------------------------------------------------------------------

def forces_lj_full_interaction(RX, RY, RZ, PBC, DBOX, RCUT):
    
    XPERIOD = PBC[0]
    YPERIOD = PBC[1]
    ZPERIOD = PBC[2]
    
    BX = DBOX[0]
    BY = DBOX[1]
    BZ = DBOX[2] 
    
    NATOMS = len(RX)
    
    a = 1
    b = 1
    N1 = 12
    N2 = 6
    
    shift=-4*(RCUT**(-N1)-RCUT**(-N2))
      
    FX = np.zeros(NATOMS)
    FY = np.zeros(NATOMS)
    FZ = np.zeros(NATOMS)
    U = np.zeros(NATOMS)
              
    for i in range(0, NATOMS):
        
        EUP=0.0
        FXP=0.0
        FYP=0.0
        FZP=0.0
    
        for j in range(i+1, NATOMS):
            
            #print(f"i={i}, j={j}")
              
            X=0 
            Y=0 
            Z=0 
            R2=0 
            R=0 
            RI=0
            X=RX[i]-RX[j] 
            Y=RY[i]-RY[j] 
            Z=RZ[i]-RZ[j]
    
            X=X-XPERIOD*(BX*round(X/BX))
            Y=Y-YPERIOD*(BY*round(Y/BY))
            Z=Z-ZPERIOD*(BZ*round(Z/BZ))
              
            R2=X*X+Y*Y+Z*Z
            R=R2**0.5
                
            if R>0:
                RI=1/R
            else:
                print(f"WARNING: R={RI}, i={i}, j={j}")
              
            if R2<RCUT*RCUT:
                    
                EUP=4*a*((RI**N1)-b*(RI**N2));
                EUP=EUP+shift;

                FXP=4*a*(N1*X*RI**(N1+2)-b*N2*X*RI**(N2+2));
                FYP=4*a*(N1*Y*RI**(N1+2)-b*N2*Y*RI**(N2+2));
                FZP=4*a*(N1*Z*RI**(N1+2)-b*N2*Z*RI**(N2+2));
              
            elif R2>=RCUT*RCUT:
                    
                EUP=0.0
                FXP=0.0
                FYP=0.0
                FZP=0.0

            FX[i]=FX[i]+FXP
            FY[i]=FY[i]+FYP
            FZ[i]=FZ[i]+FZP
            U[i]=U[i]+0.5*EUP
            
            FX[j]=FX[j]-FXP
            FY[j]=FY[j]-FYP
            FZ[j]=FZ[j]-FZP
            U[j]=U[j]+0.5*EUP

    return FX, FY, FZ, U
              
#------------------------------------------------------------------------------

def kinetic_energy(VX, VY, VZ):
    
    NATOMS = len(VX)
    EKC = 0
    for i in range(NATOMS):
        EKC=EKC+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i])
    
    return EKC

#------------------------------------------------------------------------------

def potential_energy(U):
    
    NATOMS = len(U)
    EPC = 0
    for i in range(NATOMS):
        EPC=EPC+U[i]
    
    return EPC

#------------------------------------------------------------------------------
    
def total_energy_per_atom(U, VX, VY, VZ):
    
    NATOMS = len(U)
    
    EKC = kinetic_energy(VX, VY, VZ)
    EPC = potential_energy(U)
    
    tot_en = (EPC+EKC)/NATOMS
    
    return tot_en
    
#------------------------------------------------------------------------------

def newton_rk1(RX, RY, RZ, VX, VY, VZ, FX, FY, FZ, DT):
    
    NATOMS = len(RX)
    
    for i in range(NATOMS):
        VX[i] = VX[i] + FX[i] * DT
        VY[i] = VY[i] + FY[i] * DT
        VZ[i] = VZ[i] + FZ[i] * DT
        
        RX[i] = RX[i] + VX[i] * DT + 0.5 * FX[i] * DT * DT
        RY[i] = RY[i] + VY[i] * DT + 0.5 * FY[i] * DT * DT
        RZ[i] = RZ[i] + VZ[i] * DT + 0.5 * FZ[i] * DT * DT
    
    
    return RX, RY, RZ, VX, VY, VZ

#------------------------------------------------------------------------------

def newton_verlet_initial_step(RX, RY, RZ, VX, VY, VZ, FX, FY, FZ, DT):
    
    RX_old = np.copy(RX)
    RY_old = np.copy(RY)
    RZ_old = np.copy(RZ)
    
    RX, RY, RZ, VX, VY, VZ = newton_rk1(RX, RY, RZ, VX, VY, VZ, FX, FY, FZ, DT)
    
    return RX, RY, RZ, RX_old, RY_old, RZ_old, VX, VY, VZ
    
#------------------------------------------------------------------------------

def newton_verlet_pbc(RX, RY, RZ, RX_old, RY_old, RZ_old, FX, FY, FZ, PBC, DBOX, DT):
    
    NATOMS = len(RX)
    
    XPERIOD = PBC[0]
    YPERIOD = PBC[1]
    ZPERIOD = PBC[2]
    
    BX = DBOX[0]
    BY = DBOX[1]
    BZ = DBOX[2] 
    
    RX_new = np.zeros(NATOMS)
    RY_new = np.zeros(NATOMS)
    RZ_new = np.zeros(NATOMS)
    
    VX = np.zeros(NATOMS)
    VY = np.zeros(NATOMS)
    VZ = np.zeros(NATOMS)
    
    for i in range(NATOMS):
        
        RX_new[i] = 2*RX[i]-RX_old[i]+(DT**2)*FX[i]
        RX_new[i] = apply_pcb(RX_new[i], XPERIOD, BX)
        
        RY_new[i] = 2*RY[i]-RY_old[i]+(DT**2)*FY[i]
        RY_new[i] = apply_pcb(RY_new[i], YPERIOD, BY)
        
        RZ_new[i] = 2*RZ[i]-RZ_old[i]+(DT**2)*FZ[i]
        RZ_new[i] = apply_pcb(RZ_new[i], ZPERIOD, BZ)
        
        DRX = apply_pcb(RX_new[i]-RX_old[i], XPERIOD, BX)
        DRY = apply_pcb(RY_new[i]-RY_old[i], YPERIOD, BY)
        DRZ = apply_pcb(RZ_new[i]-RZ_old[i], ZPERIOD, BZ)
        
        VX[i] = (DRX)/(2*DT)
        VY[i] = (DRY)/(2*DT)
        VZ[i] = (DRZ)/(2*DT)
        
    return RX_new, RY_new, RZ_new, RX, RY, RZ, VX, VY, VZ

#------------------------------------------------------------------------------

def newton_half_step_leap_frog(VX, VY, VZ, FX, FY, FZ, DT):
    
    NATOMS = len(VX)
    
    VXH1 = np.zeros(NATOMS)
    VYH1 = np.zeros(NATOMS)
    VZH1 = np.zeros(NATOMS)
    
    for i in range(NATOMS):
        VXH1[i]=VX[i]+0.5*DT*FX[i]
        VYH1[i]=VY[i]+0.5*DT*FY[i] 
        VZH1[i]=VZ[i]+0.5*DT*FZ[i]

    VXH2=np.copy(VXH1)
    VYH2=np.copy(VYH1)
    VZH2=np.copy(VZH1)
    
    return VXH2, VYH2, VZH2

#------------------------------------------------------------------------------

def newton_verlet_leap_frog_pbc(RX, RY, RZ, VXH2, VYH2, VZH2, FX, FY, FZ, PBC, DBOX, DT):

    NATOMS = len(RX)
    
    XPERIOD = PBC[0]
    YPERIOD = PBC[1]
    ZPERIOD = PBC[2]
    
    BX = DBOX[0]
    BY = DBOX[1]
    BZ = DBOX[2] 
    
    VXH1 = np.zeros(NATOMS)
    VYH1 = np.zeros(NATOMS)
    VZH1 = np.zeros(NATOMS)
    
    VX = np.zeros(NATOMS)
    VY = np.zeros(NATOMS)
    VZ = np.zeros(NATOMS)
    
    RX_old = np.zeros(NATOMS)
    RY_old = np.zeros(NATOMS)
    RZ_old = np.zeros(NATOMS)
    
    for i in range(NATOMS):
        
        VXH1[i]=VXH2[i] 
        VYH1[i]=VYH2[i] 
        VZH1[i]=VZH2[i]

        RX_old[i]=RX[i] 
        RY_old[i]=RY[i] 
        RZ_old[i]=RZ[i]

        VXH2[i]=VXH1[i]+DT*(FX[i])
        VYH2[i]=VYH1[i]+DT*(FY[i])
        VZH2[i]=VZH1[i]+DT*(FZ[i])    
        
        RX[i]=RX_old[i]+DT*VXH2[i]
        RX[i]=apply_pcb(RX[i], XPERIOD, BX)

        RY[i]=RY_old[i]+DT*VYH2[i]
        RY[i] = apply_pcb(RY[i], YPERIOD, BY)
        
        RZ[i]=RZ_old[i]+DT*VZH2[i]
        RZ[i] = apply_pcb(RZ[i], ZPERIOD, BZ)
        
        VX[i]=0.5*(VXH1[i]+VXH2[i])
        VY[i]=0.5*(VYH1[i]+VYH2[i])
        VZ[i]=0.5*(VZH1[i]+VZH2[i])
        
    return RX, RY, RZ, VXH2, VYH2, VZH2, VX, VY, VZ

#------------------------------------------------------------------------------

def newton_velocity_verlet_pbc(RX, RY, RZ, VX, VY, VZ, PBC, DBOX, RCUT, DT):

    NATOMS = len(RX)
    
    XPERIOD = PBC[0]
    YPERIOD = PBC[1]
    ZPERIOD = PBC[2]
    
    BX = DBOX[0]
    BY = DBOX[1]
    BZ = DBOX[2] 
    
    FX, FY, FZ, U = forces_lj_full_interaction(RX, RY, RZ, PBC, DBOX, RCUT)
    
    FX_old = FX.copy()
    FY_old = FY.copy()
    FZ_old = FZ.copy()

    
    for i in range(NATOMS):

        RX[i]=RX[i] + VX[i]*DT + 0.5*(FX[i])*DT*DT
        RX[i]=apply_pcb(RX[i], XPERIOD, BX)
        
        RY[i]=RY[i] + VY[i]*DT + 0.5*(FY[i])*DT*DT
        RY[i] = apply_pcb(RY[i], YPERIOD, BY)
        
        RZ[i]=RZ[i] + VZ[i]*DT + 0.5*(FZ[i])*DT*DT
        RZ[i] = apply_pcb(RZ[i], ZPERIOD, BZ)
        
    FX, FY, FZ, U = forces_lj_full_interaction(RX, RY, RZ, PBC, DBOX, RCUT)
    
    for i in range(NATOMS):
        
        VX[i] = VX[i] + 0.5*(FX[i] + FX_old[i])*DT
        VY[i] = VY[i] + 0.5*(FY[i] + FY_old[i])*DT
        VZ[i] = VZ[i] + 0.5*(FZ[i] + FZ_old[i])*DT

        
    return RX, RY, RZ, VX, VY, VZ, U

#------------------------------------------------------------------------------



    

    
