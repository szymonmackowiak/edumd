# Example 3.

# Intercation: Verlet list
# with PBC
# Algorithm Verlet Leap-Frog (forces computed outside the integration function)

# Simulation starts here

from edumd import *

UNITS = [4, 4, 4]
PBC = [1, 1, 1]

NUNITX=UNITS[0] # Number of FCC cells in the x direction
NUNITY=UNITS[1] # Number of FCC cells in the y direction
NUNITZ=UNITS[2] # Number of FCC cells in the z direction

NATOMS=4*NUNITX*NUNITY*NUNITZ # Number of atoms

run = 10000
DT = 0.001

RCUT=2.5
SKIN=0.5

mass=1
KB=1
g=3*NATOMS-3
rho=0.8442
T0=1
eq_dev=0.05

BX=(NATOMS/((NUNITZ/NUNITX)*rho))**(1/3)
BY=BX
BZ=(NUNITZ/NUNITX)*BX 

DBOX = [BX, BY, BZ]
INTRANGE = RCUT + SKIN


#------------------------------------------------------------------------------    


DRX = np.zeros(NATOMS)
DRY = np.zeros(NATOMS)
DRZ = np.zeros(NATOMS)
DR2 = np.zeros(NATOMS)


tot_en_list = []

EP = []
K = []


sim_time = []

RX, RY, RZ = block_fcc(UNITS, DBOX, rho, eq_dev)
VX, VY, VZ = shuffle_initial_speeds(NATOMS, T0)

v_list, marker = verlet_list(RX, RY, RZ, PBC, DBOX, RCUT, SKIN)

FX, FY, FZ, U = forces_lj_verlet_list(RX, RY, RZ, v_list, marker, PBC, DBOX, RCUT)

VXH2, VYH2, VZH2 = newton_half_step_leap_frog(VX, VY, VZ, FX, FY, FZ, DT)

tot_en = total_energy_per_atom(U, VX, VY, VZ)

save_RX = open(f"output/RX.txt", "w")
save_RY = open(f"output/RY.txt", "w")
save_RZ = open(f"output/RZ.txt", "w")
save_RX.close()
save_RY.close()
save_RZ.close()


save_RX = open(f"output/RX.txt", "a")
save_RY = open(f"output/RY.txt", "a")
save_RZ = open(f"output/RZ.txt", "a")

verlet_update = 1


for i in range(run):
        
    
# forces ----------------------------------------------------------------------
        
    FX, FY, FZ, U = forces_lj_verlet_list(RX, RY, RZ, v_list, marker, PBC, DBOX, RCUT)
    
# evolution -------------------------------------------------------------------
     
    RX_old = RX.copy()
    RY_old = RY.copy()
    RZ_old = RZ.copy()
    
    RX, RY, RZ, VXH2, VYH2, VZH2, VX, VY, VZ = newton_verlet_leap_frog_pbc(RX, RY, RZ, VXH2, VYH2, VZH2, FX, FY, FZ, PBC, DBOX, DT)
    
    DRX = DRX+(RX_old-RX)
    DRY = DRY+(RY_old-RY)
    DRZ = DRZ+(RZ_old-RZ)
    DR2 = DR2 + DRX*DRX + DRY*DRY + DRZ*DRZ
    
    
    if verlet_update == 1:
        for j in range(len(DR2)):
            if DR2[j] > (0.5*SKIN)**2:
                v_list, marker = verlet_list(RX, RY, RZ, PBC, DBOX, RCUT, SKIN)
                print("UPDATE")
                DRX = np.zeros(NATOMS)
                DRY = np.zeros(NATOMS)
                DRZ = np.zeros(NATOMS)
                DR2 = np.zeros(NATOMS)
    
    
    if i%10 == 0:
        print(f"i = {i} out of {run}, tot_en = {tot_en}")
        for i in range(NATOMS):
            save_RX.write(f"{str(round(RX[i], 4))} ")
            save_RY.write(f"{str(round(RY[i], 4))} ")
            save_RZ.write(f"{str(round(RZ[i], 4))} ")
        save_RX.write("\n")
        save_RY.write("\n")
        save_RZ.write("\n")
        
    tot_en = total_energy_per_atom(U, VX, VY, VZ)
    tot_en_list.append(tot_en)
    EP.append(potential_energy(U))
    K.append(kinetic_energy(VX, VY, VZ))
    sim_time.append(i*DT)

save_RX.close()
save_RY.close()
save_RZ.close()

plt.plot(tot_en_list)
plt.show()


save_tot_en = open(f"output/tot_en.txt", "w")
save_tot_en.write(str(tot_en_list))
save_tot_en.close()

t2 = time.time()

print(f"time: {t2-t1}")

