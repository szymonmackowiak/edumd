# Example 2.

# Intercation: Verlet list
# with PBC
# Algorithm Velocity Verlet (forces computed inside the integration function)

# Simulation starts here

from edumd import *

UNITS = [3, 3, 3]
PBC = [1, 1, 1]

NUNITX=UNITS[0] # Number of FCC cells in the x direction
NUNITY=UNITS[1] # Number of FCC cells in the y direction
NUNITZ=UNITS[2] # Number of FCC cells in the z direction

NATOMS=4*NUNITX*NUNITY*NUNITZ # Number of atoms
run = 5000
DT = 0.001

RCUT=2.5
SKIN=0.5

mass=1
KB=1
g=3*NATOMS-3
rho=0.2
T0=2
eq_dev=0.05

BX=(NATOMS/((NUNITZ/NUNITX)*rho))**(1/3)
BY=BX
BZ=(NUNITZ/NUNITX)*BX 

DBOX = [BX, BY, BZ]

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

FX, FY, FZ, U = forces_lj_full_interaction(RX, RY, RZ, PBC, DBOX, RCUT)

VXH2, VYH2, VZH2 = newton_half_step_leap_frog(VX, VY, VZ, FX, FY, FZ, DT)

RX, RY, RZ, RX_old, RY_old, RZ_old, VX, VY, VZ = newton_verlet_initial_step(RX, RY, RZ, VX, VY, VZ, FX, FY, FZ, DT)
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


for i in range(run):
        
    
# forces ----------------------------------------------------------------------
        
    #FX, FY, FZ, U = old_forces_lj_full_interaction(RX, RY, RZ, PBC, DBOX, RCUT)
    # Forcces computed inside the integration function

# evolution -------------------------------------------------------------------
    
    RX, RY, RZ, VX, VY, VZ, U = newton_velocity_verlet_pbc(RX, RY, RZ, VX, VY, VZ, PBC, DBOX, RCUT, DT)
        
    if i%20 == 0:
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

