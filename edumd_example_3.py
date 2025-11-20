# Example 1.

# Full interaction

# Simulation starts here

from edumd import *

UNITS = [4, 4, 4]
PBC = [1, 1, 1]

NUNITX=UNITS[0] # Liczba komorek FCC w kierunku x
NUNITY=UNITS[1] # Liczba komorek FCC w kierunku y
NUNITZ=UNITS[2] # Liczba komorek FCC w kierunku z

NATOMS=4*NUNITX*NUNITY*NUNITZ # Liczba atomow

run = 10000
DT = 0.001

RCUT=2.5
SKIN=0.3

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

# grid structure --------------------------------------------------------------

counter, grid_origin, grid_N, grid_dim, NEIGHBOURS_LIST, NEIGHBOURS_LIST_2 = grid_structure(DBOX, INTRANGE)

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

verlet_update = 2


for i in range(run):
        
    
# forces ----------------------------------------------------------------------
        
    #FX, FY, FZ, U = old_forces_lj_full_interaction(RX, RY, RZ, PBC, DBOX, RCUT)
    #FX, FY, FZ, U = forces_lj_full_interaction(RX, RY, RZ, PBC, DBOX, RCUT)
    FX, FY, FZ, U = forces_lj_verlet_list(RX, RY, RZ, v_list, marker, PBC, DBOX, RCUT)
    

# evolution -------------------------------------------------------------------
    
    #RX, RY, RZ, RX_old, RY_old, RZ_old, VX, VY, VZ = newton_verlet(RX, RY, RZ, RX_old, RY_old, RZ_old, FX, FY, FZ, DT) 
    
    #W algorytmie Verleta jest problem z obliczaniem prędkości i EK gdy mamy PBC.
    
    #RX, RY, RZ, VX, VY, VZ = newton_velocity_verlet(RX, RY, RZ, VX, VY, VZ, DT, PBC, DBOX, RCUT)
    #RX, RY, RZ, RX_old, RY_old, RZ_old, VX, VY, VZ = newton_verlet_pbc(RX, RY, RZ, RX_old, RY_old, RZ_old, FX, FY, FZ, PBC, DBOX, DT)    
    RX_old = RX.copy()
    RY_old = RY.copy()
    RZ_old = RZ.copy()
    
    RX, RY, RZ, VXH2, VYH2, VZH2, VX, VY, VZ = newton_verlet_leap_frog_pbc(RX, RY, RZ, VXH2, VYH2, VZH2, FX, FY, FZ, PBC, DBOX, DT, SKIN)
    
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
                
    if verlet_update == 2:
        to_update = list(np.where(DR2>(0.5*SKIN)**2)[0])
        if to_update:
            #print(f"to update: {to_update}")
            where_is_parcicle, cell_content = grid_content(RX, RY, RZ, grid_N, grid_dim, grid_origin)
            for item in to_update:
                particles_to_update, neighbouring_particles = particles_to_update_and_neighbours(item, where_is_parcicle, cell_content, NEIGHBOURS_LIST, NEIGHBOURS_LIST_2)
                v_list, marker = update_verlet_list(RX, RY, RZ, PBC, DBOX, RCUT, SKIN, particles_to_update, neighbouring_particles, v_list, marker)
                #print("UPDATE")
                for jtem in particles_to_update:
                    DRX[jtem] = 0
                    DRY[jtem] = 0
                    DRZ[jtem] = 0
                    DR2[jtem] = 0
        
        
        #for item in ind:
        #    particles_to_update, neighbouring_particles = particles_to_update_and_neighbours(item, where_is_parcicle, cell_content, NEIGHBOURS_LIST, NEIGHBOURS_LIST_2)
            
            #v_list_2, marker_2 = update_verlet_list(RX, RY, RZ, PBC, DBOX, RCUT, SKIN, particles_to_update, neighbouring_particles, v_list, marker)
           

        
    #RX, RY, RZ, VX, VY, VZ, U = newton_velocity_verlet_pbc(RX, RY, RZ, VX, VY, VZ, PBC, DBOX, RCUT, DT)
    
    
    
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

