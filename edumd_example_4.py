from edumd import *

UNITS = [4, 4, 4]
PBC = [1, 1, 1]

NUNITX=UNITS[0] # Liczba komorek FCC w kierunku x
NUNITY=UNITS[1] # Liczba komorek FCC w kierunku y
NUNITZ=UNITS[2] # Liczba komorek FCC w kierunku z

NATOMS=4*NUNITX*NUNITY*NUNITZ # Liczba atomow

run = 1000
DT = 0.001

RCUT=2.5
SKIN=0.5

mass=1
KB=1
g=3*NATOMS-3
rho=0.8
T0=2
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

verlet_update = 1

for i in range(run):
        
    
    RX_old = RX.copy()
    RY_old = RY.copy()
    RZ_old = RZ.copy()
    
    RX, RY, RZ, VXH2, VYH2, VZH2, VX, VY, VZ = newton_verlet_leap_frog_pbc(RX, RY, RZ, VXH2, VYH2, VZH2, FX, FY, FZ, PBC, DBOX, DT, SKIN)
    
    DRX = DRX+(RX_old-RX)
    DRY = DRY+(RY_old-RY)
    DRZ = DRZ+(RZ_old-RZ)
    DR2 = DR2 + DRX*DRX + DRY*DRY + DRZ*DRZ
                
    if verlet_update == 1:
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
        
            
    
    if i%10 == 0:
        print(f"i = {i} out of {run}, tot_en = {tot_en}")
        
        
    tot_en = total_energy_per_atom(U, VX, VY, VZ)
    tot_en_list.append(tot_en)
    EP.append(potential_energy(U))
    K.append(kinetic_energy(VX, VY, VZ))
    sim_time.append(i*DT)


t2 = time.time()

print(f"time: {t2-t1}")

