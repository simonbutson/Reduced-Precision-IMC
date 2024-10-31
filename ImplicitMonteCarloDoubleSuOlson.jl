# 1D Implicit Monte Carlo - Double Precision Su Olson Problem
# Simon Butson
# Adapted from FCIMC code written in Python

using Random
using Plots
using DelimitedFiles

# ========================
# Input Deck
# ========================
"""Read input file for simulation dependent parameters"""
inputs = readdlm(raw"calc1SO.in") # Change path depending on location of input files

calcnum = String(inputs[1,2])[5:end] # Calculation number - take characters in the 5th place till end of string "calcx"
dt = Float64(inputs[2,2]) # Timestep length
width = Float64(inputs[3,2]) # Slab width
dx = Float64(inputs[4,2]) # Spatial cell size
timesteps = Float64(inputs[5,2]) # Number of timesteps
n_input = Int64(inputs[6,2]) # Starting particles input
sigma_s = Float64(inputs[7,2]) # Scatter cross-section
sigma_a = Float64(inputs[8,2]) # Absorption cross-section
x_0 = Float64(inputs[9,2]) # Source width
tau_0 = Float64(inputs[10,2]) # Source active time
epsilon = Float64(inputs[11,2]) # Epsilon parameter
T_init = Float64(inputs[12,2]) # Initial slab temperature

tmax = timesteps * dt # Maximum time
shakes = dt/(1e-8) # Timesteps in shakes (1e-8)

Ncells = Int(width/dx) # Number of cells

sourceCells = Int(x_0/dx) # Number of cells with source

particleVars = 8 # Particle array variables
 
n_census = 0 # Number of particles stored in census
n_max = 500000000 # Reasonable upper limit - can change as needed

# ========================
# Mesh Arrays & Quantities
# ========================

# Particle array structure
# (time, cell-index, x-position, mu, frequency, energy, startenergy)
#particles = zeros(Nparticles,particleVars)

mesh_centers = LinRange(dx/2, width - dx/2, Ncells)
mesh_nodes = LinRange(0,width,Ncells+1)

particles = Vector{Vector{Float64}}()
mesh_tempsaved = Vector{Vector{Float64}}()
radenergy_saved = Vector{Vector{Float64}}()
materialenergy_saved = Vector{Vector{Float64}}()

mesh_tempinit = T_init .* ones(Ncells) # Initial temperature (keV)
mesh_temp = copy(mesh_tempinit) # Temperature (keV)
mesh_fleck = zeros(Ncells) # Fleck factor
mesh_beta = zeros(Ncells) # Beta factor
mesh_sigma = zeros(Ncells) # Planck mean opacity
mesh_energydep = zeros(1,Ncells) # Deposited energy in a timestep
material_energy = zeros(Ncells) # Material energy density
radenergydens = zeros(Ncells) # Radiation energy density

source_energies = Vector{Float64}()
energy_increases = Vector{Float64}()
lost_energy = Vector{Float64}()
stored_energy = Vector{Float64}()

PLE = zeros(Ncells)


# Physical constants - Scaled Units

#phys_h = 4.1356675E-18 # Planck constant [keV*s]

#phys_c = 3 * 10^10 # Speed of light [cm/s]
phys_c = 1.0

#phys_sb = 2.0 * pi^5 / (15.0 * (phys_h^3) * (phys_c^2)) # Stefan-Boltzmann constant [kev^-3*cm^-2*s^-1]
#phys_sb = 6.4092726315814575E+32 # Stefan-Boltzmann constant [kev^-3*cm^-2*s^-1]
phys_sb = 1

#phys_a = 4.0 * phys_sb / phys_c # Radiation constant [keV^-3*cm^-3]
phys_a = 1.0

alpha = 4 * phys_a / epsilon 

#phys_invh = 1/phys_h
D = zeros(Ncells)

# Material values
T0 = 1.0 # 
bee = zeros(Ncells) # Heat Capacity


# ========================
# Random Number Initialization
# ========================
Random.seed!(12345)

# ========================
# Main Timeloop
# ========================

# Start Loop Until Time Ends:

# Update Temperature Dependent Quantities:
# Opacity, Beta, Fleck Factor
function update()
    # Calculate Planck mean opacity for the mesh

    bee[:] = alpha * mesh_temp[:].^3

    mesh_sigma[:] = sigma_a * ones(Ncells)

    # Calculate beta (non-linearity) factor
    #mesh_beta[:] = 4.0 * phys_a .* (mesh_temp[:].^3) / bee

    mesh_beta[:] = ones(Ncells) # Result simplifies to 1 so no need to calculate it

    #mesh_beta[:] = 4.0 * phys_a .* (mesh_temp[:].^3)

    # Calculate Fleck factor
    mesh_fleck[:] = 1.0 ./ (1.0 .+ alpha * mesh_beta[:] * phys_c * dt .* mesh_sigma[:])

end

# Source New IMC Particles:
# Energy terms, Emission probabilities,
# Number of particles, Particle creation
function sourcing(t)    
    """Energy source terms"""
    
    sources = zeros(Ncells)
    sources[1:sourceCells+1] .= 1.0

    # Check if source is still active
    if t <= tau_0
        e_body = (sigma_a * mesh_fleck[:] .* phys_c * dx * dt * phys_a .* (mesh_temp[:].^4)) .+ (sources[:] * dt * dx)
    else
        e_body = (sigma_a * mesh_fleck[:] .* phys_c * dx * dt * phys_a .* (mesh_temp[:].^4))
    end

    e_total = sum(e_body[:]) # Total Energy

    print("The energy in sourcing is: ", e_total, "\n")

    push!(source_energies, e_total)
    
    """Calculate number of particles to source on the surface and throughout the mesh"""
    n_sourcing = n_input
    n_census = length(particles)
    if n_input + n_census > n_max
        n_sourcing = n_max - n_census - Ncells - 1
    end

    n_source = zeros(Int64, Ncells)
    n_body = zeros(Int64, Ncells)

    # Fleck & Cummings Deterministic Sourcing
    # n_surf = round(e_surf*n_source/e_total)+1

    for cellindex in 1:Ncells
        if mesh_temp[cellindex] > 0 || cellindex <= sourceCells+1
            n_body[cellindex] = round(e_body[cellindex]*n_sourcing/e_total)+1
        end
    end

    """Create source particles"""

    # Body-source particles
    for cellindex in 1:Ncells
        if n_body[cellindex] <= 0
            continue
        end
        nrg = e_body[cellindex] / Float64(n_body[cellindex])
        startnrg = nrg

        for j in 1:n_body[cellindex]
            origin = cellindex
            xpos = dx * rand()
            mu = 1.0 - 2.0 * rand()
            spawntime = t + dt * rand()
            frq = 1.0 # Frequency not used so set to arbitrary value
            push!(particles, [origin, spawntime, cellindex, xpos, mu, frq, nrg, startnrg])
        end
    end

end


# Monte Carlo Simulation:
function MC(t)
    endsteptime = t + dt 
    PLE = zeros(Ncells)
    mesh_energydep[:] = zeros(Ncells)
    absorbed_particles = 0

    for particle in 1:length(particles)
    
        currentparticle = particles[particle, :] 

        origin = currentparticle[1][1]
        currenttime = currentparticle[1][2]
        cellindex = Int64(currentparticle[1][3])
        position = currentparticle[1][4] 
        mu = currentparticle[1][5]
        freq = currentparticle[1][6]
        energy = currentparticle[1][7]
        startenergy = currentparticle[1][8]

        minenergy = 0.01 * startenergy

        iterations = 0
        while true
            iterations += 1

            # Calculate distances to boundary, collision, and census -> take minimum

            # Boundary distance
            if mu > 0.0
                dist_b = (dx - position)/mu
            else
                dist_b = abs(position/mu)
            end
        
            # Collision distance (Scattering)
            dist_col = abs(randexp()) / (sigma_s + (1-mesh_fleck[cellindex])*sigma_a) # Effective scattering includes isotropic term
            # Census distance
            dist_cen = phys_c * (endsteptime - currenttime)

            # Actual distance - closest distnace
            dist = min(dist_b, dist_col, dist_cen)

            # Calculate new energy and deposit lost energy
            newenergy = energy * exp(-mesh_fleck[cellindex] * sigma_a * dist)
            #newenergy = energy * exp(-sigma_a * dist)
            if newenergy <= minenergy
                newenergy = 0.0
            end

            # Deposit the particle's energy
            mesh_energydep[cellindex] += energy - newenergy

            # Advance position, time, and energy
            # If energy is zero or domain boundary crossed -> kill particle
            if newenergy == 0.0
                # Flag particle for later destruction
                absorbed_particles += 1
                particles[particle][8] = -1.0
                break
            end

            # Otherwise, advance the position, time and energy
            position += mu * dist
            currenttime += dist / phys_c 
            energy = newenergy

            # If the event was a boundary-crossing, and the boundary is the
            # domain boundary, then kill the particle
            if dist == dist_b
                if mu > 0
                    if cellindex == Ncells # Right-boundary
                        # Flag particle for later destruction
                        particles[particle][8] = -1.0
                        break
                    end
                    cellindex += 1
                    position = 0 
                end
                if mu < 0
                    if cellindex == 1 # Left-boundary
                        # Reflecting boundary
                        mu = -mu
                    else
                        cellindex -= 1
                        position = dx
                    end
                end
            end
         

            # If collision occured update frequency and direction -> return to cross-section calculation
            if dist == dist_col
            # Collision (i.e. absorption, but treated as pseudo-scattering)
                #freq = mesh_temp[cellindex] * randexp()
                mu = 1.0 - 2.0 * rand()
                end
            # If census event occured, finish history and update particle properties in list
            if dist == dist_cen
                # Finished with this particle
                # Update the particle's properties in the list
                # Starting energy doesn't change
                particles[particle][:] = [origin, currenttime, cellindex, position, mu, freq, energy, startenergy]
                break
            end
        end 
    # New particle history
    #print("Particle state at end of time-step ", particles[particle][:], "\n")
    end
    #print("The number of particles that were absorbed is ", absorbed_particles, "\n")
end

function clean()
    # Removes particles that died during the simulation
    energy_loss = 0
    for prtclindex in length(particles):-1:1 # Loop backwards
        if particles[prtclindex][8] == -1.0
            energy_loss += particles[prtclindex][7]
            deleteat!(particles, prtclindex)
            #global n_census -= 1
        end
    end
    #print("The energy lost during clean-up is ", energy_loss, "\n")
    push!(lost_energy, energy_loss)
    #print("The number of particles after clean-up is ", length(particles), "\n")
end


# Tally End of Timestep Quantities:
# Radiation energy density, increase Temperature
function tally(t)
    """Tally end of timestep quantities """
    # Radiation energy density
    radenergydens = zeros(Ncells)
    radenergydens[:] = phys_a * mesh_temp[:].^4 # keV/cm^3

    sourceenergy = zeros(Ncells)
    sourceenergy[1:sourceCells+1] .= 1.0
    
    # Temperature increase
    nrg_inc = zeros(Ncells)
    nrg_inc[:] = (mesh_energydep[:]/dx) .- (sigma_a * mesh_fleck[:] .* radenergydens[:] * phys_c * dt) 
  
    material_energy[:] = material_energy[:] .+ nrg_inc[:]

    for i in 1:Ncells
        if material_energy[i] < 0
            material_energy[i] = 0
        end
    end

    print("Material energy ", maximum(material_energy), "\n")

    mesh_temp[:] = material_energy[:].^(1/4) # Can invert material energy to recover temperature due to linearization in Su Olson problem

    print("The mesh temperature is: \n")
    print(sum(mesh_temp[:]), "\n")

    push!(mesh_tempsaved, copy(mesh_temp)) # Add the current mesh temperature to a saved list

    push!(materialenergy_saved, material_energy[:]) # Add the current material energy to a saved list

    #print("Mesh temperature: \n")
    #print(mesh_temp, "\n")

    total_energyincrease = sum(nrg_inc[:])/ Ncells
    #print("The total energy increase was ", total_energyincrease, "\n")
    push!(energy_increases, total_energyincrease)

    # Save radiation energy
    radnrg = zeros(Ncells)

    for ii in 1:length(particles)
        radindex = Int(particles[ii][3])
        radnrg[radindex] += particles[ii][7]/dx
    end

    print("Radiation energy ", maximum(radnrg), "\n")

    push!(radenergy_saved, copy(radnrg))
    #print("Radiation energy \n", radnrg, "\n")
end

# ========================
# Energy Check
# ========================

function energycheck()
    """ Check for energy conservation """
    # Input energy
    # radenergydens[:] = phys_a * mesh_temp[:].^4
    # input_energy = sum(radenergydens)
    # print("The input energy after is: ", input_energy, "\n")

    # # # Internal energy (Radiation and Material)
    total_matenergy = sum(mesh_energydep)
    #print("The material energy is: ", total_matenergy, "\n")
   
    total_partenergy = 0
    for ii in 1:length(particles)
        total_partenergy += particles[ii][7]
    end
    #print("The particle energy is: ", total_partenergy, "\n")

    #print("Energy check ", (stored_energy[end] - lost_energy[end] - total_partenergy) / stored_energy[end], "\n")

    if(length(source_energies)) > 1
        print("Energy check: ", 100*(source_energies[end] - source_energies[end-1] - energy_increases[end-1])/source_energies[end], "% error\n")
    end

    #print("Energy check: ", old_sourceenergy - e_total, "\n")
    #internal_energy = total_matenergy + total_partenergy
    #print("The internal energy is: ", internal_energy, "\n")

end



# ========================
# Driver script
# ========================


function main()
    t = 0.0 # Start time
    print("Starting IMC Calculations \n")
    while t <= tmax
        print("Simulating time t = ", t, "\n") 
        update()  # Update values at start of time-step
        
        sourcing(t) # Source particles
        
        MC(t) # Monte Carlo simulation of time-step

        clean() # Remove dead particles from simulation

        tally(t) # Tally end of time-step quantities
        
        energycheck() # Check for energy conservation

        # Increment time-step
        t += dt
    end
end

# Start simulation and time it
runtime = @elapsed begin
    main()
end

print("The simulation took ", runtime, " seconds to run \n")

# ========================
# Results
# ========================

"""Benchmark Solutions"""
x_bench = [0.01000, 0.10000, 0.17783, 0.31623, 0.45000, 0.50000, 0.56234, 0.75000, 1.00000, 1.33352, 1.77828, 3.16228, 5.62341, 10.00000, 17.78279]

t_bench = [0.10000, 0.31623, 1.00000, 3.16228, 10.00000, 31.6228, 100.000]

rad_benchone = [[0.09531, 0.09531, 0.09532, 0.09529, 0.08823, 0.04765, 0.00375, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000], 
                [0.27526, 0.27526, 0.27527, 0.26262, 0.20312, 0.13762, 0.06277, 0.00280, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000],
                [0.64308, 0.63585, 0.61958, 0.56187, 0.44711, 0.35801, 0.25374, 0.11430, 0.03648, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000],
                [1.20052, 1.18869, 1.16190, 1.07175, 0.90951, 0.79902, 0.66678, 0.44675, .027540, 0.14531, 0.05968, 0.00123, 0.00000, 0.00000, 0.00000],
                [2.23575, 2.21944, 2.18344, 2.06448, 1.86072, 1.73178, 1.57496, 1.27398, 0.98782, 0.70822, 0.45016, 0.09673, 0.00375, 0.00000, 0.00000],
                [0.69020, 0.68974, 0.68878, 0.68569, 0.68111, 0.67908, 0.67619, 0.66548, 0.64691, 0.61538, 0.56353, 0.36965, 0.10830, 0.00390, 0.00000],
                [0.35720, 0.35714, 0.35702, 0.35664, 0.35599, 0.35574, 0.35538, 0.35393, 0.35141, 0.34697, 0.33924, 0.30346, 0.21382, 0.07200, 0.00272]
               ] # ca = 1, cs = 0 

rad_benchtwo = [[0.09757, 0.09757, 0.09758, 0.09756, 0.09033, 0.04878, 0.00383, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000], 
                [0.29363, 0.29365, 0.29364, 0.28024, 0.21573, 0.14681, 0.06783, 0.00292, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000],
                [0.72799, 0.71888, 0.69974, 0.63203, 0.50315, 0.40769, 0.29612, 0.13756, 0.04396, 0.00324, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000],
                [1.28138, 1.26929, 1.24193, 1.15018, 0.98599, 0.87477, 0.74142, 0.51563, 0.33319, 0.18673, 0.08229, 0.00160, 0.00000, 0.00000, 0.00000],
                [2.26474, 2.24858, 2.21291, 2.09496, 1.89259, 1.76429, 1.60822, 1.30947, 1.02559, 0.74721, 0.48739, 0.11641, 0.00554, 0.00000, 0.00000],
                [0.68703, 0.68656, 0.68556, 0.68235, 0.67761, 0.67550, 0.67252, 0.66146, 0.64239, 0.61024, 0.55789, 0.36631, 0.11177, 0.00491, 0.00000],
                [0.35675, 0.35668, 0.35654, 0.35618, 0.35552, 0.35527, 0.35491, 0.35346, 0.35092, 0.34646, 0.33868, 0.30281, 0.21323, 0.07236, 0.00296]
               ] # ca = 0.5, cs = 0.5

mat_benchone = [[0.00468, 0.00468, 0.00468, 0.00468, 0.00455, 0.00234, 0.00005, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000], 
                [0.04093, 0.04093, 0.04093, 0.04032, 0.03314, 0.02046, 0.00635, 0.00005, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000],
                [0.27126, 0.26839, 0.26261, 0.23978, 0.18826, 0.14187, 0.08838, 0.03014, 0.00625, 0.00017, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000],
                [0.94670, 0.93712, 0.91525, 0.84082, 0.70286, 0.60492, 0.48843, 0.30656, 0.17519, 0.08352, 0.02935, 0.00025, 0.00000, 0.00000, 0.00000],
                [2.11186, 2.09585, 2.06052, 1.94365, 1.74291, 1.61536, 1.46027, 1.16591, 0.88992, 0.62521, 0.38688, 0.07642, 0.00253, 0.00000, 0.00000],
                [0.70499, 0.70452, 0.70348, 0.70020, 0.69532, 0.69308, 0.68994, 0.67850, 0.65868, 0.62507, 0.57003, 0.36727, 0.10312, 0.00342, 0.00000],
                [0.35914, 0.35908, 0.35895, 0.35854, 0.35793, 0.35766, 0.35728, 0.35581, 0.35326, 0.34875, 0.34086, 0.30517, 0.21377, 0.07122, 0.00261]
               ] # ca = 1, cs = 0
       
mat_benchtwo = [[0.00242, 0.00242, 0.00242, 0.00242, 0.00235, 0.00121, 0.00003, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000],
                [0.02255, 0.02253, 0.02256, 0.02223, 0.01826, 0.01128, 0.00350, 0.00003, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000],
                [0.17609, 0.17420, 0.17035, 0.15520, 0.12164, 0.09194, 0.05765, 0.01954, 0.00390, 0.00009, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000],
                [0.77654, 0.76878, 0.75108, 0.69082, 0.57895, 0.49902, 0.40399, 0.25610, 0.14829, 0.07161, 0.02519, 0.00018, 0.00000, 0.00000, 0.00000],
                [2.00183, 1.98657, 1.95286, 1.84104, 1.64778, 1.52383, 1.37351, 1.09216, 0.83248, 0.58640, 0.36629, 0.07658, 0.00290, 0.00000, 0.00000],
                [0.71860, 0.71805, 0.71687, 0.71312, 0.70755, 0.70499, 0.70144, 0.68851, 0.66637, 0.62937, 0.57001, 0.36066, 0.10181, 0.00385, 0.00000],
                [0.36067, 0.36065, 0.36047, 0.36005, 0.35945, 0.35917, 0.35876, 0.35727, 0.35465, 0.35004, 0.34200, 0.30553, 0.21308, 0.07077, 0.00273]
               ] # ca = 0.5, cs = 0.5

print(typeof(mesh_centers))




"""Create plots at different times matching the benchmark results"""

# Radiation energy density

# Benchmark values currently set to plot at times 10, 1, 0.1

plot(mesh_centers, radenergy_saved[Int(timesteps),:], marker="o", title="Double Precision IMC Results for dt = "*string(dt)*"", label="tau = "*string(round(phys_c * dt * timesteps, sigdigits=2))*" IMC", size=(600,600), margin=5Plots.mm)
plot!(x_bench, rad_benchtwo[5,:], markershape=:xcross, label="tau = "*string(phys_c * dt * timesteps)*" Bench")

xlabel!("x")
ylabel!("Radiation Energy Density")
xlims!(0.0, 4.0)
#ylims!(0.0, 1.2)
xticks!(0:0.5:4)
#yticks!(0:0.2:1.2)

plot!(mesh_centers, radenergy_saved[Int(timesteps/10),:], marker="o", label="tau = "*string(phys_c * dt * timesteps/10)*" IMC")
plot!(x_bench, rad_benchtwo[3,:], markershape=:xcross, label="tau = "*string(phys_c * dt * timesteps/10)*" Bench")

plot!(mesh_centers, radenergy_saved[Int(timesteps/100),:], marker="o", label="tau = "*string(phys_c * dt * timesteps/100)*" IMC")
display(plot!(x_bench, rad_benchtwo[1,:], markershape=:xcross, label="tau = "*string(phys_c * dt * timesteps/100)*" Bench"))


#savefig("Radiation Energy Density Double.png")

#writedlm("RadEnergyResultsDouble.csv", radenergy_saved, ',')

# Material energy density

plot(mesh_centers, materialenergy_saved[Int(timesteps),:], marker="o", xscale=:log10, yscale=:log10, minorgrid=true, title="Double Precision IMC Results for dt = "*string(dt)*"", label="tau = "*string(round(phys_c * dt * timesteps, sigdigits=2))*" IMC", size=(600,600), margin=5Plots.mm)
plot!(x_bench, mat_benchtwo[5,:], markershape=:xcross, label="tau = "*string(phys_c * dt * timesteps)*" Bench")

xlabel!("x")
ylabel!("Material Energy Density")
#xlims!(0.0, 30.0) 
ylims!(1e-3, 1e2)
#xticks!(0:0.5:4)
#yticks!(0:0.2:1.2)

plot!(mesh_centers, materialenergy_saved[Int(timesteps/10),:], marker="o", label="tau = "*string(phys_c * dt * timesteps/10)*" IMC")
plot!(x_bench, mat_benchtwo[3,:], markershape=:xcross, label="tau = "*string(phys_c * dt * timesteps/10)*" Bench")

plot!(mesh_centers, materialenergy_saved[Int(timesteps/100),:], marker="o", label="tau = "*string(phys_c * dt * timesteps/100)*" IMC")
display(plot!(x_bench, mat_benchtwo[1,:], markershape=:xcross, label="tau = "*string(phys_c * dt * timesteps/100)*" Bench"))


#savefig("Material Energy Density Double.png")

#writedlm("MatEnergyResultsDouble.csv", materialenergy_saved, ',')



