# FUNCTIONS FOR THE MOLECULAR DYNAMICS PROGRAM
import numpy as _np

class Simulator:
    def __init__(self, Nstep, Nparticles, Tstep, Temperature, Mu=None, Coupling=None, dp=None, mass=1.0):
        # Simulation data
        self.Data_pos  = _np.empty(shape=(Nstep,Nparticles))
        self.Data_vel  = _np.empty(shape=(Nstep,Nparticles))
        self.Data_time = _np.empty(shape=Nstep)
        self.Data_step = _np.empty(shape=Nstep)
        self.Data_Epot = _np.empty(shape=(Nstep,Nparticles))
        self.Data_Ekin = _np.empty(shape=(Nstep,Nparticles))
        self.Data_Etot = _np.empty(shape=(Nstep,Nparticles))
        self.Data_Temp = _np.empty(shape=Nstep)

        # Simulation variables
        self.Position = _np.zeros(Nparticles)
        self.Velocity = _np.zeros(Nparticles)
        self.Force    = _np.zeros(Nparticles)
        self.U        = _np.zeros(Nparticles)
        
        # Simulation constants
        self.Nparticles = Nparticles
        self.Tstep = Tstep
        self.m = mass
        self.Temperature = Temperature
        self.kB = 0.00831446                                    # [Dalton*nm**2 /(K*ps**2) = kJ/(mol*K)]
        self.Mu = Mu
        self.Coupling = Coupling
        self.dp = dp
        self.it = 0
        self.Accept = 0
    
    def VelocityVerlet_NVE(self, Potential):
        '''Velocity Verlet algorithm using half step with no thermostat (NVE).'''
        self.Velocity += 0.5*self.Tstep*self.Force/self.m       # [nm/ps]
        self.Position += self.Tstep*self.Velocity               # [nm]
        
        vecPotential = _np.vectorize(Potential, cache=False)
        self.U, self.Force = vecPotential(self.Position)        # [kJ/mol], [kJ/(nm*mol)]=[Dalton*nm/ps**2]
        
        self.Velocity += 0.5*self.Tstep*self.Force/self.m       # [nm/ps]

    def Langevin_BAOAB(self, Potential):
        '''Symmetric Langevin Velocity-Verlet method (NVT).'''
        self.Velocity += 0.5*self.Tstep*self.Force/self.m       # Half step (B)
        self.Position += 0.5*self.Tstep*self.Velocity           # Half step (A)                                         
        c = _np.exp(-self.Mu*self.Tstep)                        # Weak solve of Ornstein-Uhlenbeck process (O)
        self.Velocity = (c * self.Velocity) + (_np.sqrt((1-c*c)*self.kB*self.Temperature)*_np.sqrt(self.m)*_np.random.normal(0,1,1))/self.m
        self.Position += 0.5*self.Tstep*self.Velocity           # Half step (A)
        vecPotential = _np.vectorize(Potential, cache=False)
        self.U, self.Force = vecPotential(self.Position)        # [kJ/mol], [kJ/(nm*mol)]=[Dalton*nm/ps**2]
        self.Velocity += 0.5*self.Tstep*self.Force/self.m       # Half step (B)

    def Langevin_ABOBA(self, Potential):
        '''Symmetric Langevin Position-Verlet method (NVT).'''
        self.Position += 0.5*self.Tstep*self.Velocity           # Half step (A)
        vecPotential = _np.vectorize(Potential, cache=False)
        self.U, self.Force = vecPotential(self.Position)        # [kJ/mol], [kJ/(nm*mol)]=[Dalton*nm/ps**2]
        self.Velocity += 0.5*self.Tstep*self.Force/self.m       # Half step (B)                                        
        c = _np.exp(-self.Mu*self.Tstep)                        # Weak solve of Ornstein-Uhlenbeck process (O)
        self.Velocity = (c * self.Velocity) + (_np.sqrt((1-c*c)*self.kB*self.Temperature)*_np.sqrt(self.m)*_np.random.normal(0,1,1))/self.m
        self.Velocity += 0.5*self.Tstep*self.Force/self.m       # Half step (B)
        self.Position += 0.5*self.Tstep*self.Velocity           # Half step (A)

    def Langevin_OBABO(self, Potential):
        '''Bussi-Parrinello Langevin method (NVT).'''
        c = _np.exp(-self.Mu*self.Tstep)                        # Weak solve of Ornstein-Uhlenbeck process, half step (O)
        self.Velocity = (c * self.Velocity) + (_np.sqrt((1-c*c)*self.kB*self.Temperature)*_np.sqrt(self.m)*_np.random.normal(0,1,1))/self.m
        self.Velocity += 0.5*self.Tstep*self.Force/self.m       # Half step (B)
        self.Position += self.Tstep*self.Velocity               # Full step (A)
        vecPotential = _np.vectorize(Potential, cache=False)
        self.U, self.Force = vecPotential(self.Position)        # [kJ/mol], [kJ/(nm*mol)]=[Dalton*nm/ps**2]
        self.Velocity += 0.5*self.Tstep*self.Force/self.m       # Half step (B)
        self.Velocity = (c * self.Velocity) + (_np.sqrt((1-c*c)*self.kB*self.Temperature)*_np.sqrt(self.m)*_np.random.normal(0,1,1))/self.m # Half step (O)
        
    def Thermostat_And(self):
        '''Andersen thermostat (NVT).'''
        for particle in range(self.Nparticles):
            # Check for stochastic collision
            if self.Mu * self.Tstep > _np.random.uniform(0,1):
                self.Velocity[particle] = self.RandomVelocity(self.Temperature,self.m)  # [nm/ps]

    def Thermostat_Ber(self):
        '''Berendsen thermostat (NVT).'''
        TInstant = self.m/self.kB*_np.mean(self.Velocity**2)   # [K]
        Lambda = _np.sqrt(1+self.Tstep/self.Coupling*((self.Temperature/TInstant)-1))
        self.Velocity *= Lambda    # [nm/ps]

    def MonteCarlo(self, i, Potential):
        '''Monte Carlo simulation (NVT).'''
        Uold, Fold = Potential(self.Position[i])
        Position_new = self.Position[i] + (_np.random.uniform(0.0, 1.0) - 0.5) * self.dp 
        Unew, Fnew = Potential(Position_new)
        self.U[i] = Uold
        
        # Metropolis Criteria
        if _np.exp(-(Unew-Uold)/(self.kB*self.Temperature)) >= _np.random.uniform(0.0, 1.0):
            self.Position[i] = Position_new
            self.U[i] = Unew
            self.Accept+=1.0

    def RandomVelocity(self, Temperature, m, size=1):
        '''Random velocity asignment from the Maxwell-Boltzmann distribution.'''
        return _np.sqrt(self.kB*self.Temperature/m) * _np.random.normal(0, 1, size)

    def SampleData(self):
        '''Collection of data after each MD/MC itteration.'''
        self.Data_pos[self.it] = self.Position                          # [nm]
        self.Data_vel[self.it] = self.Velocity                          # [nm/ps]
        self.Data_time[self.it] = self.it*self.Tstep                          # [ps]
        self.Data_step[self.it] = self.it                                # Unitless
        self.Data_Epot[self.it] = self.U                                # [kJ/mol]
        self.Data_Ekin[self.it] = 0.5*self.m*self.Velocity**2                # [kJ/mol]
        self.Data_Etot[self.it] = self.U+0.5*self.m*self.Velocity**2              # [kJ/mol]
        self.Data_Temp[self.it] = self.m/self.kB*_np.mean(self.Velocity**2)        # [K]
        self.it += 1
