# FUNCTIONS FOR THE MOLECULAR DYNAMICS PROGRAM
import numpy as _np

class Simulator:
    def __init__(self, Nstep, Nparticles, Tstep, Temperature, MU=None, Coupling=None, dp=None):
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
        self.m = 1
        self.Temperature = Temperature
        self.kB = 0.00831446                             # [Dalton*nm**2 /(K*ps**2) = kJ/(mol*K)]
        self.MU = MU
        self.Coupling = Coupling
        self.dp = dp
        self.it = 0
        self.Accept = 0
    
    def VelocityVerlet_NVE(self, i, Potential):
        '''Velocity Verlet algorithm using half step with no thermostat (NVE).'''
        self.Velocity[i] += 0.5*self.Tstep*self.Force[i]/self.m                 # [nm/ps]
        self.Position[i] += self.Tstep*self.Velocity[i]                    # [nm]
        self.U[i], self.Force[i] = Potential(self.Position[i])           # [kJ/mol], [kJ/(nm*mol)]=[Dalton*nm/ps**2]
        self.Velocity[i] += 0.5*self.Tstep*self.Force[i]/self.m                 # [nm/ps]

    def Langevin_BAOAB(self, Position, Velocity, Force, Potential):
        '''Symmetric Langevin Velocity-Verlet method (NVT).'''
        Velocity += 0.5*self.Tstep*Force/self.m                # Half step (B)
        Position += 0.5*self.Tstep*Velocity               # Half step (A)
                                                        
        c = _np.exp(-self.MU*self.Tstep)                        # Weak solve of Ornstein-Uhlenbeck process (O)
        Velocity = (c * Velocity) + (_np.sqrt((1-c*c)*self.kB*self.Temperature)*_np.sqrt(self.m)*_np.random.normal(0,1,1))/self.m
        
        Position += 0.5*self.Tstep*Velocity               # Half step (A)
        Energy, Force = Potential(Position)          # [kJ/mol], [kJ/(nm*mol)]=[Dalton*nm/ps**2]
        Velocity += 0.5*self.Tstep*Force/self.m                # Half step (B)
        
        return Position, Velocity, Energy, Force

    def Langevin_ABOBA(self, Position, Velocity, Force, Potential):
        '''Symmetric Langevin Position-Verlet method (NVT).'''
        Position += 0.5*self.Tstep*Velocity               # Half step (A)
        Energy, Force = Potential(Position)          # [kJ/mol], [kJ/(nm*mol)]=[Dalton*nm/ps**2]
        Velocity += 0.5*self.Tstep*Force/self.m                # Half step (B)

                                                    
        c = _np.exp(-self.MU*self.Tstep)                        # Weak solve of Ornstein-Uhlenbeck process (O)
        Velocity = (c * Velocity) + (_np.sqrt((1-c*c)*self.kB*self.Temperature)*_np.sqrt(m)*_np.random.normal(0,1,1))/m
        
        Velocity += 0.5*self.Tstep*Force/self.m                # Half step (B)
        Position += 0.5*self.Tstep*Velocity               # Half step (A)
        
        return Position, Velocity, Energy, Force

    def Langevin_OBABO(self, Position, Velocity, Force, Potential):
        '''Bussi-Parrinello Langevin method (NVT).'''
        c = _np.exp(-self.MU*self.Tstep)                        # Weak solve of Ornstein-Uhlenbeck process, half step (O)
        Velocity = (c * Velocity) + (_np.sqrt((1-c*c)*self.kB*self.Temperature)*_np.sqrt(self.m)*_np.random.normal(0,1,1))/self.m
        Velocity += 0.5*self.Tstep*Force/self.m                # Half step (B)
        
        Position += self.Tstep*Velocity                   # Full step (A)
        Energy, Force = Potential(Position)          # [kJ/mol], [kJ/(nm*mol)]=[Dalton*nm/ps**2]
        
        Velocity += 0.5*self.Tstep*Force/self.m                # Half step (B)
        Velocity = (c * Velocity) + (_np.sqrt((1-c*c)*self.kB*self.Temperature)*_np.sqrt(self.m)*_np.random.normal(0,1,1))/self.m # Half step (O)
        
        return Position, Velocity, Energy, Force
        
    def Thermostat_And(self):
        '''Andersen thermostat (NVT).'''
        for particle in range(self.Nparticles):
            # Check for stochastic collision
            if self.MU * self.Tstep > _np.random.uniform(0,1):
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

    def RandomVelocity(self, Temperature, m):
        '''Random velocity asignment from the Maxwell-Boltzmann distribution.'''
        return _np.sqrt(self.kB*self.Temperature/m) * _np.random.normal(0,1,1)

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
