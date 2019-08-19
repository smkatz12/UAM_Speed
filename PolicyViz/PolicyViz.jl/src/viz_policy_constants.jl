#Regression Constants
export PRAS,INTRSPEEDS,THETAS,PSIS,OWNSPEEDS,ACTIONS,VERTICAL_TAUS,RANGES, nstates

mps2fps = 3.28084

### STATE CUTPOINTS ###
RANGES = [0.0,25.0,50.0,75.0,100.0,150.0,200.0,300.0,400.0,500.0,510.0,750.0,1000.0,1500.0,2000.0,3000.0,4000.0,5000.0,7000.0,9000.0,11000.0,13000.0,15000.0] #ft
THETAS = Array(LinRange(-π,π,21))
PSIS   = Array(LinRange(-π,π,21))
OWNSPEEDS = [0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 28.0, 32.0, 36.0].*mps2fps   #ft/s
INTRSPEEDS = [0.0, 8.0, 16.0, 24.0, 32.0, 40.0].*mps2fps #ft/s
const VERTICAL_TAUS   = collect(0:dt:nTau_max)
const PRAS    = [0, 1, 2, 3, 4]

const RANGEMAX = maximum(RANGES)/1000.0 #kft
const NSTATES = length(RANGES)*length(THETAS)*length(PSIS)*length(OWNSPEEDS)*length(INTRSPEEDS)
const ACTIONS = deg2rad.([0.0, -1.5, 1.5, -3.0, 3.0])