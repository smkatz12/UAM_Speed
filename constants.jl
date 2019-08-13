export COC,WD,WA,SD,SA, stateType, actType, ACTIONS, discount_f, RANGES,THETAS,PSIS,OWNSPEEDS,INTRPSEEDS, interp, turns

g = 9.81

# ADVISORY INDICES
COC=0
WD=1
WA=2
SD=3
SA=4

# State Type:
stateType = Tuple{Float64,Float64,Float64,Float64,Float64,Int}
actType = Int
ACTIONS = [COC,WD,WA,SD,SA]

# Default parameters
discount_f = 1.0

### STATE CUTPOINTS ###
RANGES = [0.0,25.0,50.0,75.0,100.0,150.0,200.0,300.0,400.0,500.0,510.0,750.0,1000.0,1500.0,2000.0,3000.0,4000.0,5000.0,7000.0,9000.0,11000.0,13000.0,15000.0]
THETAS = Array(LinRange(-pi,pi,41))
PSIS   = Array(LinRange(-pi,pi,41))
OWNSPEEDS = [0.0, 25.0, 50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0, 300.0]   #ft/s
INTRSPEEDS = [0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0] #ft/s

interp = LocalGIFunctionApproximator(RectangleGrid(RANGES,THETAS,PSIS,OWNSPEEDS,INTRSPEEDS,ACTIONS)) # Create the local function approximator using the grid


### Dictionaries to define transitions ###
probs = [0.5,0.25,0.25]
accels = Dict(COC=>([0.34,0.33,0.33],[0.0,-0.01g,0.01g]),
              WD=>(probs,[-0.04g,-0.035g,-0.045g]),
              WA=>(probs,[0.04g,0.035g,0.045g]),
              SD=>(probs,[-0.08g,-0.075g,-0.085g]),
              SA=>(probs,[0.08g,0.075g,0.085g]),
              -1=>([0.34,0.33,0.33],[0.0,-0.01g,0.01g])) # FOR v5, 0, 1, -1