# 2460066.500000000 = A.D. 2023-May-02 00:00:00.0000 TDB
{sun23 = "X =-1.324316846741195E+06 Y =-1.446765512621427E+05 Z = 3.202738040537933E+04
 VX= 4.021534124461124E-03 VY=-1.472929900704240E-02 VZ= 2.948517133131729E-05"
mer23 = "X =-5.140210130150759E+07 Y =-4.381189053847280E+07 Z = 1.056887671139676E+06
 VX= 2.207720296247787E+01 VY=-3.453839718779295E+01 VZ=-4.845914556700924E+00"
ven23 = "X =-9.848082136897813E+07 Y = 4.552713536799312E+07 Z = 6.265256376415793E+06
 VX=-1.504882769689358E+01 VY=-3.187376111911025E+01 VZ= 4.312134646762669E-01"
ear23 = "X =-1.150223276469004E+08 Y =-9.909841274895248E+07 Z = 3.755860068777204E+04
 VX= 1.907614272241964E+01 VY=-2.258721645907844E+01 VZ= 2.018298546438757E-03"
mar23 = "X =-2.010819716042931E+08 Y = 1.476265433762218E+08 Z = 8.029003272268474E+06
 VX=-1.349413909238045E+01 VY=-1.742820367074673E+01 VZ=-3.381872148509668E-02"
jup23 = "X = 6.792848124845287E+08 Y = 2.926757289590889E+08 Z =-1.641170746619791E+07
 VX=-5.318068089405075E+00 VY= 1.261385267190711E+01 VZ= 6.658478320447347E-02"
sat23 = "X = 1.265581460039990E+09 Y =-7.382676489447224E+08 Z =-3.755222856885275E+07
 VX= 4.325577300311080E+00 VY= 8.324810463771845E+00 VZ=-3.167483277652900E-01"
ura23 = "X = 1.945406877724633E+09 Y = 2.202709703800068E+09 Z =-1.702215196928215E+07
 VX=-5.154055480935001E+00 VY= 4.190735053981157E+00 VZ= 8.237035843708407E-02"
nep23 = "X = 4.455662623587698E+09 Y =-3.833067895772491E+08 Z =-9.479194396546717E+07
 VX= 4.297305233379821E-01 VY= 5.447227229797033E+00 VZ=-1.222099521045787E-01"}

# 2459701.500000000 = A.D. 2022-May-02 00:00:00.0000 TDB
{sun22 = "X =-1.333512273509328E+06 Y = 3.445488642388455E+05 Z = 2.832491561644843E+04
 VX=-3.688954247616953E-03 VY=-1.535801956098796E-02 VZ= 2.095356489194861E-04"
  mer22 = "X =-5.787964463498133E+07 Y = 5.001123068678558E+06 Z = 5.595649085666763E+06
 VX=-1.412786024627931E+01 VY=-4.647259111665712E+01 VZ=-2.500663888122642E+00"
  ven22 = "X = 3.587612704045174E+07 Y =-1.019279819610429E+08 Z =-3.522707892272152E+06
 VX= 3.267015721969536E+01 VY= 1.184076411660992E+01 VZ=-1.722399554824109E+00"
  ear22 = "X =-1.146121503189677E+08 Y =-9.911837382505153E+07 Z = 3.390113022846729E+04
 VX= 1.917565337627882E+01 VY=-2.252086620472453E+01 VZ= 1.293657910466095E-04"
  mar22 = "X = 1.168429386885249E+08 Y =-1.723122647793255E+08 Z =-6.489026128990121E+06
 VX= 2.090480902530879E+01 VY= 1.575033363242685E+01 VZ=-1.822586427078310E-01"
  jup22 = "X = 7.305933493985869E+08 Y =-1.301920721371184E+08 Z =-1.580507858542991E+07
 VX= 2.137303769034409E+00 VY= 1.347648922277490E+01 VZ=-1.037029936830756E-01"
  sat22 = "X = 1.103882677050047E+09 Y =-9.837166526117289E+08 Z =-2.684592972538310E+07
 VX= 5.884999887995806E+00 VY= 7.191307123896053E+00 VZ=-3.588311673691735E-01"
  ura22 = "X = 2.102764999640690E+09 Y = 2.064949098017363E+09 Z =-1.957237631821322E+07
 VX=-4.821150907798146E+00 VY= 4.541582713253250E+00 VZ= 7.931778087597818E-02"
  nep22 = "X = 4.438820542662540E+09 Y =-5.547741023157761E+08 Z =-9.087265524840316E+07
 VX= 6.385999387923026E-01 VY= 5.425372291426413E+00 VZ=-1.270328159123488E-01"}

# This code turns the strings of data into usable vectors
stringtovector = function(string) {
  words = paste(c("Y =", "Z =", "\n VX=", "VY=", "VZ="), collapse = "|")
  string = trimws(gsub(words, "split", string))
  string = gsub("X =", "", string)
  string = as.numeric(strsplit(string, "split")[[1]])
  return(string)
}

{sun23 = stringtovector(sun23)
  mer23 = stringtovector(mer23)
  ven23 = stringtovector(ven23)
  ear23 = stringtovector(ear23)
  mar23 = stringtovector(mar23)
  jup23 = stringtovector(jup23)
  sat23 = stringtovector(sat23)
  ura23 = stringtovector(ura23)
  nep23 = stringtovector(nep23)
  sun22 = stringtovector(sun22)
  mer22 = stringtovector(mer22)
  ven22 = stringtovector(ven22)
  ear22 = stringtovector(ear22)
  mar22 = stringtovector(mar22)
  jup22 = stringtovector(jup22)
  sat22 = stringtovector(sat22)
  ura22 = stringtovector(ura22)
  nep22 = stringtovector(nep22)}

{SSMass = c(19885000, 3.302, 48.685, 59.721, 6.4171, 18981, 5683.4,
            868.13, 1024.1)
  SSMass = SSMass * (10^23)
  # The PV vectors are of the form x1 y1 z1 vx1 vy1 vz1 x2 y2
  # z2...
  SSPV23 = c(sun23, mer23, ven23, ear23, mar23, jup23, sat23, ura23,
             nep23)
  SSPV23 = SSPV23 * 1000 #Converts from km to m
  SSPV22 = c(sun22, mer22, ven22, ear22, mar22, jup22, sat22, ura22,
             nep22)
  SSPV22 = SSPV22 * 1000
  SSP22 = c(sun22[1:3], mer22[1:3], ven22[1:3], ear22[1:3], mar22[1:3],
            jup22[1:3], sat22[1:3], ura22[1:3], nep22[1:3])
  SSV22 = c(sun22[4:6], mer22[4:6], ven22[4:6], ear22[4:6], mar22[4:6],
            jup22[4:6], sat22[4:6], ura22[4:6], nep22[4:6])
  SSP22 = SSP22 * 1000
  SSV22 = SSV22 * 1000
  PlanetNames = c("Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter",
                  "Saturn", "Uranus", "Neptune")
  Year = 365 * 24 * 60 * 60 #Approximately one earth year in seconds
  G = 6.6743e-11 #Gravitational Constant
  ChargeConstant = 8987600000 #The proportionality coefficent in Coulomb's Law (1 / 4*pi*epsilon)
  ME = SSMass[4] #Mass of Earth
  RE = 6371000 #Radius of Earth
}
# Calculates the effect particle 2 has on particle 1
grav_2_on_1 = function(m2, p1, p2, constant) {
  r = p1- p2
  mag_r = norm(r, type = "2")
  unit_vec = r/mag_r
  g =-constant * m2 * unit_vec/(mag_r^2)
  return(g)
}

#TestofgravityonEarth'ssurface,(-9.81,0,0)is
#expected
grav_2_on_1(ME,c(RE,0,0),c(0,0,0),G)

#Takes the vector mass and positions of a set of particles
#and returns a vector of the resultant acceleration of
#each particles. We will mostly be using dim=3
#(3-dimensional) and constant=G.
gravity_of_all=function(m,p,dim,constant){
  l=length(m)
  gxyz=rep(0,l*dim)
  for (i in 1:l){
    m_without=m[-i]
    p_without=p[-(((dim*i)-(dim-1)):(dim*i))]
    p_cur=p[((dim*i)-(dim-1)):(dim*i)]
    for(j in 1:(l-1)){
      gxyz[((dim*i)-(dim-1)):(dim*i)]=gxyz[((dim*
                                               i)-(dim-1)):(dim *i)]+grav_2_on_1(m_without[j],
                                                                                 p_cur,p_without[((dim*j)-(dim-1)):(dim*
                                                                                                                      j)],constant)
    }
  }
  return(gxyz)
}

{
  testM=c(ME,70) #A test vector with the mass of the Earth and the mass of an average human
  testP=c(0,0,0,RE,0,0) #A test vector setting Earth at the origin and a human
  #on the surface of Earth
  #Earth should have a very low acceleration in the x-axis
  #and the human should have an acceleration of -g in the
  #x-axis
  gravity_of_all(testM,testP,3,G)
}

#Returns the acceleration of each planet due to the other
#planets at the initial frame of reference
gravity_of_all(SSMass,SSP22,3,G)

#These next 3 blocks of code use the respective
#approximation to output a matrix with each position in
#the first row and each velocity in the second, for
#example : x1 y1 z1 x2 y2 z2 vx1 vy1 vz1 vx2 vy2 vz2
Euler_Method= function(m,x,v,delta,constant,dim){
  g=gravity_of_all(m,x,dim,constant)
  x1=x+(v*delta)
  v1=v+(g*delta)
  return(matrix(c(x1,v1),nrow=2,byrow=TRUE))
}
Euler_Cromer_Method=function(m,x,v,delta,constant,dim){
  g=gravity_of_all(m,x,dim,constant)
  v1=v+(g*delta)
  x1=x+(v1*delta)
  return(matrix(c(x1,v1),nrow=2,byrow=TRUE))
}
Euler_Richardson_Method= function(m,x,v,delta,constant,dim){
  g=gravity_of_all(m,x,dim,constant)
  vmid=v+(0.5*g*delta)
  xmid=x+(0.5*v*delta)
  gmid=gravity_of_all(m,xmid, dim,constant)
  x1=x+(vmid*delta)
  v1=v+(gmid*delta)
  return(matrix(c(x1,v1),nrow=2,byrow=TRUE))
}


# This code takes an input of the initial conditions of the
# system, the timespan to simulate for and the
# approximation to use. It will return a data frame
# consisting of the position and velocity of each particle
# at each time interval
TimeEvolutionSystem = function(m, x0, v0, t, n, constant, process, names) {
  delta = t/n
  l = length(m)
  dim = length(x0)/l
  X = V = matrix(rep(0, dim * l * (n + 1)), nrow = (n + 1))
  X[1, ] = x0
  V[1, ] = v0
  for (i in 1:n) {
    xv = process(m, X[i, ], V[i, ], delta, constant, dim)
    X[(i + 1), ] = xv[1, ]
    V[(i + 1), ] = xv[2, ]
  }
  Xdata = as.data.frame(X[, (1:dim)])
  Vdata = as.data.frame(V[, (1:dim)])
  for (i in 1:dim) {
    names(Xdata)[i] = paste(names[1], "x", i)
    names(Vdata)[i] = paste(names[1], "v", i)
  }
  XV = cbind(Xdata, Vdata)
  for (j in 2:l) {
    Xdata = as.data.frame(X[, (((dim * j)- (dim- 1)):(dim *
                                                        j))])
    Vdata = as.data.frame(V[, (((dim * j)- (dim- 1)):(dim *
                                                        j))])
    for (i in 1:dim) {
      names(Xdata)[i] = paste(names[j], "x", i)
      names(Vdata)[i] = paste(names[j], "v", i)
    }
    XV = cbind(XV, Xdata, Vdata)
  }
  return(XV)
}
# A simple test of a human skydiving on Earth from 4000m.
# This will not be too accurate as the simulation cannot
# account for resistance or terminal velocity, and the
# human would also fall through the surface of the Earth.
# However, it should be expected that the human accelerates
# towards the Earth at-g, whilst the Earth barely moves.

TimeEvolutionSystem(c(ME, 60), c(0, 0, 0, RE + 4000, 0, 0), c(0, 0, 0, 0, 0, 0), 10, 10, G, Euler_Method, c("Earth", "Human"))
{
  E_S_Time=Sys.time()
  E_DF=TimeEvolutionSystem(SSMass,SSP22,SSV22,Year,2000,G,Euler_Method,PlanetNames)
  E_F_Time=Sys.time()
  E_F_Time-E_S_Time #Time taken to run the Euler method
}

{
  EC_S_Time=Sys.time()
  EC_DF=TimeEvolutionSystem(SSMass,SSP22,SSV22,Year,2000,G,Euler_Cromer_Method,PlanetNames)
  EC_F_Time=Sys.time()
  EC_F_Time-EC_S_Time #Time taken to run the Euler-Cromer method
}

{
  ER_S_Time=Sys.time()
  ER_DF=TimeEvolutionSystem(SSMass,SSP22,SSV22,Year,2000,G,Euler_Richardson_Method,PlanetNames)
  ER_F_Time=Sys.time()
  ER_F_Time-ER_S_Time #Time taken to run the Euler-Richardson method
}

{
  ER_half_S_Time=Sys.time()
  ER_DF=TimeEvolutionSystem(SSMass,SSP22,SSV22,Year,1000,G,Euler_Richardson_Method,PlanetNames)
  ER_half_F_Time=Sys.time()
  #Time taken to run the Euler-Richardson method with half
  #the intervals
  ER_half_F_Time-ER_half_S_Time
}

#This code can reduce a dataframe about an entire system
#to only include data about specific particles involved.
#This way, it is possible to evaluate data about specific
#particles whilst including the effect caused by other
#particles in the system
LimitParticlesDF= function(DF,dim,limit){
  veclimit=rep(0,(dim*2*length(limit)))
  for (i in 1:length(limit)){
    veclimit[((1+(i*2*dim)-(2*dim)):(i*2*dim))]=(2*
                                                   dim*limit[i])-seq((2*dim)-1,0,by=-1)
  }
  NewDF = subset(DF,select=veclimit)
  return(NewDF)
}
#These dataframes reduce the dataframe to just the inner
#Solar System (Sun to Mars)
{
  I_E_DF=LimitParticlesDF(E_DF,3,(1:5))
  I_EC_DF=LimitParticlesDF(EC_DF,3,(1:5))
  I_ER_DF=LimitParticlesDF(ER_DF,3,(1:5))
}
#This code reduces the dataframe to a select amount of
#evenly-spaced time intervals
LimitIterationsDF=function(DF,iterations){
  interval=(nrow(DF)-1)/iterations
  intseq=seq(1,by=interval,length=iterations+1)
  intseq=round(intseq)
  NewDF=DF[intseq,]
  return(NewDF)
}
#Here, the time intervals are reduced to being one month
#apart
{
  R_E_DF=LimitIterationsDF(E_DF,12)
  R_EC_DF=LimitIterationsDF(EC_DF,12)
  R_ER_DF=LimitIterationsDF(ER_DF,12)
}

#This code plots a dataframe with the desired
#coordinates. We will use dimplot=c(1,2) to plot the
#first coordinate (x) against the second (y)
DataFramePlot= function(DF,dim,dimplot,title,xlabel,ylabel){
  nparticle=ncol(DF)/(2*dim) #How many particles in the system
  maxx=minx=DF[1,dimplot[1]] #Helps to set the size of the graph
  maxy=miny=DF[1,dimplot[2]]
  particlenames="null"
  for (i in 1:nparticle){
    maxx=max(maxx,DF[,(i* 2 *dim)-((2*dim)-dimplot[1])])
    maxy=max(maxy,DF[,(i* 2 *dim)-((2*dim)-dimplot[2])])
    minx=min(minx,DF[,(i* 2 *dim)-((2*dim)-dimplot[1])])
    miny=min(miny,DF[,(i* 2 *dim)-((2*dim)-dimplot[2])])
    name=gsub(" x 1","",names(DF)[(i*2*dim)-((2*dim)-dimplot[1])])
    particlenames=c(particlenames,name)
  }
  particlenames=particlenames[-1]
  plot(DF[,dimplot[1]],DF[,dimplot[2]],type="b",xlim=c((1.25*minx),(1.25*maxx)),ylim=c((1.25*miny),(1.25*maxy)),main=title,xlab = xlabel,ylab=ylabel,asp=1,cex=0.2,cex.main=0.75)
  for (i in 2:nparticle){
    lines(DF[,(i*2*dim)-((2*dim)-dimplot[1])],
          DF[,(i*2*dim)-((2 *dim)-dimplot[2])],type="b",
          col=i,cex=0.2)
  }
  legend("bottomright",legend= particlenames,fill=(1:nparticle))
}

DataFramePlot(R_E_DF,3,c(1,2),"Simulation of the Solar System over an Earth year using the Euler method","X (m)","Y (m)")

DataFramePlot(I_E_DF,3,c(1,2),"Simulation of the Inner Solar System over an Earth year using the Euler method","X (m)","Y (m)")

DataFramePlot(R_EC_DF,3,c(1,2),"Simulation of the Solar System over an Earth year using the Euler-Cromer method","X (m)","Y (m)")

DataFramePlot(I_EC_DF,3,c(1,2),"Simulation of the Inner Solar System over an Earth year using the Euler-Cromer method","X (m)","Y (m)")

DataFramePlot(R_ER_DF,3,c(1,2),"Simulation of the Solar System over an Earth year using the Euler-Richardson method","X (m)","Y (m)")

DataFramePlot(I_ER_DF,3,c(1,2),"Simulation of the Inner Solar System over an Earth year using the Euler-Richardson method","X (m)","Y (m)")

# Checks the total linear momentum in each coordinate
L_Momentum_Conservation = function(DF, m, dim) {
  for (i in 1:dim) {
    sequence = seq(((2 * dim)- (dim- i)), by = (2 * dim), length = length(m))
    for (j in 1:nrow(DF)) {
      DF[j, i] = sum(DF[j, sequence] * m)
    }
    names(DF)[i] = paste("Momentum X", i)
  }
  DF = DF[1:dim]
  return(DF)
}

# Total Kinetic Energy of the System
Kinetic_Energy = function(DF, m, dim) {
  K = rep(0, nrow(DF))
  nparticles = ncol(DF)/(2 * dim)
  V2 = rep(0, nparticles)
  for (i in 1:nrow(DF)) {
    for (j in 1:nparticles) {
      V2[j] = sum((DF[i, (((2 * dim * j) - (dim - 1)):(2 * dim * j))])^2)
    }
    K[i] = sum(0.5 * m * V2)
  }
  return(K)
}

# Total Potential Energy of the System
Potential_Energy = function(DF, m, dim, constant) {
  U = rep(0, nrow(DF))
  nparticles = ncol(DF)/(2 * dim)
  for (i in 1:nrow(DF)) {
    Ucur = rep(0, dim)
    for (j in (1:nparticles)[-nparticles]) {
      for (k in ((j + 1):nparticles)) {
        r = (DF[i, (((2 * dim * j)- ((2 * dim)- 1)):((2 *dim * j)- (2 * dim)))]- DF[i, (((2 * dim * k)- ((2 * dim)- 1)):((2 * dim * k)- (2 * dim)))])
        Ucur = Ucur + 2 * constant * m[j] * m[k]/sum(abs(r))
      }
    }
    U[i] = sum(Ucur)
  }
  return(U)
}

Total_Energy = function(DF, m, dim, constant) {
  K = Kinetic_Energy(DF, m, dim)
  U = Potential_Energy(DF, m, dim, constant)
  Total = K + U
  return(Total)
}

{
  L_Momentum_Conservation(R_E_DF, SSMass, 3) #Total linear momentum of E method
  
  L_Momentum_Conservation(R_ER_DF, SSMass, 3) #Total linear momentum of EC method
  
  L_Momentum_Conservation(R_ER_DF, SSMass, 3) #Total linear momentum of ER method
  
  
  Total_Energy(R_E_DF,SSMass,3,G) #Evolution of total energy in Euler method
  
  Total_Energy(R_EC_DF,SSMass,3,G) #Evolution of total energy in EC method
  
  Total_Energy(R_ER_DF,SSMass,3,G) #Evolution of total energy in ER method
}

#Finds the relative error between the experimental values
#DF1 and the accepted values DF2
relative_error=function(DF1,DF2){
  final_result=DF1[nrow(DF1),]
  error=(final_result-DF2)/DF2
  return(error)
}

#Finds the root mean square percent error for each
#particle
RMSPE_particles= function(error,dim){
  nparticles=length(error)/(dim*2)
  mat= matrix(rep(0,nparticles*2),nrow=1)
  ans= as.data.frame(mat)
  for (i in 1:(nparticles*2)){
    ans[,i]=RMSPE_system(error[(((dim*i)-(dim-1)):(dim*
                                                     i))])
    name=names(error)[((dim*i)-(dim-1))]
    name=gsub("1","",name)
    names(ans)[i]=name
  }
  return(ans)
}

#Finds the RMSPE of the whole system
RMSPE_system= function(error){
  error = error^2
  sum_error=sum(error)
  mean = sum_error/length(error)
  root_mean = 100 * sqrt(mean)
  return(root_mean)
}

E_RE = relative_error(E_DF, SSPV23) #Relative error in the Euler Method
RMSPE_system(E_RE) #Root mean square percentage error in the Euler Method

RMSPE_particles(E_RE, 3) #RMSPE for each particle in the Euler Method

EC_RE = relative_error(EC_DF, SSPV23) #Relative error in the Euler-Cromer Method
RMSPE_system(EC_RE) #RMSPE in the Euler-Cromer Method

RMSPE_particles(EC_RE, 3) #RMSPE for each particle in the Euler-Cromer Method

ER_RE = relative_error(ER_DF, SSPV23) #Relative error in the Euler-Richardson Method
RMSPE_system(ER_RE) #RMSPE in the Euler-Richardson Method

RMSPE_particles(ER_RE, 3) #RMSPE for each particle in the Euler-Richardson Method

100 * ER_RE[24]

