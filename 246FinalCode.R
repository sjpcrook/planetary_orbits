#Constants
{
#2460066.500000000 = A.D. 2023-May-02 00:00:00.0000 TDB 
#       Mass               X-Position             Y-Position           Z-Position             X-Velocity                Y-Velocity              Z-Velocity

sun = c(1988500E+24, -1.324316846741195E+06, -1.446765512621427E+05, 3.202738040537933E+04, 4.021534124461124E-03, -1.472929900704240E-02, 2.948517133131729E-05)
mer = c(3.302E+23, -5.140210130150759E+07, -4.381189053847280E+07, 1.056887671139676E+06, 2.207720296247787E+01, -3.453839718779295E+01, -4.845914556700924E+00)
ven = c(48.685E+23, -9.848082136897813E+07, 4.552713536799312E+07, 6.265256376415793E+06, -1.504882769689358E+01, -3.187376111911025E+01, 4.312134646762669E-01)
ear = c(5.9722E+24, -1.150223276469004E+08, -9.909841274895248E+07, 3.755860068777204E+04, 1.907614272241964E+01, -2.258721645907844E+01, 2.018298546438757E-03)
mar = c(6.4171E+23, -2.010819716042931E+08, 1.476265433762218E+08, 8.029003272268474E+06, -1.349413909238045E+01, -1.742820367074673E+01, -3.381872148509668E-02)
jup = c(18.982E+26, 6.792848124845287E+08, 2.926757289590889E+08, -1.641170746619791E+07, -5.318068089405075E+00, 1.261385267190711E+01, 6.658478320447347E-02)
sat = c(5.6834E+26, 1.265581460039990E+09, -7.382676489447224E+08, -3.755222856885275E+07, 4.325577300311080E+00, 8.324810463771845E+00, -3.167483277652900E-01)
ura = c(86.813E+24, 1.945406877724633E+09, 2.202709703800068E+09, -1.702215196928215E+07, -5.154055480935001E+00, 4.190735053981157E+00, 8.237035843708407E-02)
nep = c(102.409E+24, 4.455662623587698E+09, -3.833067895772491E+08, -9.479194396546717E+07, 4.297305233379821E-01, 5.447227229797033E+00, -1.222099521045787E-01)
} #Data taken from Horizons System
{
SSMass = c(sun[1], mer[1], ven[1], ear[1], mar[1], jup[1], sat[1], ura[1], nep[1])
SSPosition = c(sun[2:4], mer[2:4], ven[2:4], ear[2:4], mar[2:4], jup[2:4], sat[2:4], ura[2:4], nep[2:4])
SSVelocity = c(sun[5:7], mer[5:7], ven[5:7], ear[5:7], mar[5:7], jup[5:7], sat[5:7], ura[5:7], nep[5:7])
SSPosition = SSPosition*1000 # Convert from km to m
VelocityVectorSolar = VelocityVector*1000 # Convert from km/s to m/s
} #Code to convert solar system data into vectors for use
{
  InnerSSMass = SSMass[1:5]
  InnerSSPosition = SSPosition[1:15]
  InnerSSVelocity = SSVelocity[1:15]
} #Vectors only including the inner Solar System (Sun - Mars)
Year = 365*24*60*60 #Approximately one earth year in seconds

G = 6.6743E-11 #Gravitational Constant 
ChargeConstant = 8.9876E+9 #The proportionality coefficent in Coulomb's Law (1 / 4*pi*epsilon)

ME = SSMass[4]
RE = 6371000

#Gravity Calculation Systems

grav_2_on_1 = function(m2, p1, p2, constant) # Calculates the effect particle 2 has on particle 1
{
  r = p1 - p2
  mag_r = norm(r, type = "2")
  unit_vec = r / mag_r
  g = -constant*m2*unit_vec/(mag_r^2)
  g[is.na(g)] = 0
  return(g)
}
grav_2_on_1(MassVector[4], c(6371000, 0, 0), c(0,0,0), G)
grav_2_on_1(MassVector[4], c(6371000, -6371000, 10000000), c(0,0,0), G)

gravity_of_all = function(m, p, dim, constant) # Takes the vector mass and positions of a set of particles and returns a vector of the resultant effect on each particles
{
  l = length(m)
  gxyz = rep(0, l*dim)
  for(i in 1:l)
  {
    m_without = m[-i]
    p_without = p[-(((dim*i)-(dim-1)):(dim*i))]
    p_cur = p[((dim*i)-(dim-1)):(dim*i)]
    for(j in 1:(l-1))
    {
      gxyz[((dim*i)-(dim-1)):(dim*i)] = gxyz[((dim*i)-(dim-1)):(dim*i)] + grav_2_on_1(m_without[j], p_cur, p_without[((dim*j)-(dim-1)):(dim*j)], constant)
    }
  }
  return(gxyz)
}

testM = c(MassVector[4], 70)
testP = c(0, 0, 0, 6371000, 0, 0)
gravity_of_all(testM, testP, 3, G)
#gravity_of_all(test_mom_m, test_mom_p, 3, G)
gravity_of_all(MassVector, SSPosition, 3, G)

#Systems
SystemEuler = function(m, p0, v0, t, n, constant) #x[n+1] = x[n] + v[n]t, v[n+1]  = v[n] + a[n]t
{
  delta = t/n
  l = length(m)
  dim = length(p0) / l
  A = matrix(rep(0, dim*l*(n+1)), nrow = (n+1))
  A[1, ] = p0
  for(i in 1:n)
  {
    gxyz = gravity_of_all(m, A[i, ], dim, constant)
    A[(i+1), ] = A[i, ] + (v0*delta)
    v0 = v0 + (gxyz*delta)
  }
  return(A)
}
A = SystemEuler(MassVector, SSPosition, VelocityVector, Year, 365, G)
AMars = SystemEuler(MassVector[1:5], SSPosition[1:15], VelocityVector[1:15], Year, 10000, G)

SystemEulerCromer = function(m, p0, v0, t, n, constant) #x[n+1] = x[n] + v[n+1]t, v[n+1]  = v[n] + a[n]t
{
  delta = t/n
  l = length(m)
  dim = length(p0) / l
  A = matrix(rep(0, dim*l*(n+1)), nrow = (n+1))
  A[1, ] = p0
  for(i in 1:n)
  {
    gxyz = gravity_of_all(m, A[i, ], dim, constant)
    v0 = v0 + (gxyz*delta)
    A[(i+1), ] = A[i, ] + (v0*delta)
  }
  return(A)
}
B = SystemEulerCromer(MassVector, SSPosition, VelocityVector, Year, 365, G)
BMars = SystemEulerCromer(MassVector[1:5], SSPosition[1:15], VelocityVector[1:15], Year, 1000, G)

SystemEulerRichardson = function(m, p0, v0, t, n, constant) #x[n+1] = x[n] + v[mid]t, v[n+1]  = v[n] + a[mid]t
{
  delta = t/n
  l = length(m)
  dim = length(p0) / l
  A = matrix(rep(0, dim*l*(n+1)), nrow = (n+1))
  A[1, ] = p0
  for(i in 1:n)
  {
    gxyz = gravity_of_all(m, A[i, ], dim, constant)
    vmid = v0 + (0.5*gxyz*delta)
    pmid = A[i, ] + (0.5*v0*delta)
    amid = gravity_of_all(m, pmid, dim, constant)
    A[(i+1), ] = A[i, ] + (vmid*delta)
    v0 = v0 + (amid*delta)
  }
  return(A)
}
C = SystemEulerRichardson(MassVector, SSPosition, VelocityVector, Year, 365, G)
CMars = SystemEulerRichardson(MassVector[1:5], PositionVector[1:15], VelocityVector[1:15], Year, 1000, G)

#Plotting
MatrixPlot = function(A, n, dimplot)
{
  nx = length(A[, 1])
  dim = length(A[1, ])/n
  X = rep(0, nx*n)
  Y = X
  for(i in 1:n)
  {
    X[(((i*nx)-(nx-1)):(i*nx))] = A[, ((i*dim)-(dim-dimplot[1]))]
    Y[(((i*nx)-(nx-1)):(i*nx))] = A[, ((i*dim)-(dim-dimplot[2]))]
  }
  maxx = max(abs(X))
  maxy = max(abs(Y))
  plot(X[1:nx], Y[1:nx], type="b", xlim = c(-maxx, maxx), ylim = c(-maxy, maxy))
  for(i in 2:n)
  {
    lines(X[(((i*nx)-(nx-1)):(i*nx))], Y[(((i*nx)-(nx-1)):(i*nx))], type = "b", col = i)
  }
}

MatrixPlot(A, 9, c(1, 2))
MatrixPlot(AMars, 5, c(1, 2))

MatrixPlot(B, 9, c(1, 2))
MatrixPlot(BMars, 5, c(1, 2))

MatrixPlot(C, 9, c(1, 2))
MatrixPlot(CMars, 5, c(1, 2))
