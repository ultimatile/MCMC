using ProgressMeter

#classical J1-J2 Heisenberg model
#Heatbath method
#Annealing

abstract type ModelParameters end
struct mp <: ModelParameters
  L::Int64
  J1::Float64 # neareast neighbor exchange coupling
  J2::Float64 # next-neareast neighbor exchange coupling
  runs::Int64
  Tmin::Float64
  Tmax::Float64
  Tsteps::Int64
  deltaT::Float64
  mcs_max::Int64 #Maximum number of MC steps
  discard::Int64 #Thermalization
  frac::Int64
  inext::Array{Int64,1}
  iprev::Array{Int64,1}

  function mp(L::Int64, J1::Float64, J2::Float64, runs::Int64, Tmin::Float64, Tmax::Float64, Tsteps::Int64, mcs_max::Int64, discard::Int64)
    deltaT = (Tmax - Tmin) / (Tsteps - 1)
    inext, iprev = neighbor_indices(L)
    frac = mcs_max - discard
    new(L, J1, J2, runs, Tmin, Tmax, Tsteps, deltaT, mcs_max, discard, frac, inext, iprev)
  end
end

sintheta(c) = sqrt(1 - clamp(c, -1.0, 1.0)^2)

SiSj(pi, pj, ci, cj) = cos(pi - pj) * sintheta(ci) * sintheta(cj) + ci * cj

function neighbor_indices(L::Int)
  inext = [i < L ? i + 1 : 1 for i in 1:L]
  iprev = [i > 1 ? i - 1 : L for i in 1:L]
  return inext, iprev
end

function compute_hlocal(phi, costheta, ix, iy, inextx, inexty, iprevx, iprevy, J1, J2)
  hlocalx = 0.0
  hlocaly = 0.0
  hlocalz = 0.0
  J1isnotzero = !isapprox(J1, 0.0; atol=eps(Float64))
  J2isnotzero = !isapprox(J2, 0.0; atol=eps(Float64))
  if J1isnotzero
    sinthetapx = sintheta(costheta[inextx, iy])
    sinthetamx = sintheta(costheta[iprevx, iy])
    sinthetapy = sintheta(costheta[ix, inexty])
    sinthetamy = sintheta(costheta[ix, iprevy])
    hlocalx += J1 * (cos(phi[inextx, iy]) * sinthetapx + cos(phi[iprevx, iy]) * sinthetamx
                     + cos(phi[ix, inexty]) * sinthetapy + cos(phi[ix, iprevy]) * sinthetamy)
    hlocaly += J1 * (sin(phi[inextx, iy]) * sinthetapx + sin(phi[iprevx, iy]) * sinthetamx
                     + sin(phi[ix, inexty]) * sinthetapy + sin(phi[ix, iprevy]) * sinthetamy)
    hlocalz += J1 * (costheta[inextx, iy] + costheta[iprevx, iy] + costheta[ix, inexty] + costheta[ix, iprevy])
  end
  if J2isnotzero
    sinthetapxpy = sintheta(costheta[inextx, inexty])
    sinthetapxmy = sintheta(costheta[inextx, iprevy])
    sinthetamxpy = sintheta(costheta[iprevx, inexty])
    sinthetamxmy = sintheta(costheta[iprevx, iprevy])
    hlocalx += J2 * (cos(phi[inextx, inexty]) * sinthetapxpy + cos(phi[inextx, iprevy]) * sinthetapxmy
                     + cos(phi[iprevx, inexty]) * sinthetamxpy + cos(phi[iprevx, iprevy]) * sinthetamxmy)
    hlocaly += J2 * (sin(phi[inextx, inexty]) * sinthetapxpy + sin(phi[inextx, iprevy]) * sinthetapxmy
                     + sin(phi[iprevx, inexty]) * sinthetamxpy + sin(phi[iprevx, iprevy]) * sinthetamxmy)
    hlocalz += J2 * (costheta[inextx, inexty] + costheta[inextx, iprevy] + costheta[iprevx, inexty] + costheta[iprevx, iprevy])
  end
  return hlocalx, hlocaly, hlocalz
end

function observe_ene(mp::ModelParameters, phi, costheta)
  tmpE = 0.0
  inext = mp.inext
  iprev = mp.iprev
  L = mp.L
  J1 = mp.J1
  J2 = mp.J2
  J1isnotzero = !iszero(J1)
  J2isnotzero = !iszero(J2)
  for ix in 1:L, iy in 1:L
    phixy = phi[ix, iy]
    costhetaxy = costheta[ix, iy]
    ixnext = inext[ix]
    iynext = inext[iy]
    ixprev = iprev[ix]
    iyprev = iprev[iy]
    if J1isnotzero
      tmpE += J1 * (
        SiSj(phixy, phi[ix, iynext], costhetaxy, costheta[ix, iynext])
        + SiSj(phixy, phi[ix, iyprev], costhetaxy, costheta[ix, iyprev])
        + SiSj(phixy, phi[ixnext, iy], costhetaxy, costheta[ixnext, iy])
        + SiSj(phixy, phi[ixprev, iy], costhetaxy, costheta[ixprev, iy])
      )
    end
    if J2isnotzero
      tmpE += J2 * (
        SiSj(phixy, phi[ixnext, iynext], costhetaxy, costheta[ixnext, iynext])
        + SiSj(phixy, phi[ixnext, iyprev], costhetaxy, costheta[ixnext, iyprev])
        + SiSj(phixy, phi[ixprev, iynext], costhetaxy, costheta[ixprev, iynext])
        + SiSj(phixy, phi[ixprev, iyprev], costhetaxy, costheta[ixprev, iyprev])
      )
    end
  end
  return tmpE
end

function update_HB!(mp::ModelParameters, phi, costheta, T)
  L = mp.L
  inext = mp.inext
  iprev = mp.iprev
  for ix in 1:L, iy in 1:L
    hlocalx, hlocaly, hlocalz = compute_hlocal(phi, costheta, ix, iy, inext[ix], inext[iy], iprev[ix], iprev[iy], mp.J1, mp.J2)
    hlocalxy = hlocalx^2 + hlocaly^2
    sqrt_hlocalxy = sqrt(hlocalxy)
    sqrt_hlocalxyz = sqrt(hlocalxy + hlocalz^2)
    betaH = sqrt_hlocalxyz / T
    sznew = -1 - log1p(rand() * expm1(-2 * betaH)) / betaH
    sinthetanew = sintheta(sznew)
    cpsi = hlocalz / sqrt_hlocalxyz
    spsi = sintheta(cpsi)
    cphi = hlocalx / sqrt_hlocalxy
    sphi = hlocaly / sqrt_hlocalxy
    phi_rand = 2pi * rand()
    sxnew = cos(phi_rand) * sinthetanew
    synew = sin(phi_rand) * sinthetanew
    sx = cpsi * cphi * sxnew - sphi * synew + spsi * cphi * sznew
    sy = cpsi * sphi * sxnew + cphi * synew + spsi * sphi * sznew
    costheta[ix, iy] = -spsi * sxnew + cpsi * sznew
    # mod is necessary to transform the angle from [-pi,pi] to [0,2pi]
    phi[ix, iy] = mod(atan(sy, sx), 2pi)
  end
end

function main(mp::ModelParameters)
  Tsteps = mp.Tsteps
  L = mp.L
  runs = mp.runs
  deltaT = mp.deltaT
  frac = mp.frac
  discard = mp.discard
  ave_ene = zeros(Tsteps)
  var_ene = zeros(Tsteps)
  ave_spec = zeros(Tsteps)
  var_spec = zeros(Tsteps)
  Ts = [mp.Tmax - deltaT * (Tstep - 1) for Tstep in 1:Tsteps]
  @showprogress for run in 1:runs
    phi = 2pi * rand(L, L)
    costheta = 2rand(L, L) .- 1
    #phi=zeros(L,L);costheta=zeros(L,L)
    for Tstep in 1:Tsteps
      energy = 0.0
      energy2 = 0.0
      T = Ts[Tstep]
      for mcs in 1:mp.mcs_max
        update_HB!(mp, phi, costheta, T)
        if mcs > discard
          ene_per_mcs = observe_ene(mp, phi, costheta)
          energy += ene_per_mcs
          energy2 += ene_per_mcs^2
        end
      end
      energy /= frac * 2
      energy2 /= frac * 4
      spec = (energy2 - energy^2) / T^2 / L^2
      ave_spec[Tstep] += spec
      var_spec[Tstep] += spec^2
      ave_ene[Tstep] += energy
      var_ene[Tstep] += energy^2
    end
  end
  ave_spec /= runs
  var_spec /= runs
  ave_ene /= runs
  var_ene /= runs
  open("E_HB.dat", "w") do fp_E
    open("C_HB.dat", "w") do fp_C
      for Tstep in 1:Tsteps
        T = Ts[Tstep]
        write(fp_E, "$T $(ave_ene[Tstep]) $(sqrt((var_ene[Tstep]-ave_ene[Tstep]^2)/(runs-1)))\n")
        write(fp_C, "$T $(ave_spec[Tstep]) $(sqrt((var_spec[Tstep]-ave_spec[Tstep]^2)/(runs-1)))\n")
      end
    end
  end
end


function main(ARGS)
  L = 4
  J1 = 1.0
  J2 = 0.0
  runs = 5
  Tmin = 0.1
  Tmax = 2.0
  Tsteps = 50
  #Tmin=2;Tmax=Tmin;Tsteps=1
  mcs_max = 20000
  discard = 10000
  return main(mp(L, J1, J2, runs, Tmin, Tmax, Tsteps, mcs_max, discard))
end

@main
