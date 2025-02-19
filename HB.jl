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
  ip::Array{Int64,1}
  im::Array{Int64,1}

  function mp(L::Int64, J1::Float64, J2::Float64, runs::Int64, Tmin::Float64, Tmax::Float64, Tsteps::Int64, mcs_max::Int64, discard::Int64)
    deltaT = (Tmax - Tmin) / (Tsteps - 1)
    ip, im = neighbor_indices(L)
    frac = mcs_max - discard
    new(L, J1, J2, runs, Tmin, Tmax, Tsteps, deltaT, mcs_max, discard, frac, ip, im)
  end
end

sintheta(c) = sqrt(1 - clamp(c, -1.0, 1.0)^2)

SiSj(pi, pj, ci, cj) = cos(pi - pj) * sintheta(ci) * sintheta(cj) + ci * cj

function neighbor_indices(L::Int)
  ip = [i < L ? i + 1 : 1 for i in 1:L]
  im = [i > 1 ? i - 1 : L for i in 1:L]
  return ip, im
end

function compute_hlocal(phi, costheta, ix, iy, ipx, ipy, imx, imy, J1, J2)
  hlocalx = 0.0
  hlocaly = 0.0
  hlocalz = 0.0
  J1isnotzero = !isapprox(J1, 0.0; atol=eps(Float64))
  J2isnotzero = !isapprox(J2, 0.0; atol=eps(Float64))
  if J1isnotzero
    sinthetapx = sintheta(costheta[ipx, iy])
    sinthetamx = sintheta(costheta[imx, iy])
    sinthetapy = sintheta(costheta[ix, ipy])
    sinthetamy = sintheta(costheta[ix, imy])
    hlocalx += J1 * (cos(phi[ipx, iy]) * sinthetapx + cos(phi[imx, iy]) * sinthetamx
                     + cos(phi[ix, ipy]) * sinthetapy + cos(phi[ix, imy]) * sinthetamy)
    hlocaly += J1 * (sin(phi[ipx, iy]) * sinthetapx + sin(phi[imx, iy]) * sinthetamx
                     + sin(phi[ix, ipy]) * sinthetapy + sin(phi[ix, imy]) * sinthetamy)
    hlocalz += J1 * (costheta[ipx, iy] + costheta[imx, iy] + costheta[ix, ipy] + costheta[ix, imy])
  end
  if J2isnotzero
    sinthetapxpy = sintheta(costheta[ipx, ipy])
    sinthetapxmy = sintheta(costheta[ipx, imy])
    sinthetamxpy = sintheta(costheta[imx, ipy])
    sinthetamxmy = sintheta(costheta[imx, imy])
    hlocalx += J2 * (cos(phi[ipx, ipy]) * sinthetapxpy + cos(phi[ipx, imy]) * sinthetapxmy
                     + cos(phi[imx, ipy]) * sinthetamxpy + cos(phi[imx, imy]) * sinthetamxmy)
    hlocaly += J2 * (sin(phi[ipx, ipy]) * sinthetapxpy + sin(phi[ipx, imy]) * sinthetapxmy
                     + sin(phi[imx, ipy]) * sinthetamxpy + sin(phi[imx, imy]) * sinthetamxmy)
    hlocalz += J2 * (costheta[ipx, ipy] + costheta[ipx, imy] + costheta[imx, ipy] + costheta[imx, imy])
  end
  return hlocalx, hlocaly, hlocalz
end

function calc_ene(mp::ModelParameters, phi, costheta)
  tmpE = 0.0
  ip = mp.ip
  im = mp.im
  L = mp.L
  J1 = mp.J1
  J2 = mp.J2
  J1isnotzero = !isapprox(J1, 0.0; atol=eps(Float64))
  J2isnotzero = !isapprox(J2, 0.0; atol=eps(Float64))
  for ix in 1:L, iy in 1:L
    if J1isnotzero
      tmpE += J1 * (
        SiSj(phi[ix, iy], phi[ix, ip[iy]], costheta[ix, iy], costheta[ix, ip[iy]])
        + SiSj(phi[ix, iy], phi[ix, im[iy]], costheta[ix, iy], costheta[ix, im[iy]])
        + SiSj(phi[ix, iy], phi[ip[ix], iy], costheta[ix, iy], costheta[ip[ix], iy])
        + SiSj(phi[ix, iy], phi[im[ix], iy], costheta[ix, iy], costheta[im[ix], iy])
      )
    end
    if J2isnotzero
      tmpE += J2 * (
        SiSj(phi[ix, iy], phi[ip[ix], ip[iy]], costheta[ix, iy], costheta[ip[ix], ip[iy]])
        + SiSj(phi[ix, iy], phi[ip[ix], im[iy]], costheta[ix, iy], costheta[ip[ix], im[iy]])
        + SiSj(phi[ix, iy], phi[im[ix], ip[iy]], costheta[ix, iy], costheta[im[ix], ip[iy]])
        + SiSj(phi[ix, iy], phi[im[ix], im[iy]], costheta[ix, iy], costheta[im[ix], im[iy]])
      )
    end
  end
  return tmpE
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
  ip = mp.ip
  im = mp.im
  for run in 1:runs
    phi = 2pi * rand(L, L)
    costheta = 2 * rand(L, L) .- 1
    #phi=zeros(L,L);costheta=zeros(L,L)
    for Tstep in 1:Tsteps
      energy = 0.0
      energy2 = 0.0
      T = Ts[Tstep]
      for mcs in 1:mp.mcs_max
        for ix in 1:L, iy in 1:L
          ipx = ip[ix]
          imx = im[ix]
          ipy = ip[iy]
          imy = im[iy]
          hlocalx, hlocaly, hlocalz = compute_hlocal(phi, costheta, ix, iy, ipx, ipy, imx, imy, mp.J1, mp.J2)
          betaH = sqrt(hlocalx^2 + hlocaly^2 + hlocalz^2) / T
          sznew = -1 - log1p(rand() * expm1(-2 * betaH)) / betaH
          sinthetanew = sintheta(sznew)
          cpsi = hlocalz / sqrt(hlocalx^2 + hlocaly^2 + hlocalz^2)
          spsi = sintheta(cpsi)
          cphi = hlocalx / sqrt(hlocalx^2 + hlocaly^2)
          sphi = hlocaly / sqrt(hlocalx^2 + hlocaly^2)
          phi_rand = 2pi * rand()
          sxnew = cos(phi_rand) * sinthetanew
          synew = sin(phi_rand) * sinthetanew
          sx = cpsi * cphi * sxnew - sphi * synew + spsi * cphi * sznew
          sy = cpsi * sphi * sxnew + cphi * synew + spsi * sphi * sznew
          costheta[ix, iy] = -spsi * sxnew + cpsi * sznew
          # mod is necessary to transform the angle from [-pi,pi] to [0,2pi]
          phi[ix, iy] = mod(atan(sy, sx), 2pi)
        end
        if mcs > discard
          tmpE = calc_ene(mp, phi, costheta)
          energy += tmpE
          energy2 += tmpE^2
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
        write(fp_E, "$(mp.Tmax-deltaT*(Tstep-1)) $(ave_ene[Tstep]) $(sqrt((var_ene[Tstep]-ave_ene[Tstep]^2)/(runs-1)))\n")
        write(fp_C, "$(mp.Tmax-deltaT*(Tstep-1)) $(ave_spec[Tstep]) $(sqrt((var_spec[Tstep]-ave_spec[Tstep]^2)/(runs-1)))\n")
      end
    end
  end
end

const L = 4
const J1 = 1.0
const J2 = 0.0
const runs = 5
const Tmin = 0.1
const Tmax = 2.0
const Tsteps = 50
#Tmin=2;Tmax=Tmin;Tsteps=1
const mcs_max = 20000
const discard = 10000

modpara = mp(L, J1, J2, runs, Tmin, Tmax, Tsteps, mcs_max, discard)
@time main(modpara)
