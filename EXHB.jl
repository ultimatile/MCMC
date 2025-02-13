#classical J1-J2 Heisenberg model
#Heatbath method
#Annealing
struct mp
  L::Int64
  J1::Float64 #exchange coupling
  J2::Float64 #hopping
  runs::Int64
  Tmin::Float64
  Tmax::Float64
  Tsteps::Int64
  mcs_max::Int64 #Maximum number of MC steps
  discard::Int64 #Thermalization
  frac::Int64
  ip::Array{Int64,1}
  im::Array{Int64,1}
end

sintheta(c) = sqrt(abs(1 - c^2))
SiSj(pi, pj, ti, tj) = cos(pi - pj) * sintheta(ti) * sintheta(tj) + ti * tj

function update_HB!(mp, phi, costheta, T)
  for ix in 1:mp.L
    for iy in 1:mp.L
      sinthetapx = sintheta(costheta[mp.ip[ix], iy])
      sinthetamx = sintheta(costheta[mp.im[ix], iy])
      sinthetapy = sintheta(costheta[ix, mp.ip[iy]])
      sinthetamy = sintheta(costheta[ix, mp.im[iy]])
      sinthetapxpy = sintheta(costheta[mp.ip[ix], mp.ip[iy]])
      sinthetapxmy = sintheta(costheta[mp.ip[ix], mp.im[iy]])
      sinthetamxpy = sintheta(costheta[mp.im[ix], mp.ip[iy]])
      sinthetamxmy = sintheta(costheta[mp.im[ix], mp.im[iy]])
      hlocalx = (mp.J1 * (cos(phi[mp.ip[ix], iy]) * sinthetapx + cos(phi[mp.im[ix], iy]) * sinthetamx
                          + cos(phi[ix, mp.ip[iy]]) * sinthetapy + cos(phi[ix, mp.im[iy]]) * sinthetamy)
                 +
                 mp.J2 * (cos(phi[mp.ip[ix], mp.ip[iy]]) * sinthetapxpy + cos(phi[mp.ip[ix], mp.im[iy]]) * sinthetapxmy
                          + cos(phi[mp.im[ix], mp.ip[iy]]) * sinthetamxpy + cos(phi[mp.im[ix], mp.im[iy]]) * sinthetamxmy))
      hlocaly = (mp.J1 * (sin(phi[mp.ip[ix], iy]) * sinthetapx + sin(phi[mp.im[ix], iy]) * sinthetamx
                          + sin(phi[ix, mp.ip[iy]]) * sinthetapy + sin(phi[ix, mp.im[iy]]) * sinthetamy)
                 +
                 mp.J2 * (sin(phi[mp.ip[ix], mp.ip[iy]]) * sinthetapxpy + sin(phi[mp.ip[ix], mp.im[iy]]) * sinthetapxmy
                          + sin(phi[mp.im[ix], mp.ip[iy]]) * sinthetamxpy + sin(phi[mp.im[ix], mp.im[iy]]) * sinthetamxmy))
      hlocalz = (mp.J1 * (costheta[mp.ip[ix], iy] + costheta[mp.im[ix], iy] + costheta[ix, mp.ip[iy]] + costheta[ix, mp.im[iy]])
                 +
                 mp.J2 * (costheta[mp.ip[ix], mp.ip[iy]] + costheta[mp.ip[ix], mp.im[iy]] + costheta[mp.im[ix], mp.ip[iy]] + costheta[mp.im[ix], mp.im[iy]]))
      betaH = sqrt(hlocalx^2 + hlocaly^2 + hlocalz^2) / T
      rndm = rand()
      sznew = -1 - log1p(rndm * expm1(-2 * betaH)) / betaH
      sinthetanew = sintheta(sznew)
      cpsi = hlocalz / sqrt(hlocalx^2 + hlocaly^2 + hlocalz^2)
      spsi = sintheta(cpsi)
      cphi = hlocalx / sqrt(hlocalx^2 + hlocaly^2)
      sphi = hlocaly / sqrt(hlocalx^2 + hlocaly^2)
      rndm = rand()
      sxnew = cos(2pi * rndm) * sinthetanew
      synew = sin(2pi * rndm) * sinthetanew
      sx = cpsi * cphi * sxnew - sphi * synew + spsi * cphi * sznew
      sy = cpsi * sphi * sxnew + cphi * synew + spsi * sphi * sznew
      costheta[ix, iy] = -spsi * sxnew + cpsi * sznew
      sinthetanew = sintheta(costheta[ix, iy])
      if sy / sinthetanew > 1
        phi[ix, iy] = 0.5pi
      elseif sx / sinthetanew > 1
        phi[ix, iy] = 0
      elseif sx / sinthetanew < -1
        phi[ix, iy] = pi
      elseif sy / sinthetanew < -1
        phi[ix, iy] = 1.5pi
      elseif sy / sinthetanew > 0
        phi[ix, iy] = acos(sx / sinthetanew)
      else
        phi[ix, iy] = 2pi - acos(sx / sinthetanew)
      end
    end
  end
  return phi, costheta
end

function calc_ene(mp, phi, costheta)
  tmpE = 0
  for ix in 1:mp.L
    for iy in 1:mp.L
      tmpE += mp.J1 * (
        SiSj(phi[ix, iy], phi[ix, mp.ip[iy]], costheta[ix, iy], costheta[ix, mp.ip[iy]])
        + SiSj(phi[ix, iy], phi[ix, mp.im[iy]], costheta[ix, iy], costheta[ix, mp.im[iy]])
        + SiSj(phi[ix, iy], phi[mp.ip[ix], iy], costheta[ix, iy], costheta[mp.ip[ix], iy])
        + SiSj(phi[ix, iy], phi[mp.im[ix], iy], costheta[ix, iy], costheta[mp.im[ix], iy])
      ) + mp.J2 * (
        SiSj(phi[ix, iy], phi[mp.ip[ix], mp.ip[iy]], costheta[ix, iy], costheta[mp.ip[ix], mp.ip[iy]])
        + SiSj(phi[ix, iy], phi[mp.ip[ix], mp.im[iy]], costheta[ix, iy], costheta[mp.ip[ix], mp.im[iy]])
        + SiSj(phi[ix, iy], phi[mp.im[ix], mp.ip[iy]], costheta[ix, iy], costheta[mp.im[ix], mp.ip[iy]])
        + SiSj(phi[ix, iy], phi[mp.im[ix], mp.im[iy]], costheta[ix, iy], costheta[mp.im[ix], mp.im[iy]])
      )
    end
  end
  return tmpE
end

function main(mp)
  ave_ene = zeros(mp.Tsteps)
  var_ene = zeros(mp.Tsteps)
  ave_spec = zeros(mp.Tsteps)
  var_spec = zeros(mp.Tsteps)
  odd_group = [2i - 1 for i in 1:div(mp.Tsteps, 2)]
  even_group = [2i for i in 1:div(mp.Tsteps - 1, 2)]
  r = (mp.Tmax - mp.Tmin) / (mp.Tsteps - 1)
  T = [mp.Tmin + r * (i - 1) for i in 1:mp.Tsteps]
  progress = Progress(mp.runs)
  for run in 1:mp.runs
    ireplica = collect(1:mp.Tsteps)
    phi = 2pi * rand(mp.Tsteps, mp.L, mp.L)
    costheta = 2 * rand(mp.Tsteps, mp.L, mp.L) .- 1
    #phi = zeros(L, L); costheta = zeros(L, L)
    energy = zeros(mp.Tsteps)
    energy2 = zeros(mp.Tsteps)
    open("T_HB_r$(run).dat", "w") do fp_T
      for mcs in 1:mp.mcs_max
        for i in ireplica
          print(fp_T, "$(T[i]) ")
        end
        print(fp_T, "\n")
        tmpE = zeros(mp.Tsteps)
        for Tstep in 1:mp.Tsteps
          phi[ireplica[Tstep], :, :], costheta[ireplica[Tstep], :, :] = update_HB!(mp, phi[ireplica[Tstep], :, :], costheta[ireplica[Tstep], :, :], T[Tstep])
          tmpE[Tstep] = calc_ene(mp, phi[ireplica[Tstep], :, :], costheta[ireplica[Tstep], :, :])
          if mcs > mp.discard
            energy[Tstep] += tmpE[Tstep]
            energy2[Tstep] += tmpE[Tstep]^2
          end
        end

        if mcs % 2 == 1 #odd group
          Tgroup = odd_group
        else #even group
          Tgroup = even_group
        end

        for Tstep in Tgroup
          irep1 = ireplica[Tstep]
          irep2 = ireplica[Tstep+1]
          if exp((1 / T[Tstep] - 1 / T[Tstep+1]) * (tmpE[Tstep+1] - tmpE[Tstep])) > rand()
            ireplica[Tstep], ireplica[Tstep+1] = irep2, irep1
          end
        end #for Tstep in Tgroup
      end #for mcs in 1:mp.mcs_max
    end
    #mcs average 
    energy ./= mp.frac * 2
    energy2 ./= mp.frac * 4
    spec = (energy2 - energy .^ 2) ./ T .^ 2 ./ mp.L^2
    ave_spec += spec
    var_spec += spec .^ 2
    ave_ene += energy
    var_ene += energy .^ 2
    next!(progress)
  end #for run in 1:mp.runs
  #run average
  ave_spec ./= mp.runs
  var_spec ./= mp.runs
  ave_ene ./= mp.runs
  var_ene ./= mp.runs
  open("E_HB.dat", "w") do fp_E
    open("C_HB.dat", "w") do fp_C
      for Tstep in 1:mp.Tsteps
        write(fp_E, "$(T[Tstep]) $(ave_ene[Tstep]) $(sqrt((var_ene[Tstep] - ave_ene[Tstep] ^ 2) / (mp.runs - 1)))\n")
        write(fp_C, "$(T[Tstep]) $(ave_spec[Tstep]) $(sqrt((var_spec[Tstep] - ave_spec[Tstep] ^ 2) / (mp.runs - 1)))\n")
      end
    end
  end
end

using ProgressMeter
L = 4
J1 = 1
J2 = 0
runs = 5
Tmin = 0.5
Tmax = 2
Tsteps = 100
#Tmin=2;Tmax=Tmin;Tsteps=1
mcs_max = 10000
discard = 8000
frac = mcs_max - discard
ip = zeros(Int64, L);
im = zeros(Int64, L);
for i in 1:L
  ip[i] = i + 1
  im[i] = i - 1
end;
ip[L] = 1;
im[1] = L;
modpara = mp(L, J1, J2, runs, Tmin, Tmax, Tsteps, mcs_max, discard, frac, ip, im)
@time main(modpara)
