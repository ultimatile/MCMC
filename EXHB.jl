include("HB.jl")
#classical J1-J2 Heisenberg model
#Heatbath method
#replica exchange

# function observe_ene(mp, phi, costheta)
#   tmpE = 0
#   for ix in 1:mp.L
#     for iy in 1:mp.L
#       tmpE += mp.J1 * (
#         SiSj(phi[ix, iy], phi[ix, mp.inext[iy]], costheta[ix, iy], costheta[ix, mp.inext[iy]])
#         + SiSj(phi[ix, iy], phi[ix, mp.iprev[iy]], costheta[ix, iy], costheta[ix, mp.iprev[iy]])
#         + SiSj(phi[ix, iy], phi[mp.inext[ix], iy], costheta[ix, iy], costheta[mp.inext[ix], iy])
#         + SiSj(phi[ix, iy], phi[mp.iprev[ix], iy], costheta[ix, iy], costheta[mp.iprev[ix], iy])
#       ) + mp.J2 * (
#         SiSj(phi[ix, iy], phi[mp.inext[ix], mp.inext[iy]], costheta[ix, iy], costheta[mp.inext[ix], mp.inext[iy]])
#         + SiSj(phi[ix, iy], phi[mp.inext[ix], mp.iprev[iy]], costheta[ix, iy], costheta[mp.inext[ix], mp.iprev[iy]])
#         + SiSj(phi[ix, iy], phi[mp.iprev[ix], mp.inext[iy]], costheta[ix, iy], costheta[mp.iprev[ix], mp.inext[iy]])
#         + SiSj(phi[ix, iy], phi[mp.iprev[ix], mp.iprev[iy]], costheta[ix, iy], costheta[mp.iprev[ix], mp.iprev[iy]])
#       )
#     end
#   end
#   return tmpE
# end

function main(mp::ModelParameters)
  runs = mp.runs
  L = mp.L
  frac = mp.frac
  Tsteps = mp.Tsteps
  ave_ene = zeros(Tsteps)
  var_ene = zeros(Tsteps)
  ave_spec = zeros(Tsteps)
  var_spec = zeros(Tsteps)
  odd_group = [2i - 1 for i in 1:div(Tsteps, 2)]
  even_group = [2i for i in 1:div(Tsteps - 1, 2)]
  r = (mp.Tmax - mp.Tmin) / (Tsteps - 1)
  T = [mp.Tmin + r * (i - 1) for i in 1:Tsteps]
  @showprogress for run in 1:runs
    ireplica = collect(1:Tsteps)
    phi = 2pi * rand(Tsteps, L, L)
    costheta = 2 * rand(Tsteps, L, L) .- 1
    #phi = zeros(L, L); costheta = zeros(L, L)
    energy = zeros(Tsteps)
    energy2 = zeros(Tsteps)
    open("T_HB_r$(run).dat", "w") do fp_T
      for mcs in 1:mp.mcs_max
        for i in ireplica
          print(fp_T, "$(T[i]) ")
        end
        print(fp_T, "\n")
        tmpE = zeros(Tsteps)
        for Tstep in 1:Tsteps
          # phi[ireplica[Tstep], :, :], costheta[ireplica[Tstep], :, :] = 
          update_HB!(mp, phi[ireplica[Tstep], :, :], costheta[ireplica[Tstep], :, :], T[Tstep])
          tmpE[Tstep] = observe_ene(mp, phi[ireplica[Tstep], :, :], costheta[ireplica[Tstep], :, :])
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
    energy ./= frac * 2
    energy2 ./= frac * 4
    spec = (energy2 - energy .^ 2) ./ T .^ 2 ./ L^2
    ave_spec += spec
    var_spec += spec .^ 2
    ave_ene += energy
    var_ene += energy .^ 2
  end #for run in 1:mp.runs
  #run average
  ave_spec ./= runs
  var_spec ./= runs
  ave_ene ./= runs
  var_ene ./= runs
  open("E_EXHB.dat", "w") do fp_E
    open("C_EXHB.dat", "w") do fp_C
      for Tstep in 1:Tsteps
        write(fp_E, "$(T[Tstep]) $(ave_ene[Tstep]) $(sqrt((var_ene[Tstep] - ave_ene[Tstep] ^ 2) / (mp.runs - 1)))\n")
        write(fp_C, "$(T[Tstep]) $(ave_spec[Tstep]) $(sqrt((var_spec[Tstep] - ave_spec[Tstep] ^ 2) / (mp.runs - 1)))\n")
      end
    end
  end
end

# L = 4
# J1 = 1
# J2 = 0
# runs = 5
# Tmin = 0.5
# Tmax = 2
# Tsteps = 100
# #Tmin=2;Tmax=Tmin;Tsteps=1
# mcs_max = 10000
# discard = 8000
# frac = mcs_max - discard
# inext = zeros(Int64, L);
# # iprev = zeros(Int64, L);
# # for i in 1:L
# #   inext[i] = i + 1
# #   iprev[i] = i - 1
# # end;
# # inext[L] = 1;
# # iprev[1] = L;
# modpara = mp(L, J1, J2, runs, Tmin, Tmax, Tsteps, mcs_max, discard, frac, inext, iprev)
# main(modpara)

function main(ARGS)
  L = 4
  J1 = 1.0
  J2 = 0.0
  runs = 5
  Tmin = 0.5
  Tmax = 2.0
  Tsteps = 50
  #Tmin=2;Tmax=Tmin;Tsteps=1
  mcs_max = 20000
  discard = 10000
  return main(mp(L, J1, J2, runs, Tmin, Tmax, Tsteps, mcs_max, discard))
end

@main
