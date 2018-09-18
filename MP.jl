#classical J1-J2 Heisenberg model
#Metropolis method
#Annealing
struct mp
    L::Int64
    J1::Float64 #exchange coupling
    J2::Float64 #hopping
    runs::Int64
    Tmin :: Float64
    Tmax :: Float64
    Tsteps :: Int64
    deltaT :: Float64
    mcs_max::Int64 #Maximum number of MC steps
    discard::Int64 #Thermalization
    frac :: Int64
end
sintheta(c)=sqrt(abs(1-c^2))
SiSj(pi,pj,ti,tj)=cos(pi-pj)*sintheta(ti)*sintheta(tj)+ti*tj
function main(mp)
    ave_ene=zeros(mp.Tsteps);var_ene=zeros(mp.Tsteps)
    ave_spec=zeros(mp.Tsteps);var_spec=zeros(mp.Tsteps)
    ip=zeros(Int,mp.L);im=zeros(Int,mp.L);for i in 1:mp.L;ip[i]=i+1;im[i]=i-1;end;ip[mp.L]=1;im[1]=mp.L
    for run in 1:mp.runs
        phi=2pi*rand(mp.L,mp.L);costheta=2*rand(mp.L,mp.L).-1
        #phi=zeros(L,L);costheta=zeros(L,L)
        for Tstep in 1:mp.Tsteps
            energy=0;energy2=0
            T=mp.Tmax-mp.deltaT*(Tstep-1)
            #auto=0
            for mcs in 1:mp.mcs_max
                tmpH=0
                for ix in 1:mp.L
                    for iy in 1:mp.L
                        tmpE=mp.J1*(
                        SiSj(phi[ix,iy],phi[ix,ip[iy]],costheta[ix,iy],costheta[ix,ip[iy]])
                        +SiSj(phi[ix,iy],phi[ix,im[iy]],costheta[ix,iy],costheta[ix,im[iy]])
                        +SiSj(phi[ix,iy],phi[ip[ix],iy],costheta[ix,iy],costheta[ip[ix],iy])
                        +SiSj(phi[ix,iy],phi[im[ix],iy],costheta[ix,iy],costheta[im[ix],iy])
                        )+mp.J2*(
                        SiSj(phi[ix,iy],phi[ip[ix],ip[iy]],costheta[ix,iy],costheta[ip[ix],ip[iy]])
                        +SiSj(phi[ix,iy],phi[ip[ix],im[iy]],costheta[ix,iy],costheta[ip[ix],im[iy]])
                        +SiSj(phi[ix,iy],phi[im[ix],ip[iy]],costheta[ix,iy],costheta[im[ix],ip[iy]])
                        +SiSj(phi[ix,iy],phi[im[ix],im[iy]],costheta[ix,iy],costheta[im[ix],im[iy]])
                        )
                        newphi=2pi*rand();newcostheta=2*rand()-1
                        newE=mp.J1*(
                        SiSj(newphi,phi[ix,ip[iy]],newcostheta,costheta[ix,ip[iy]])
                        +SiSj(newphi,phi[ix,im[iy]],newcostheta,costheta[ix,im[iy]])
                        +SiSj(newphi,phi[ip[ix],iy],newcostheta,costheta[ip[ix],iy])
                        +SiSj(newphi,phi[im[ix],iy],newcostheta,costheta[im[ix],iy])
                        )+mp.J2*(
                        SiSj(newphi,phi[ip[ix],ip[iy]],newcostheta,costheta[ip[ix],ip[iy]])
                        +SiSj(newphi,phi[ip[ix],im[iy]],newcostheta,costheta[ip[ix],im[iy]])
                        +SiSj(newphi,phi[im[ix],ip[iy]],newcostheta,costheta[im[ix],ip[iy]])
                        +SiSj(newphi,phi[im[ix],im[iy]],newcostheta,costheta[im[ix],im[iy]])
                        )
                        deltaE=newE-tmpE
                        prob=rand()
                        if exp(-deltaE/T)>prob
                            phi[ix,iy]=newphi;costheta[ix,iy]=newcostheta
                            if mcs>mp.discard;tmpH+=newE;end
                        else
                            if mcs>mp.discard;tmpH+=tmpE;end
                        end
                    end
                end
                energy+=tmpH;energy2+=tmpH^2
            end
            energy/=mp.frac*2;energy2/=mp.frac*4
            spec=(energy2-energy^2)/T^2/mp.L^2
            ave_spec[Tstep]+=spec
            var_spec[Tstep]+=spec^2
            ave_ene[Tstep]+=energy
            var_ene[Tstep]+=energy^2
           end
    end
    ave_spec/=mp.runs;var_spec/=mp.runs
    ave_ene/=mp.runs;var_ene/=mp.runs
    open( "E_MP.dat", "w" ) do fp_E
        open( "C_MP.dat", "w" ) do fp_C
            for Tstep in 1:mp.Tsteps
                write( fp_E, "$(mp.Tmax-mp.deltaT*(Tstep-1)) $(ave_ene[Tstep]) $(sqrt((var_ene[Tstep]-ave_ene[Tstep]^2)/(mp.runs-1)))\n" )
                write( fp_C, "$(mp.Tmax-mp.deltaT*(Tstep-1)) $(ave_spec[Tstep]) $(sqrt((var_spec[Tstep]-ave_spec[Tstep]^2)/(mp.runs-1)))\n" )
            end
        end
    end
end

const L=4
const J1=1
const J2=0
const runs=5
const Tmin=1
const Tmax=2
const Tsteps=50
#Tmin=2;Tmax=Tmin;Tsteps=1
const deltaT=(Tmax-Tmin)/(Tsteps-1)
const mcs_max=10000
const discard=8000
const frac=mcs_max-discard
modpara=mp(L,J1,J2,runs,Tmin,Tmax,Tsteps,deltaT,mcs_max,discard,frac)
@time main(modpara)
