#classical J1-J2 Heisenberg model
#Heatbath method
#Annealing
L=4
Tmin=1;Tmax=2;Tsteps=50
#Tmin=2;Tmax=Tmin;Tsteps=1
deltaT=(Tmax-Tmin)/Tsteps
J1=1;J2=0
sintheta(c)=sqrt(abs(1-c^2))
SiSj(pi,pj,ti,tj)=cos(pi-pj)*sintheta(ti)*sintheta(tj)+ti*tj
mcs_max=10000
runs=5
discard=8000
frac=mcs_max-discard
ave_ene=zeros(Tsteps);var_ene=zeros(Tsteps)
ave_spec=zeros(Tsteps);var_spec=zeros(Tsteps)
ip=zeros(Int,L);im=zeros(Int,L)
for i in 1:L;ip[i]=i+1;im[i]=i-1;end
ip[L]=1;im[1]=L
for run in 1:runs
    phi=2pi*rand!(iphi);costheta=rand!(icostheta)
    #phi=zeros(L,L);costheta=zeros(L,L)
    for Tstep in 1:Tsteps
        energy=0;energy2=0
        T=Tmax-deltaT*(Tstep-1)
        for mcs in 1:mcs_max
            for ix in 1:L
                for iy in 1:L
                    sinthetapx=sintheta(costheta[ip[ix],iy])
                    sinthetamx=sintheta(costheta[im[ix],iy])
                    sinthetapy=sintheta(costheta[ix,ip[iy]])
                    sinthetamy=sintheta(costheta[ix,im[iy]])
                    sinthetapxpy=sintheta(costheta[ip[ix],ip[iy]])
                    sinthetapxmy=sintheta(costheta[ip[ix],im[iy]])
                    sinthetamxpy=sintheta(costheta[im[ix],ip[iy]])
                    sinthetamxmy=sintheta(costheta[im[ix],im[iy]])
                    hlocalx=(J1*(cos(phi[ip[ix],iy])*sinthetapx+cos(phi[im[ix],iy])*sinthetamx
                    +cos(phi[ix,ip[iy]])*sinthetapy+cos(phi[ix,im[iy]])*sinthetamy)
                    +J2*(cos(phi[ip[ix],ip[iy]])*sinthetapxpy+cos(phi[ip[ix],im[iy]])*sinthetapxmy
                    +cos(phi[im[ix],ip[iy]])*sinthetamxpy+cos(phi[im[ix],im[iy]])*sinthetamxmy))
                    hlocaly=(J1*(sin(phi[ip[ix],iy])*sinthetapx+sin(phi[im[ix],iy])*sinthetamx
                    +sin(phi[ix,ip[iy]])*sinthetapy+sin(phi[ix,im[iy]])*sinthetamy)
                    +J2*(sin(phi[ip[ix],ip[iy]])*sinthetapxpy+sin(phi[ip[ix],im[iy]])*sinthetapxmy
                    +sin(phi[im[ix],ip[iy]])*sinthetamxpy+sin(phi[im[ix],im[iy]])*sinthetamxmy))
                    hlocalz=(J1*(costheta[ip[ix],iy]+costheta[im[ix],iy]+costheta[ix,ip[iy]]+costheta[ix,im[iy]])
                    +J2*(costheta[ip[ix],ip[iy]]+costheta[ip[ix],im[iy]]+costheta[im[ix],ip[iy]]+costheta[im[ix],im[iy]]))
                    betaH=sqrt(hlocalx^2+hlocaly^2+hlocalz^2)/T
                    rndm=rand()
                    sznew=-1-log1p(rndm*expm1(-2*betaH))/betaH
                    sinthetanew=sintheta(sznew)
                    cpsi=hlocalz/sqrt(hlocalx^2+hlocaly^2+hlocalz^2)
                    spsi=sintheta(cpsi)
                    cphi=hlocalx/sqrt(hlocalx^2+hlocaly^2)
                    sphi=hlocaly/sqrt(hlocalx^2+hlocaly^2)
                    rndm=rand()
                    sxnew=cos(2pi*rndm)*sinthetanew
                    synew=sin(2pi*rndm)*sinthetanew
                    sx=cpsi*cphi*sxnew-sphi*synew+spsi*cphi*sznew
                    sy=cpsi*sphi*sxnew+cphi*synew+spsi*sphi*sznew
                    costheta[ix,iy]=-spsi*sxnew+cpsi*sznew
                    sinthetanew=sintheta(costheta[ix,iy])
                    if sy/sinthetanew > 1
                        phi[ix,iy]=0.5pi
                    elseif sx/sinthetanew > 1
                        phi[ix,iy]=0
                    elseif sx/sinthetanew < -1
                        phi[ix,iy]=pi
                    elseif sy/sinthetanew < -1
                        phi[ix,iy]=1.5pi
                    elseif sy/sinthetanew > 0
                        phi[ix,iy]=acos(sx/sinthetanew)
                    else
                        phi[ix,iy]=2pi-acos(sx/sinthetanew)
                    end
                end
            end
            if mcs>discard
                tmpE=0
                for ix in 1:L
                    for iy in 1:L
                        tmpE+=J1*(
                        SiSj(phi[ix,iy],phi[ix,ip[iy]],costheta[ix,iy],costheta[ix,ip[iy]])
                        +SiSj(phi[ix,iy],phi[ix,im[iy]],costheta[ix,iy],costheta[ix,im[iy]])
                        +SiSj(phi[ix,iy],phi[ip[ix],iy],costheta[ix,iy],costheta[ip[ix],iy])
                        +SiSj(phi[ix,iy],phi[im[ix],iy],costheta[ix,iy],costheta[im[ix],iy])
                        )+J2*(
                        SiSj(phi[ix,iy],phi[ip[ix],ip[iy]],costheta[ix,iy],costheta[ip[ix],ip[iy]])
                        +SiSj(phi[ix,iy],phi[ip[ix],im[iy]],costheta[ix,iy],costheta[ip[ix],im[iy]])
                        +SiSj(phi[ix,iy],phi[im[ix],ip[iy]],costheta[ix,iy],costheta[im[ix],ip[iy]])
                        +SiSj(phi[ix,iy],phi[im[ix],im[iy]],costheta[ix,iy],costheta[im[ix],im[iy]])
                        )
                    end
                end
                energy+=tmpE;energy2+=tmpE^2
            end
        end
        energy/=frac*2;energy2/=frac*4
        spec=(energy2-energy^2)/T^2/L^2
        ave_spec[Tstep]+=spec
        var_spec[Tstep]+=spec^2
        ave_ene[Tstep]+=energy
        var_ene[Tstep]+=energy^2
    end
end
ave_spec/=runs;var_spec/=runs
ave_ene/=runs;var_ene/=runs
open( "E_HB.dat", "w" ) do fp_E
    open( "C_HB.dat", "w" ) do fp_C
        for Tstep in 1:Tsteps
            write( fp_E, "$(Tmax-deltaT*(Tstep-1)) $(ave_ene[Tstep]) $(sqrt((var_ene[Tstep]-ave_ene[Tstep]^2)/(runs-1)))\n" )
            write( fp_E, "$(Tmax-deltaT*(Tstep-1)) $(ave_spec[Tstep]) $(sqrt((var_spec[Tstep]-ave_spec[Tstep]^2)/(runs-1)))\n" )
        end
    end
end
