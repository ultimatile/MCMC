#classical J1-J2 Heisenberg model
#Metropolis method
#Annealing
L=4
Tmin=1;Tmax=2;Tsteps=50
#Tmin=20;Tmax=Tmin;Tsteps=1
deltaT=(Tmax-Tmin)/Tsteps
J1=1;J2=0
sintheta(c)=sqrt(abs(1-c^2))
SiSj(pi,pj,ti,tj)=cos(pi-pj)*sintheta(ti)*sintheta(tj)+ti*tj
mcs_max=20000
runs=5
discard=18000
frac=mcs_max-discard
ave_spec=zeros(Tsteps);var_spec=zeros(Tsteps)
ave_ene=zeros(Tsteps);var_ene=zeros(Tsteps)
#phi=zeros(L,L);costheta=zeros(L,L)
ip=zeros(Int,L);im=zeros(Int,L)
for i in 1:L;ip[i]=i+1;im[i]=i-1;end
ip[L]=1;im[1]=L
for run in 1:runs
    phi=2pi*rand(L,L);costheta=2*rand(L,L).-1
    #phi=zeros(L,L);costheta=zeros(L,L)
    for Tstep in 1:Tsteps
        energy=0;energy2=0
        T=Tmax-deltaT*(Tstep-1)
        #auto=0
        for mcs in 1:mcs_max
            tmpH=0
            for ix in 1:L
                for iy in 1:L
                    tmpE=J1*(
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
                    newphi=2pi*rand();newcostheta=2*rand()-1
                    newE=J1*(
                    SiSj(newphi,phi[ix,ip[iy]],newcostheta,costheta[ix,ip[iy]])
                    +SiSj(newphi,phi[ix,im[iy]],newcostheta,costheta[ix,im[iy]])
                    +SiSj(newphi,phi[ip[ix],iy],newcostheta,costheta[ip[ix],iy])
                    +SiSj(newphi,phi[im[ix],iy],newcostheta,costheta[im[ix],iy])
                    )+J2*(
                    SiSj(newphi,phi[ip[ix],ip[iy]],newcostheta,costheta[ip[ix],ip[iy]])
                    +SiSj(newphi,phi[ip[ix],im[iy]],newcostheta,costheta[ip[ix],im[iy]])
                    +SiSj(newphi,phi[im[ix],ip[iy]],newcostheta,costheta[im[ix],ip[iy]])
                    +SiSj(newphi,phi[im[ix],im[iy]],newcostheta,costheta[im[ix],im[iy]])
                    )
                    deltaE=newE-tmpE
                    prob=rand()
                    if exp(-deltaE/T)>prob
                        phi[ix,iy]=newphi;costheta[ix,iy]=newcostheta
                        if mcs>discard;tmpH+=newE;end
                    else
                        if mcs>discard;tmpH+=tmpE;end
                    end
                end
            end
            energy+=tmpH;energy2+=tmpH^2
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
open( "E_MP.dat", "w" ) do fp_E
    open( "C_MP.dat", "w" ) do fp_C
        for Tstep in 1:Tsteps
            write( fp_E, "$(Tmax-deltaT*(Tstep-1)) $(ave_ene[Tstep]) $(sqrt((var_ene[Tstep]-ave_ene[Tstep]^2)/(runs-1)))\n" )
            write( fp_E, "$(Tmax-deltaT*(Tstep-1)) $(ave_spec[Tstep]) $(sqrt((var_spec[Tstep]-ave_spec[Tstep]^2)/(runs-1)))\n" )
        end
    end
end
