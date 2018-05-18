       integer iflag,ipdc,nmol,iloop,iand
      parameter(nmol=2)
          real*8 acc1(nmol,3,3),acc(nmol,3,3),fmass(3)
           real*8 sumdx,sumdy,sumdz
         real*8  potcc(nmol,3),potcd(nmol,3),force(nmol,3,3)
          real*8 dipole(nmol,3,3),dipole_p(nmol,3,3)
       real*8 rat(nmol,3,3),rf(nmol,3,3)
        real*8 cvec(3,3),civec(3,3)
         real*8 qat(3,nmol),qatd(3,nmol),alph(nmol,3) 
       real*8 efdd(nmol,3,3),efdc(nmol,3,3),fac(4)
           real*8 fq(nmol,3,3,3),pot(nmol,3)
       double precision dipdip(0:3),chgchg(0:3),chgdip(0:3)
         real*8 kcal_j,debfac,pi,rcut,ewfac,sqpii!,fac(4)
         real*8 force1(nmol,3,3),force2(nmol,3,3)
          real*8 dforce(nmol,3,3),dforce1(nmol,3,3)
          real*8 dforce2(nmol,3,3),diporig(nmol,3,3)
           real*8 diptens1(nmol,nmol,3,3),force_c(nmol,3,3)
           real*8 diptens2(nmol,nmol,3,3)
           real*8 diptens3(nmol,nmol,3,3)
           real*8 diptens4(nmol,nmol,3,3)
           real*8 diptens5(nmol,nmol,3,3)
           real*8 diptens6(nmol,nmol,3,3),dipole_t(nmol,3,3)
             real*8 efext(3),echarge
             parameter( echarge=1.60217733d-19)
             parameter(ewfac=.4d0)
             parameter(SQPII=0.56418958354D0)
             real*8 c_16,c_14,c_12,c_6
             real*8 c16,c14,c12,c6,epso,esconvert
                   real*8 dpdr(128,3,3,3,3)
             character*30, xyzFILE,velFILE,cvecFILE,xyZOUT,velOUT,tempOUT,dipOUT,peOUT,keOUT,totOUT,md_trajOUT,md_velOUT,grooOUT,gromOUT
             common/files/xyzFILE,velFILE,cvecFILE,xyZOUT,velOUT,tempOUT,dipOUT,peOUT,keOUT,totOUT,md_trajOUT,md_velOUT,grooOUT,gromOUT
             common/for_dp/dpdr
         common/dipwork/efdc,efdd,efext,pot,potcc,potcd,force
          common/dip/dipole,dipole_p,dforce,dipole_t
        common/part/fq
        common/pos/rat,rf
         common/dipp2/sumdx,sumdy,sumdz
        common/dipole_tensor1/diptens1,diptens2,diptens3
        common/dipole_tensor2/diptens4,diptens5,diptens6
        common/box/cvec,civec
          common/charge/qat,qatd
       common/loop/iloop,iflag,iand,ipdc
        common/polar/alph
         common/constant/fac,pi,c16,c14,c12,c6
         common/constant/c_16,c_14,c_12,c_6
         common/force_cou/force_c        
 
