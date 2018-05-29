

      subroutine real(utot,upair,ucc,ucd,udd,uspring,uintra,ucx,udx)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       real*8 axes(3,3),raxes(3,3),waxes(3,3)
        real*8 gradq(3,3,3)
       include 'parameter.i'
c       include 'parameter2.i'

       real*8 r1(3,3),bond(nmol,3,6),dr1(3,3)                  
        debfac=echarge*1d-10/3.33d-30

         efext(1)=0.d0
         efext(2)=0.d0
         efext(3)=0.d0
      uintra=0.d0
      utot=0.d0
      ucc=0.d0
      ucd=0.d0
      udd=0.d0
      uspring=0.d0
      upair=0.d0
      ucx=0.d0
      udx=0.d0
      
      rcutc=cvec(1,1)/2.d0
      rcut2=rcutc*rcutc*100.d0
       iltt=0                     
       iltt1=0                     
       iltt2=0                     
       iltt3=0                     

c     sum over all intra interactions which don't involve m-sites
         if(iand.eq.1.and.iflag.eq.1)then           
      do i=2,nmol

         do ia=1,3
                  do  idim=1,3
            r1(idim,ia)=rat(i,ia,idim)-rat(i,1,idim)
              r1(idim,ia)=r1(idim,ia)-anint(r1(idim,ia)/cvec(1,1))
     : *cvec(1,1)
                 enddo
                 enddo
                 call pot_nasa(r1,dr1,eng1)
c               eng1=eng1!*kcal_j
               do ia=1,3
             do idim=1,3
          force(i,ia,idim)=-dr1(idim,ia)+
     :  force(i,ia,idim)
             enddo
             enddo

                uintra=uintra+eng1
               enddo
           endif
c     sum over molecules to determine intramolecular potential energy
c     and determine charges from Partridge dipole surface

      do i=2,nmol
         r1x=rat(i,2,1)-rat(i,1,1)
         r1x=r1x-anint(r1x/cvec(1,1))*cvec(1,1) 
         r1y=rat(i,2,2)-rat(i,1,2)
         r1y=r1y-anint(r1y/cvec(2,2))*cvec(2,2) 
         r1z=rat(i,2,3)-rat(i,1,3)
         r1z=r1z-anint(r1z/cvec(3,3))*cvec(3,3) 

         r2x=rat(i,3,1)-rat(i,1,1)
         r2x=r2x-anint(r2x/cvec(1,1))*cvec(1,1) 
         r2y=rat(i,3,2)-rat(i,1,2)
         r2y=r2y-anint(r2y/cvec(2,2))*cvec(2,2) 
         r2z=rat(i,3,3)-rat(i,1,3)
         r2z=r2z-anint(r2z/cvec(3,3))*cvec(3,3) 

c        r1x=r1x*1d+10
c        r1y=r1y*1d+10
c        r1z=r1z*1d+10
c        r2x=r2x*1d+10
c        r2y=r2y*1d+10
c        r2z=r2z*1d+10
             bond(i,1,1)=r1x
             bond(i,1,2)=r1y
             bond(i,1,3)=r1z
             bond(i,2,1)=r2x
             bond(i,2,2)=r2y
             bond(i,2,3)=r2z
              bond(i,1,4)=r1x**2+r1y**2+r1z**2
              bond(i,2,4)=r2x**2+r2y**2+r2z**2
              r1mag=dsqrt(bond(i,1,4))
              r2mag=dsqrt(bond(i,2,4))
            call dip2_h2o(r1x,r1y,r1z,r2x,r2y,r2z,r1mag,r2mag
     $           ,qh1,qh2,qo
     $           ,gradq)
            if(r1mag.gt.1.2d0.or.r1mag.lt.0.8d0)then
            print*,'prob',i,r1mag,r2mag
             xx1=rat(i,1,1)-rat(1,2,1)
                xx1=xx1-anint(xx1/cvec(1,1))*cvec(1,1)
             yy1=rat(i,1,2)-rat(1,2,2)
                yy1=yy1-anint(yy1/cvec(1,1))*cvec(1,1)
             zz1=rat(i,1,3)-rat(1,2,3)
 
              zz1=zz1-anint(zz1/cvec(1,1))*cvec(1,1)
              print*,dsqrt(xx1**2+yy1**2+zz1**2)
            endif
            if(r2mag.gt.1.2d0.or.r2mag.lt.0.8d0)then
            print*,'prob1',i,r1mag,r2mag
             xx1=rat(i,1,1)-rat(1,3,1)
                xx1=xx1-anint(xx1/cvec(1,1))*cvec(1,1)
             yy1=rat(i,1,2)-rat(1,3,2)
                yy1=yy1-anint(yy1/cvec(1,1))*cvec(1,1)
             zz1=rat(i,1,3)-rat(1,3,3)
 
              zz1=zz1-anint(zz1/cvec(1,1))*cvec(1,1)
              print*,dsqrt(xx1**2+yy1**2+zz1**2)
            endif

          qat(1,i)=qo         
          qat(2,i)=qh1         
          qat(3,i)=qh2        
c     store charge gradients in array
         
      do  ivec=1,3
         do  inuc=1,3  
         fq(i,inuc,1,ivec)=
     $        gradq(inuc,1,ivec)!*echarge/1d-10
         fq(i,inuc,2,ivec)=
     $        gradq(inuc,2,ivec)!*echarge/1d-10
         
         fq(i,inuc,3,ivec)=
     $        -(fq(i,inuc,1,ivec)+fq(i,inuc,2,ivec))
          enddo     
          enddo     
         
           enddo
                  
         
         permhz=0.00000000d0
         perm_oz=0.2962054d0
c     H permanent dipole in e A
         permhx=0.10672099d0
         permhy=0.13501661d0
     
c     P_para = Px*cos(theta) + Py*sin(theta)
c     P_perp = -Px*sin(theta) + Py*cos(theta)
            
         permr=(permhx*0.61207927d0 + permhy*0.79079641d0)
         permw=(-permhx*0.79079641d0 + permhy*0.61207927d0)
         
c     P_perp = Px (w1 . x) + Py (w1 . y)
c     sum over molecules
          dipole_p(1,1,1)=0.d0
          dipole_p(1,1,2)=0.d0
          dipole_p(1,1,3)=0.d0
          dipole_p(1,2,1)=0.d0
          dipole_p(1,2,2)=0.d0
          dipole_p(1,2,3)=0.d0
          dipole_p(1,3,1)=0.d0
          dipole_p(1,3,2)=0.d0
          dipole_p(1,3,3)=0.d0
         do i=2,nmol

             imol=i
            r1dotr2 = 0.d0
 
c     1st H coordinates relative to O site
         r1x=bond(i,1,1)
         r1y=bond(i,1,2)
         r1z=bond(i,1,3)
         r1mag=dsqrt(bond(i,1,4))
         raxes(1,1)=r1x
         raxes(1,2)=r1y
         raxes(1,3)=r1z
c     2nd H coordinates relative to O site
         r2x=bond(i,2,1)
         r2y=bond(i,2,2)  
         r2z=bond(i,2,3)
         r2mag=dsqrt(bond(i,2,4))
         
         raxes(2,1)=r2x
         raxes(2,2)=r2y
         raxes(2,3)=r2z
                    
c     define bisector axis (body-centered x-axis xb)
         xx=r1x+r2x
         yy=r1y+r2y
         zz=r1z+r2z
         dis=dsqrt(xx*xx+yy*yy+zz*zz)

         xx=xx/dis
         yy=yy/dis
         zz=zz/dis
     
         r1r2dis=dis
         r1r2dis3=dis**3
     
         axes(1,1)=xx
         axes(1,2)=yy
         axes(1,3)=zz
      
c     define out of plane axis from cross product (body-centered y-axis yb = r1^r2)
         xx=r1y*r2z-r2y*r1z
         yy=r1z*r2x-r2z*r1x
         zz=r1x*r2y-r2x*r1y
         dis=sqrt(xx*xx+yy*yy+zz*zz)
         r1cr2mag=dis

         xx=xx/dis
         yy=yy/dis
         zz=zz/dis
         
c     axes(2,i) is the cross product vector unit vector
         
         axes(2,1)=xx
         axes(2,2)=yy
         axes(2,3)=zz

c     define w1, w2 axes
         waxes(1,1) = (r1y*zz - r1z*yy)/r1mag
         waxes(1,2) = (r1z*xx - r1x*zz)/r1mag
         waxes(1,3) = (r1x*yy - r1y*xx)/r1mag
         
         waxes(2,1) = (r2y*zz - r2z*yy)/r2mag
         waxes(2,2) = (r2z*xx - r2x*zz)/r2mag
         waxes(2,3) = (r2x*yy - r2y*xx)/r2mag
         
         w1dotx=waxes(1,1)*axes(1,1)+waxes(1,2)*axes(1,2)
     $        +waxes(1,3)*axes(1,3)
         sign1=-w1dotx/abs(w1dotx)
         w2dotx=waxes(2,1)*axes(1,1)+waxes(2,2)*axes(1,2)
     $        +waxes(2,3)*axes(1,3)
         sign2=-w2dotx/(abs(w2dotx))
         
         r1dotr2 = 0.d0
         r1dotr1 = 0.d0
         
         do kvec=1,3   
            r1dotr2 = r1dotr2 + raxes(1,kvec)*raxes(2,kvec)
            r1dotr1 = r1dotr1 + raxes(1,kvec)*raxes(1,kvec)
         enddo
         
         
c     O permanent dipole (pointing in direction of first axis)
c     O permanent dipole (pointing in direction of xb axis)

         dipole_p(i,1,1)=axes(1,1)*perm_oz
         dipole_p(i,1,2)=axes(1,2)*perm_oz
         dipole_p(i,1,3)=axes(1,3)*perm_oz
         
         dipole_p(i,2,1)=waxes(1,1)*sign1*permw
     $        + permr*raxes(1,1)/r1mag
         dipole_p(i,2,2)=waxes(1,2)*sign1*permw
     $        + permr*raxes(1,2)/r1mag
         dipole_p(i,2,3)=waxes(1,3)*sign1*permw
     $        + permr*raxes(1,3)/r1mag
         
         dipole_p(i,3,1)=waxes(2,1)*sign2*permw
     $        + permr*raxes(2,1)/r2mag
         dipole_p(i,3,2)=waxes(2,2)*sign2*permw
     $        + permr*raxes(2,2)/r2mag
         dipole_p(i,3,3)=waxes(2,3)*sign2*permw
     $        + permr*raxes(2,3)/r2mag
c     calculate dpdr terms 
         if(iflag.eq.1) then
         
c     zero dpdr terms
         do iatom=1,3
            do jatom=1,3
               do ivec=1,3
                  do jvec=1,3
                   dpdr(imol,iatom,jatom,ivec,jvec) = 0.d0
                  enddo
               enddo
            enddo
         enddo
         
         do ivec=1,3
            do jvec=1,3 
               if(ivec.eq.jvec) then
c     gradients for O
                  dpdr(imol,3,1,ivec,jvec) =  perm_oz/r1r2dis
         
c     gradients for w1
                  dpdr(imol,1,1,ivec,jvec) = dpdr(imol,1,1,ivec,jvec)
     $                 + sign1*permw*r1dotr2/(r1mag*r1cr2mag)
     $                 + permr/r1mag
         
                  dpdr(imol,1,2,ivec,jvec) = dpdr(imol,1,2,ivec,jvec)
     $                 - sign1*permw*r1mag**2/(r1mag*r1cr2mag)
         
c     gradients for w2
         
                  dpdr(imol,2,1,ivec,jvec) = dpdr(imol,2,1,ivec,jvec)
     $                 + sign2*permw*r2mag**2/(r2mag*r1cr2mag)
         
                  dpdr(imol,2,2,ivec,jvec) = dpdr(imol,2,2,ivec,jvec)
     $                 - sign2*permw*r1dotr2/(r2mag*r1cr2mag)
     $                 + permr/r2mag
            
               endif
         
c     gradients for O
         
              fac1 = (bond(imol,1,ivec)+bond(imol,2,ivec))
              fac2 = (bond(imol,1,jvec)+bond(imol,2,jvec))
              dpdr(imol,3,1,ivec,jvec) = dpdr(imol,3,1,ivec,jvec)
     $             - perm_oz*(fac1*fac2)/r1r2dis3
              dpdr(imol,3,2,ivec,jvec) =  dpdr(imol,3,1,ivec,jvec)

c     gradients for w1
              firstterm = (raxes(1,ivec)*raxes(2,jvec)
     $             - 2.d0*raxes(2,ivec)*raxes(1,jvec)) /
     $             (r1mag*r1cr2mag)
         
              secterm = (waxes(1,ivec)/(r1mag*r1cr2mag))*
     $             ((raxes(1,jvec)*r1cr2mag)/r1mag +
     $             (r1mag*r2mag*waxes(2,jvec)))
     
              thirdterm = permr*raxes(1,ivec)*raxes(1,jvec)/
     $             r1mag**3
         
              dpdr(imol,1,1,ivec,jvec) = dpdr(imol,1,1,ivec,jvec)
     $             + sign1*permw*(firstterm - secterm)
     $             - thirdterm
     
              firstterm = (raxes(1,ivec)*raxes(1,jvec)) /
     $             (r1mag*r1cr2mag)   
         
              secterm = (waxes(1,ivec)*waxes(1,jvec)*r1mag**2) /
     $              (r1mag*r1cr2mag)
         
              dpdr(imol,1,2,ivec,jvec) = dpdr(imol,1,2,ivec,jvec)
     $             + sign1*permw*(firstterm + secterm)
            
c     gradients for w2
              firstterm = (raxes(2,ivec)*raxes(2,jvec)) /
     $             (r2mag*r1cr2mag)
                  
              secterm = (waxes(2,ivec)*waxes(2,jvec)*r2mag**2) /
     $             (r2mag*r1cr2mag)
         
              dpdr(imol,2,1,ivec,jvec) = dpdr(imol,2,1,ivec,jvec)
     $             + sign2*permw*(- firstterm - secterm)
            
               
              firstterm = (-raxes(2,ivec)*raxes(1,jvec)
     $             + 2.d0*raxes(1,ivec)*raxes(2,jvec)) /
     $             (r2mag*r1cr2mag)

              secterm = (waxes(2,ivec)/(r2mag*r1cr2mag))*
     $             ((raxes(2,jvec)*r1cr2mag)/r2mag -
     $             (r2mag*r1mag*waxes(1,jvec)))
         
              thirdterm = permr*raxes(2,ivec)*raxes(2,jvec)/
     $             r2mag**3
         
              dpdr(imol,2,2,ivec,jvec) = dpdr(imol,2,2,ivec,jvec)
     $             + sign2*permw*(firstterm - secterm)
     $             - thirdterm
     

           enddo

        enddo
     
c     obtain dpdr on r=3 (oxygen) sites from negative of dpdr on other
c     two sites.
                              

        do ivec=1,3
           do jvec=1,3
              do jj=1,3
                 dpdr(imol,jj,3,ivec,jvec) = dpdr(imol,jj,3,ivec,jvec)
     $                - dpdr(imol,jj,1,ivec,jvec)
     $                - dpdr(imol,jj,2,ivec,jvec)
              enddo
           enddo
        enddo
        !print*,dpdr(1,3,3,1,1)
       !stop
      endif

       enddo
         
c     sum over all intra interactions

      do i=2,nmol
           imol=i
         do iia=1,2
             j=i
            jmol=j
         do jja=iia+1,3
                  qpe=0.0
                  gcfac=0.0
                  gfac=0.0
                  dfacx=0.0
                  dfacy=0.0
                  dfacz=0.0
                  ecc=0
                  ecd=0
                  edd=0

                  xdif=rat(i,iia,1)-rat(j,jja,1)
                 xdif=xdif-anint(xdif/cvec(1,1))*cvec(1,1)  
                  ydif=rat(i,iia,2)-rat(j,jja,2)
                 ydif=ydif-anint(ydif/cvec(2,2))*cvec(2,2)  
                  zdif=rat(i,iia,3)-rat(j,jja,3)
                 zdif=zdif-anint(zdif/cvec(3,3))*cvec(3,3)  

              
              rsq=xdif*xdif+ydif*ydif+zdif*zdif
              if(rsq.le.rcut2) then    !2
                 rsqi=1.d0/rsq
                 rr=dsqrt(rsq)
             iltt=iltt+1 
                !print*,rr
                 ia=0
              if(rr.le.rcutc)then          !3 
           call smear(dipdip,chgchg,chgdip,rr,rsq,i,j,
     $            iia,jja,ia)
            call field(dipdip,chgchg,chgdip,i,j,iia,jja,ia,
     $           xdif,ydif,zdif,rr,gcfac,dfacx,dfacy,dfacz,
     $           ecc,ecd,edd)
           !print*,rr!,rcutc!ecd,i,j,iia,jja,rr
        endif   !33
     
       if(iflag.eq.1) then  !5
         ucc=ucc+ecc
         ucd=ucd+ecd
         udd=udd+edd
           !print*,ucc,ucd,udd
         !stop
                  
         FORCE(i,iia,1)=FORCE(I,iia,1)+(GFAC+GCFAC)*XDIF+dfacx
         FORCE(i,iia,2)=FORCE(I,iia,2)+(GFAC+GCFAC)*YDIF+dfacy
         FORCE(i,iia,3)=FORCE(I,iia,3)+(GFAC+GCFAC)*ZDIF+dfacz
         FORCE(j,jja,1)=FORCE(J,jja,1)-(GFAC+GCFAC)*XDIF-dfacx
         FORCE(j,jja,2)=FORCE(J,jja,2)-(GFAC+GCFAC)*YDIF-dfacy
         FORCE(j,jja,3)=FORCE(J,jja,3)-(GFAC+GCFAC)*ZDIF-dfacz
         
      endif  !55

       endif   !22
                     
              
         enddo         
         enddo         
         enddo         
         !stop
      if(iand.eq.1)then
      do i=1,nmol-1
           imol=i
           if(i.eq.1)i_a=1
           if(i.gt.1)i_a=3
         do iia=1,i_a
       do j=i+1,nmol
           if(j.eq.9999)j_a=1
           if(j.gt.0)j_a=3

            jmol=j
         do jja=1,j_a
                  qpe=0.0
                  gcfac=0.0
                  gfac=0.0
                  dfacx=0.0
                  dfacy=0.0
                  dfacz=0.0
                  ecc=0
                  ecd=0
                  edd=0
                  
                  xdif=rat(i,iia,1)-rat(j,jja,1)
                 xdif=xdif-anint(xdif/cvec(1,1))*cvec(1,1)
                  ydif=rat(i,iia,2)-rat(j,jja,2)
                 ydif=ydif-anint(ydif/cvec(2,2))*cvec(2,2)
                  zdif=rat(i,iia,3)-rat(j,jja,3)
                 zdif=zdif-anint(zdif/cvec(3,3))*cvec(3,3)
                 
              
              rsq=xdif*xdif+ydif*ydif+zdif*zdif
              if(rsq.le.rcut2) then    !2
                 rsqi=1.d0/rsq
                 rr=dsqrt(rsq)
             iltt=iltt+1 
        
                 ia=1
              if(rr.le.rcutc)then          !3 
               iltt3=iltt3+1
           call smear(dipdip,chgchg,chgdip,rr,rsq,i,j,
     $            iia,jja,ia)
            call field(dipdip,chgchg,chgdip,i,j,iia,jja,ia,
     $           xdif,ydif,zdif,rr,gcfac,dfacx,dfacy,dfacz,
     $           ecc,ecd,edd)
           !print*,rr!,rcutc!ecd,i,j,iia,jja,rr
        endif   !33
     
       if(iflag.eq.1) then  !5
         ucc=ucc+ecc
         ucd=ucd+ecd
         udd=udd+edd
           !print*,ucc,ucd,udd
         !stop
            !if(i.eq.1.and.jja.gt.1)then
           !call lj2(rsqi,rr,qpe,gfac)
           !endif

              if(iia.eq.1.and.jja.eq.1)then !8
               if(i.gt.1)call lj(rsqi,rr,qpe,gfac)
               if(i.eq.1)call lj1(rsqi,rr,qpe,gfac)
            endif       !88
         upair=upair+qpe   
                  
         FORCE(i,iia,1)=FORCE(I,iia,1)+(GFAC+GCFAC)*XDIF+dfacx
         FORCE(i,iia,2)=FORCE(I,iia,2)+(GFAC+GCFAC)*YDIF+dfacy
         FORCE(i,iia,3)=FORCE(I,iia,3)+(GFAC+GCFAC)*ZDIF+dfacz
         FORCE(j,jja,1)=FORCE(J,jja,1)-(GFAC+GCFAC)*XDIF-dfacx
         FORCE(j,jja,2)=FORCE(J,jja,2)-(GFAC+GCFAC)*YDIF-dfacy
         FORCE(j,jja,3)=FORCE(J,jja,3)-(GFAC+GCFAC)*ZDIF-dfacz
         
      endif  !55

      endif   !22
                     
              
         enddo         
         enddo         
         enddo         
         enddo         
        endif
c     and find the dipole forces du/dD

      if(iflag.eq.1) then
         do  i=1,nmol
            if(i.eq.1)i_a=1
            if(i.gt.1)i_a=3
          do iat=1,i_a
                  d2=dipole(i,iat,1)*dipole(i,iat,1)  
     :          +dipole(i,iat,2)*dipole(i,iat,2)
     :        +dipole(i,iat,3)*dipole(i,iat,3)
                  uspring=uspring+fac(4)*d2/(2.d0
     :                 *alph(i,iat))
         enddo
         enddo

      endif
c         print*,iltt3
c       stop
c     evaluate the external field contribution
                 
c      call external(ucx,udx)
              
c     determine the extra forces from the variable charges
c     above those which arise from static charges
                 
           
                 
      
            
      end
      



      subroutine loopy(utot,upair,ucc,ucd,udd,uspring,uintra)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       include 'parameter.i'

       parameter(NMAX=400)
c      parameter(NMol=128)
c      parameter(msp=385,mst=6,mtypemax=100,nspecies=4)

      iflag=0
      iloop=1
      call clear
      call iter
       
      call energy(utot,upair,ucc,ucd,udd,uspring,uintra)
 10   iloop=iloop+1
      call clear
      call iter
          !print*,iflag
      if(iloop.ge.NMAX) then
      write(*,*) 'ERROR.  Maximum number of dipole iterations exceeded.'
      write(*,*) 'You may be using unrealistic geometries.'
      write(*,*) 'If not then convergence may be improved by'
      write(*,*) 'decreasing the dipole mixing parameter newfac.'
      write(*,*)
      write(*,*) 'Current status follows:'
      write(*,*) '(Energies and forces not calculated yet)'
      stop
      endif
      if(iflag.eq.0) goto 10   
      call clear
      call energy(utot,upair,ucc,ucd,udd,uspring,uintra)  

         
     
      end


      subroutine self
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         include 'parameter.i'

    
c     this subroutine calculates the Ewald 'self correction' term.

      sfac=fac(4)*sqpii*ewfac
      do  i=1,nmol
          if(i.eq.1)i_a=1
          if(i.gt.1)i_a=3
         do  ia=1,i_a
            isp=isp+1
         potcc(i,ia)=potcc(i,ia)-2.d0*qat(ia,i)*sfac
         efdd(i,ia,1)=efdd(i,ia,1)+4.d0*ewfac*ewfac*sfac*
     $           (dipole(i,ia,1)+dipole_p(i,ia,1))/3.d0
         efdd(i,ia,2)=efdd(i,ia,2)+4.d0*ewfac*ewfac*sfac*
     $        (dipole(i,ia,2)+dipole_p(i,ia,2))/3.d0
         efdd(i,ia,3)=efdd(i,ia,3)+4.d0*ewfac*ewfac*sfac*
     $        (dipole(i,ia,3)+dipole_p(i,ia,3))/3.d0
         enddo
         enddo
      end


       subroutine ewald
           implicit real*8(a-h,o-z)
       real*8 axes(3,3),raxes(3,3),waxes(3,3)
       real*8 r1(3,3)                  
        include 'parameter.i'
         real*8 ckr(nmol,3),skr(nmol,3),xm(3)


        
cc  ewfac=4.0d0, SQPII/0.56418958354D0
      call recip(vol,voli,rc)

      call fracr(1)

c     this subroutine evaluates the reciprocal space contribution to
c     the ewald sum
      qpec=0.0
c     magnitude of real space and reciprocal space lattice vectors
            
      ak1mag=dsqrt(civec(1,1)**2+civec(1,2)**2+civec(1,3)**2)
      ak2mag=dsqrt(civec(2,1)**2+civec(2,2)**2+civec(2,3)**2)
      ak3mag=dsqrt(civec(3,1)**2+civec(3,2)**2+civec(3,3)**2)
            
      r1mag=dsqrt(cvec(1,1)**2+cvec(1,2)**2+cvec(1,3)**2)
      r2mag=dsqrt(cvec(2,1)**2+cvec(2,2)**2+cvec(2,3)**2)
      r3mag=dsqrt(cvec(3,1)**2+cvec(3,2)**2+cvec(3,3)**2)
      
c     The following comments and looping system
c     are adapted from Smith[]
      
      
C     M1 RANGES OVER THE VALUES 0 TO M1MAX ONLY.
C
C     M2 RANGES OVER 0 TO M2MAX WHEN M1=0 AND OVER
C     -M2MAX TO M2MAX OTHERWISE.
C     M3 RANGES OVER 1 TO M3MAX WHEN M1=M2=0 AND OVER
C     -M3MAX TO M3MAX OTHERWISE.
C     HENCE THE RESULT OF THE SUMMATION MUST BE DOUBLED AT THE END.
C
c     kcut is the reciprocal space cut-off in Ang-1
      akcut=2.5d0
      m1max=int(akcut*r1mag/(2.d0*pi))+1
      m2max=int(akcut*r2mag/(2.d0*pi))+1
      m3max=int(akcut*r3mag/(2.d0*pi))+1
c     Sum over k vectors

            m1min=0
       do m1=m1min,m1max
         if(m1.eq.0) then
            m2min=0
         else
            m2min=-m2max
         endif
     
         do  m2=m2min,m2max
            if((m1.eq.0).and.(m2.eq.0)) then
               m3min=1
            else
               m3min=-m3max   
            endif
      
            do  m3=m3min,m3max
            xm(1)=float(m1)
            xm(2)=float(m2)
            xm(3)=float(m3)
     
            akx=xm(1)*civec(1,1)+xm(2)*civec(2,1)+xm(3)*civec(3,1)
            aky=xm(1)*civec(1,2)+xm(2)*civec(2,2)+xm(3)*civec(3,2)
            akz=xm(1)*civec(1,3)+xm(2)*civec(2,3)+xm(3)*civec(3,3)
            akk2=akx*akx+aky*aky+akz*akz
            if(akk2.le.akcut*akcut) then
c     ckr=cos(kr)
c     cskrs=sum of over j of cos(k.r[j])*chg(j)
c     sum over sites
         
            cckrs=0.d0
            dckrsx=0.d0
            dckrsy=0.d0
            dckrsz=0.d0

            cskrs=0.d0
            dskrsx=0.d0
            dskrsy=0.d0
            dskrsz=0.d0
      
            kdckrs=0.d0
            kdskrs=0.d0
      
            rcckrsx=0.d0
            rcckrsy=0.d0
         rcckrsz=0.d0

            rkdckrsx=0.d0
            rkdckrsy=0.d0
            rkdckrsz=0.d0

            rcskrsx=0.d0
            rcskrsy=0.d0
            rcskrsz=0.d0

            rkdskrsx=0.d0
            rkdskrsy=0.d0
            rkdskrsz=0.d0

               do  i=1,nmol
          if(i.eq.1)i_a=1
          if(i.gt.1)i_a=3

                  do  ia=1,i_a
      thet=2.d0*pi*(xm(1)*rf(i,ia,1)+xm(2)*rf(i,ia,2)+xm(3)*rf(i,ia,3))

      ckr(i,ia)=dcos(thet)
      skr(i,ia)=dsin(thet)
         
c     increment sums
         
      cckr=ckr(i,ia)*qat(ia,i)
      cskr=skr(i,ia)*qat(ia,i)
                  cckrs=cckrs+cckr
                  cskrs=cskrs+cskr
            
            
          dckrsx=dckrsx+ckr(i,ia)*(dipole(i,ia,1)+dipole_p(i,ia,1))
         dckrsy=dckrsy+ckr(i,ia)*(dipole(i,ia,2)+dipole_p(i,ia,2))
        dckrsz=dckrsz+ckr(i,ia)*(dipole(i,ia,3)+dipole_p(i,ia,3))
            
         dskrsx=dskrsx+skr(i,ia)*(dipole(i,ia,1)+dipole_p(i,ia,1))
        dskrsy=dskrsy+skr(i,ia)*(dipole(i,ia,2)+dipole_p(i,ia,2))
        dskrsz=dskrsz+skr(i,ia)*(dipole(i,ia,3)+dipole_p(i,ia,3))
           enddo            
           enddo            
            
c     take dot product of k with the dipole sums
            
               akdckrs=akx*dckrsx+aky*dckrsy+akz*dckrsz
               akdskrs=akx*dskrsx+aky*dskrsy+akz*dskrsz

      aaa=exp(-akk2/(4.d0*ewfac*ewfac))*(fac(4)*4.d0*pi)/(vol*akk2)
          ddd=-2.d0*aaa*((1.d0/akk2)+1.d0/(4.d0*ewfac*ewfac))
                  sum1=(cckrs-akdskrs)
                  sum2=(cskrs+akdckrs)
            
c     increment coulombic energy
            
                  qpec=qpec+aaa*(sum1*sum1+sum2*sum2)
                  ddd=ddd*(sum1*sum1+sum2*sum2)
            
c     sum over sites to increment forces/fields etc.
               do i=1,nmol
          if(i.eq.1)i_a=1
          if(i.gt.1)i_a=3

                do ia=1,i_a            
            
            djk=(dipole(i,ia,1)+dipole_p(i,ia,1))*akx
     $         +(dipole(i,ia,2)+dipole_p(i,ia,2))*aky
     $         +(dipole(i,ia,3)+dipole_p(i,ia,3))*akz
            ss=sum1*(-skr(i,ia)*qat(ia,i)-ckr(i,ia)*djk)
     $           +sum2*(-skr(i,ia)*djk+ckr(i,ia)*qat(ia,i))
            
            tt=-skr(i,ia)*cckrs+ckr(i,ia)*cskrs
            vv=skr(i,ia)*akdskrs+ckr(i,ia)*akdckrs

c     increment EFDC
               if(iloop.eq.1) then
                  efdc(i,ia,1)=efdc(i,ia,1)-2.d0*akx*aaa*tt
                  efdc(i,ia,2)=efdc(i,ia,2)-2.d0*aky*aaa*tt
                  efdc(i,ia,3)=efdc(i,ia,3)-2.d0*akz*aaa*tt
                  
               endif
      
c     increment EFDD
               efdd(i,ia,1)=efdd(i,ia,1)-2.d0*akx*aaa*vv
               efdd(i,ia,2)=efdd(i,ia,2)-2.d0*aky*aaa*vv
               efdd(i,ia,3)=efdd(i,ia,3)-2.d0*akz*aaa*vv
         
      
            if(iflag.eq.1) then
c     increment POTCC
               potcc(i,ia)=potcc(i,ia)+2.d0*aaa*
     $              (ckr(i,ia)*cckrs+skr(i,ia)*cskrs)
c     increment POTCD
               potcd(i,ia)=potcd(i,ia)+2.d0*aaa*
     $              (-ckr(i,ia)*akdskrs+skr(i,ia)*akdckrs)
c     increment forces
               force(i,ia,1)=force(i,ia,1)-2.d0*akx*aaa*ss
               force(i,ia,2)=force(i,ia,2)-2.d0*aky*aaa*ss
               force(i,ia,3)=force(i,ia,3)-2.d0*akz*aaa*ss
            endif      
                  
         enddo            
         enddo            
c     reciprocal space cut-off endif
      endif
            
          enddo          
          enddo          
          enddo          
      end

      subroutine fracr(ifrac)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         include 'parameter.i'
      

c     if(ifrac.eq.0) converts fractional to cartesian coordinates
c     if(ifrac.eq.1) converts cartesian to fractional coordinates
           do i=1,nmol
          if(i.eq.1)i_a=1
          if(i.gt.1)i_a=3

             do ia=1,i_a
              do j=1,3  
               if(ifrac.eq.0) then
      rat(i,ia,j)=cvec(1,j)*rf(i,ia,1)
     $                 +cvec(2,j)*rf(i,ia,2)+cvec(3,j)*rf(i,ia,3)
      else
      rf(i,ia,j)=(civec(j,1)*rat(i,ia,1)+civec(j,2)*rat(i,ia,2)
     $     +civec(j,3)*rat(i,ia,3))/(2.*pi)
         endif
         enddo
         enddo
         enddo
      end

      subroutine recip(vol,voli,rcut)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'parameter.i'     

      !COMMON/constants/debfac,PI,rcut
c     this subroutine evaluates the reciprocal vectors
c     to cvec and also finds the real space cut-off.
      
c     evaluate cross products
      
      
      c12x=cvec(1,2)*cvec(2,3)-cvec(1,3)*cvec(2,2)
      c12y=cvec(1,3)*cvec(2,1)-cvec(1,1)*cvec(2,3)
      c12z=cvec(1,1)*cvec(2,2)-cvec(1,2)*cvec(2,1)
            
      c23x=cvec(2,2)*cvec(3,3)-cvec(2,3)*cvec(3,2)
      c23y=cvec(2,3)*cvec(3,1)-cvec(2,1)*cvec(3,3)
      c23z=cvec(2,1)*cvec(3,2)-cvec(2,2)*cvec(3,1)
      
      c31x=cvec(3,2)*cvec(1,3)-cvec(3,3)*cvec(1,2)
      c31y=cvec(3,3)*cvec(1,1)-cvec(3,1)*cvec(1,3)
      c31z=cvec(3,1)*cvec(1,2)-cvec(3,2)*cvec(1,1)
      
      c12mag=sqrt(c12x*c12x+c12y*c12y+c12z*c12z)
      c23mag=sqrt(c23x*c23x+c23y*c23y+c23z*c23z)
      c31mag=sqrt(c31x*c31x+c31y*c31y+c31z*c31z)

      vol=cvec(1,1)*c23x+cvec(1,2)*c23y+cvec(1,3)*c23z
      voli=1./vol
         
c     evaluate the reciprocal cell vectors
      
      civec(1,1)=2.d0*pi*c23x*voli
      civec(1,2)=2.d0*pi*c23y*voli
      civec(1,3)=2.d0*pi*c23z*voli
      civec(2,1)=2.d0*pi*c31x*voli
      civec(2,2)=2.d0*pi*c31y*voli
      civec(2,3)=2.d0*pi*c31z*voli
            
      civec(3,1)=2.d0*pi*c12x*voli  
      civec(3,2)=2.d0*pi*c12y*voli   
      civec(3,3)=2.d0*pi*c12z*voli   
            
c     find 3 radii for determination of real space cut-off
         
      r12=0.5d0*vol/c12mag
      r23=0.5d0*vol/c23mag   
      r31=0.5d0*vol/c31mag
         
c     the cut-off is the smallest of the 3 radii
      
      rcut=r12
      if(r23.ge.r12) then
         if(r31.le.r12) rcut=r31
         else
            rcut=r23
            if(r31.le.r23) rcut=r31
            endif   
            
      end



      subroutine iter
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         include 'parameter.i'     
        toler=.10D-6
        oldfac=0.5d0
         !print*,'iter'
      fac4i=1.d0/fac(4)
        debfac=echarge*1d-10/3.33d-30
      if(iloop.eq.1) then
         goto 9000
      endif
      
      call ddfield
      deltadip=0.d0
      
      do  i=1,nmol
          if(i.eq.1)i_a=1
          if(i.gt.1)i_a=3

          do ia=1,i_a
               oldipx=dipole(i,ia,1)
               oldipy=dipole(i,ia,2)
               oldipz=dipole(i,ia,3)
         
c     mixing ratio
c     d(n+1)=newfac*alpha(E+T.d(n))+oldfac*d(n-1)
               dnewfac=1.-oldfac
        dipole(i,ia,1)=oldfac*oldipx+(dnewfac*alph(i,ia)*fac4i*
     $       (efdd(i,ia,1)+efdc(i,ia,1)+efext(1)))
        dipole(i,ia,2)=oldfac*oldipy+(dnewfac*alph(i,ia)*fac4i*
     $       (efdd(i,ia,2)+efdc(i,ia,2)+efext(2)))   
        dipole(i,ia,3)=oldfac*oldipz+(dnewfac*alph(i,ia)*fac4i*
     $       (efdd(i,ia,3)+efdc(i,ia,3)+efext(3)))
        dipmag=dsqrt(dipole(i,ia,1)**2+dipole(i,ia,2)**2
     : +dipole(i,ia,3)**2)
        deltadip=deltadip+((dipole(i,ia,1)-oldipx)**2.d0
     a                    +(dipole(i,ia,2)-oldipy)**2.d0
     a                    +(dipole(i,ia,3)-oldipz)**2.d0)
         !print*,dipmag*debfac
          !print*,efdc(i,ia,1) 
c           print*,dipole(1,1,1)*debfac,oldipx
c           print*,'field'
c          print*,
c     $       (efdd(i,ia,1))!*debfac         

c          print*,
c     $       (efdc(i,ia,1))!*debfac

c         stop         
c     dipmol=sqrt(dipole(isp,1)**2+dipole(isp,2)**2+dipole(isp,3)**2)
c      write(884,*) isp,dipole(isp,1)*debfac
        enddo      
        enddo      
      deltadip=dsqrt(deltadip/dfloat(nmol))*debfac
c         print*,deltadip      
c       stop
      iflag=0
      if(deltadip.gt.1.e+10) then
        print*,'Trouble convering dipoles'
         stop
      endif
      if(deltadip.le.toler) then 
         iflag=1
      endif
 
 9000 continue
         
      end  


      subroutine clear 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         include 'parameter.i'      

      
      dpdr=0.0
         isp=0
         do i=1,nmol
          if(i.eq.1)i_a=1
          if(i.gt.1)i_a=3

               do ia=1,i_a
               potcc(i,ia)=0.d0
               potcd(i,ia)=0.d0
               pot(i,ia)=0.d0
c     set initial value of charges
               !qat(ia,i)=qatd(ia,i)
         do ivec=1,3
            efdd(i,ia,ivec)=0.d0
            if(iloop.eq.1) then
               efdc(i,ia,ivec)=0.d0
               force(i,ia,ivec)=0.d0
               endif
        enddo
        enddo
        enddo
     
      end

      subroutine energy(utot,upair,ucc,ucd,udd,uspring,uintra)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         include 'parameter.i'      

      
      call real(utot,upair,ucc,ucd,udd,uspring,uintra,ucx,udx)
         if(ipdc.eq.1) then
         if(iflag.eq.1) then
            call ewald
            call self
         endif
         endif
      
      call ddfield
      if(iflag.eq.1) call dipforce
      
      call coulomb(ucc,ucd,udd,uspring)
      utot=upair+ucc+ucd+udd+uspring+ucx+udx+uintra
c$$$         write(*,*) 'utot=',utot*pefac
c$$$         write(*,*) 'upair=',upair*pefac
c$$$         write(*,*) 'uintra=',uintra*pefac
c$$$         write(*,*) 'ucc=',ucc*pefac
c$$$         write(*,*) 'ucd=',ucd*pefac
c$$$         write(*,*) 'udd=',udd*pefac
      end



c      subroutine printparam(a)
c          include 'parameter.i'
c      print*,alph(1,1),c_16,c_14,c_12,c_6
c      print*,alph(1,1),c_16,c_14,c_12,c_6
c      end

      subroutine lj1(rsqi,rr,qpe,gfac)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          include 'parameter.i'
            RR6=RSQI**3.d0
            RR8=RSQI**4.d0
            RR12=RR6*RR6
            RR14=RR6*RR8
            RR16=RR8*RR8
            rr18=rr16*rsqi
            RR10=RSQI**5.d0
            RR20=RR10*RR10
            GFAC=(16.d0*c_16*RR18+14.d0*c_14*RR16
     :    +12.d0*c_12*RR14
     $           +6.d0*c_6*RR8)
            qpe=(c_16*RR16+c_14*RR14+c_12*RR12+c_6*RR6)
           !print*,qpe
         end

      subroutine lj(rsqi,rr,qpe,gfac)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          include 'parameter.i'              
            RR6=RSQI**3.0
            RR8=RSQI**4.0
            RR12=RR6*RR6  
            RR14=RR6*RR8 
            RR16=RR8*RR8  
            rr18=rr16*rsqi
            RR10=RSQI**5.0
            RR20=RR10*RR10
            
            GFAC=(16.d0*c16*RR18+14.d0*c14*RR16+12.d0*c12*RR14
     $           +6.d0*c6*RR8)
            qpe=(c16*RR16+c14*RR14+c12*RR12+c6*RR6)
          ! print*,qpe
         end

      
      subroutine dip2_h2o(r1x,r1y,r1z,r2x,r2y,r2z,r1mag,r2mag
     $     ,p1,p2,p3
     $     ,gradq)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      real*8 x(1,3,3),d(1,3),gradq(3,3,3),gradqt(3,3,3)
      real*8 gradqr(3,3,3),gradqct(3,3,3)
         include 'parameter.i'
      
c     gradq are the charge gradients in Cartesians
c     gradq(i,j,k) is the gradient with respect to the
c     ith nuclear coordinate of the jth charge in the
c     kth direction.
c     i=1,2,3=H1,H2,O
c     j=1,2,3=qH1,qH2,qO
c     k=1,2,3=x,y,z
      
      r0=0.9572d0
      c1= 0.12d0
      
      dp1dr1=c1
      dp1dr2=c1
      dp2dr1=c1
      dp2dr2=c1

c     determine gradients of the charges
      
      f1q1r13=(dp1dr1)/r1mag
      f1q1r23=0.0
      f2q1r23=(dp1dr2)/r2mag
      f2q1r13=0.0
      f1q2r13=(dp2dr1)/r1mag
      f1q2r23=0.0
      f2q2r23=(dp2dr2)/r2mag
      f2q2r13=0.0

c     gradient of charge h1 wrt displacement of h1
      gradqr(1,1,1)=f1q1r13*r1x+f1q1r23*r2x
      gradqr(1,1,2)=f1q1r13*r1y+f1q1r23*r2y
      gradqr(1,1,3)=f1q1r13*r1z+f1q1r23*r2z

c     gradient of charge h1 wrt displacement of h2
      gradqr(2,1,1)=f2q1r13*r1x+f2q1r23*r2x
      gradqr(2,1,2)=f2q1r13*r1y+f2q1r23*r2y
      gradqr(2,1,3)=f2q1r13*r1z+f2q1r23*r2z
     
c     gradient of charge h1 wrt displacement of O
      gradqr(3,1,1)=-(gradqr(1,1,1)+gradqr(2,1,1))
      gradqr(3,1,2)=-(gradqr(1,1,2)+gradqr(2,1,2))
      gradqr(3,1,3)=-(gradqr(1,1,3)+gradqr(2,1,3))
     
c     gradient of charge h2 wrt displacement of h1
      gradqr(1,2,1)=f1q2r13*r1x+f1q2r23*r2x
      gradqr(1,2,2)=f1q2r13*r1y+f1q2r23*r2y
      gradqr(1,2,3)=f1q2r13*r1z+f1q2r23*r2z
      
c     gradient of charge h2 wrt displacement of h2
      gradqr(2,2,1)=f2q2r13*r1x+f2q2r23*r2x
      gradqr(2,2,2)=f2q2r13*r1y+f2q2r23*r2y
      gradqr(2,2,3)=f2q2r13*r1z+f2q2r23*r2z
      
c     gradient of charge h2 wrt displacement of O
      gradqr(3,2,1)=-(gradqr(1,2,1)+gradqr(2,2,1))
      gradqr(3,2,2)=-(gradqr(1,2,2)+gradqr(2,2,2))
      gradqr(3,2,3)=-(gradqr(1,2,3)+gradqr(2,2,3))
      
     
c     gradient of charge O wrt displacement of h1
      gradqr(1,3,1)=-(gradqr(1,1,1)+gradqr(1,2,1))
      gradqr(1,3,2)=-(gradqr(1,1,2)+gradqr(1,2,2))
      gradqr(1,3,3)=-(gradqr(1,1,3)+gradqr(1,2,3))
      
c     gradient of charge O wrt displacement of h2
      gradqr(2,3,1)=-(gradqr(2,1,1)+gradqr(2,2,1))
      gradqr(2,3,2)=-(gradqr(2,1,2)+gradqr(2,2,2))
      gradqr(2,3,3)=-(gradqr(2,1,3)+gradqr(2,2,3))

c     gradient of charge O wrt displacement of O
      gradqr(3,3,1)=-(gradqr(3,1,1)+gradqr(3,2,1))
      gradqr(3,3,2)=-(gradqr(3,1,2)+gradqr(3,2,2))
      gradqr(3,3,3)=-(gradqr(3,1,3)+gradqr(3,2,3))

      
      theta0=104.52d0*pi/180.00
      c2=(0.00473d0*180.d0/pi)*0.d0
      r1r2i=(1.d0/(r1mag*r2mag))
      r1dotr2=(r1x*r2x+r1y*r2y+r1z*r2z)
      costh=(r1dotr2)*r1r2i
      theta=acos(costh)
      
c     write(*,*)r1r2i,costh,theta,pi

      fac1=(c2*r1r2i)/(-sin(theta))
      
      x1fac=(r2x)-(r1dotr2*r1x/(r1mag*r1mag))
      y1fac=(r2y)-(r1dotr2*r1y/(r1mag*r1mag))
      z1fac=(r2z)-(r1dotr2*r1z/(r1mag*r1mag))
      
      x2fac=(r1x)-(r1dotr2*r2x/(r2mag*r2mag))
      y2fac=(r1y)-(r1dotr2*r2y/(r2mag*r2mag))
      z2fac=(r1z)-(r1dotr2*r2z/(r2mag*r2mag))

c     gradient of charge h1 wrt displacement of h1
      gradqt(1,1,1)=fac1*x1fac
      gradqt(1,1,2)=fac1*y1fac
      gradqt(1,1,3)=fac1*z1fac

c     gradient of charge h1 wrt displacement of h2
      gradqt(2,1,1)=fac1*x2fac
      gradqt(2,1,2)=fac1*y2fac
      gradqt(2,1,3)=fac1*z2fac
     
c     gradient of charge h1 wrt displacement of O
      gradqt(3,1,1)=-(gradqt(1,1,1)+gradqt(2,1,1))
      gradqt(3,1,2)=-(gradqt(1,1,2)+gradqt(2,1,2))
      gradqt(3,1,3)=-(gradqt(1,1,3)+gradqt(2,1,3))
     
c     gradient of charge h2 wrt displacement of h1
      gradqt(1,2,1)=gradqt(1,1,1)
      gradqt(1,2,2)=gradqt(1,1,2)
      gradqt(1,2,3)=gradqt(1,1,3)
      
c     gradient of charge h2 wrt displacement of h2
      gradqt(2,2,1)=gradqt(2,1,1)
      gradqt(2,2,2)=gradqt(2,1,2)
      gradqt(2,2,3)=gradqt(2,1,3)
      
c     gradient of charge h2 wrt displacement of O
      gradqt(3,2,1)=-(gradqt(1,2,1)+gradqt(2,2,1))
      gradqt(3,2,2)=-(gradqt(1,2,2)+gradqt(2,2,2))
      gradqt(3,2,3)=-(gradqt(1,2,3)+gradqt(2,2,3))
      
c     gradient of charge O wrt displacement of h1
      gradqt(1,3,1)=-(gradqt(1,1,1)+gradqt(1,2,1))
      gradqt(1,3,2)=-(gradqt(1,1,2)+gradqt(1,2,2))
      gradqt(1,3,3)=-(gradqt(1,1,3)+gradqt(1,2,3))
      
c     gradient of charge O wrt displacement of h2
      gradqt(2,3,1)=-(gradqt(2,1,1)+gradqt(2,2,1))
      gradqt(2,3,2)=-(gradqt(2,1,2)+gradqt(2,2,2))
      gradqt(2,3,3)=-(gradqt(2,1,3)+gradqt(2,2,3))
      
c     gradient of charge O wrt displacement of O
      gradqt(3,3,1)=-(gradqt(3,1,1)+gradqt(3,2,1))
      gradqt(3,3,2)=-(gradqt(3,1,2)+gradqt(3,2,2))
      gradqt(3,3,3)=-(gradqt(3,1,3)+gradqt(3,2,3))
      
      c3=-.1360000000d0
      
      dp1dr1 = c3
      dp1dr2 = c3
      dp2dr1 = c3
      dp2dr2 = c3
      
      f1q1r13=(dp1dr1)/r1mag
      f1q1r23=0.0
      f2q1r23=(dp1dr2)/r2mag
      f2q1r13=0.0
      f1q2r13=(dp2dr1)/r1mag
      f1q2r23=0.0
      f2q2r23=(dp2dr2)/r2mag
      f2q2r13=0.0
      
c     gradient of charge h1 wrt displacement of h1
      gradqct(1,1,1)=f1q1r13*r1x-f1q1r23*r2x 
      gradqct(1,1,2)=f1q1r13*r1y-f1q1r23*r2y 
      gradqct(1,1,3)=f1q1r13*r1z-f1q1r23*r2z 

c     gradient of charge h1 wrt displacement of h2
      gradqct(2,1,1)=f2q1r13*r1x-f2q1r23*r2x
      gradqct(2,1,2)=f2q1r13*r1y-f2q1r23*r2y
      gradqct(2,1,3)=f2q1r13*r1z-f2q1r23*r2z

c     gradient of charge h1 wrt displacement of O 
      gradqct(3,1,1)=-(gradqct(1,1,1)+gradqct(2,1,1))
      gradqct(3,1,2)=-(gradqct(1,1,2)+gradqct(2,1,2))
      gradqct(3,1,3)=-(gradqct(1,1,3)+gradqct(2,1,3))
     
c     gradient of charge h2 wrt displacement of h1
      gradqct(1,2,1)=-gradqct(1,1,1)
      gradqct(1,2,2)=-gradqct(1,1,2)
      gradqct(1,2,3)=-gradqct(1,1,3)
     
c     gradient of charge h2 wrt displacement of h2
      gradqct(2,2,1)=-gradqct(2,1,1)
      gradqct(2,2,2)=-gradqct(2,1,2)
      gradqct(2,2,3)=-gradqct(2,1,3)
      
c     gradient of charge h2 wrt displacement of O 
      gradqct(3,2,1)=-gradqct(3,1,1)
      gradqct(3,2,2)=-gradqct(3,1,2)
      gradqct(3,2,3)=-gradqct(3,1,3)
      
c     gradient of charge O wrt displacement of h1
      gradqct(1,3,1)=-(gradqct(1,1,1)+gradqct(1,2,1))
      gradqct(1,3,2)=-(gradqct(1,1,2)+gradqct(1,2,2))
      gradqct(1,3,3)=-(gradqct(1,1,3)+gradqct(1,2,3))
      
c     gradient of charge O wrt displacement of h2
      gradqct(2,3,1)=-(gradqct(2,1,1)+gradqct(2,2,1))
      gradqct(2,3,2)=-(gradqct(2,1,2)+gradqct(2,2,2))
      gradqct(2,3,3)=-(gradqct(2,1,3)+gradqct(2,2,3))
      
c     gradient of charge O wrt displacement of O 
      gradqct(3,3,1)=-(gradqct(3,1,1)+gradqct(3,2,1))
      gradqct(3,3,2)=-(gradqct(3,1,2)+gradqct(3,2,2))
      gradqct(3,3,3)=-(gradqct(3,1,3)+gradqct(3,2,3))
      
      gradq(1,1,1)= gradqr(1,1,1)+ gradqt(1,1,1)+ gradqct(1,1,1)
      gradq(1,1,2)= gradqr(1,1,2)+ gradqt(1,1,2)+ gradqct(1,1,2)
      gradq(1,1,3)= gradqr(1,1,3)+ gradqt(1,1,3)+ gradqct(1,1,3)
      
      gradq(2,1,1)= gradqr(2,1,1)+ gradqt(2,1,1)+ gradqct(2,1,1)
      gradq(2,1,2)= gradqr(2,1,2)+ gradqt(2,1,2)+ gradqct(2,1,2)
      gradq(2,1,3)= gradqr(2,1,3)+ gradqt(2,1,3)+ gradqct(2,1,3)
      
      gradq(3,1,1)= gradqr(3,1,1)+ gradqt(3,1,1)+ gradqct(3,1,1)
      gradq(3,1,2)= gradqr(3,1,2)+ gradqt(3,1,2)+ gradqct(3,1,2)
      gradq(3,1,3)= gradqr(3,1,3)+ gradqt(3,1,3)+ gradqct(3,1,3)
      
      gradq(1,2,1)= gradqr(1,2,1)+ gradqt(1,2,1)+ gradqct(1,2,1)
      gradq(1,2,2)= gradqr(1,2,2)+ gradqt(1,2,2)+ gradqct(1,2,2)
      gradq(1,2,3)= gradqr(1,2,3)+ gradqt(1,2,3)+ gradqct(1,2,3)
      
      gradq(2,2,1)= gradqr(2,2,1)+ gradqt(2,2,1)+ gradqct(2,2,1)
      gradq(2,2,2)= gradqr(2,2,2)+ gradqt(2,2,2)+ gradqct(2,2,2)
      gradq(2,2,3)= gradqr(2,2,3)+ gradqt(2,2,3)+ gradqct(2,2,3)
      
      gradq(3,2,1)= gradqr(3,2,1)+ gradqt(3,2,1)+ gradqct(3,2,1)
      gradq(3,2,2)= gradqr(3,2,2)+ gradqt(3,2,2)+ gradqct(3,2,2)
      gradq(3,2,3)= gradqr(3,2,3)+ gradqt(3,2,3)+ gradqct(3,2,3)
     
      gradq(1,3,1)= gradqr(1,3,1)+ gradqt(1,3,1)+ gradqct(1,3,1)
      gradq(1,3,2)= gradqr(1,3,2)+ gradqt(1,3,2)+ gradqct(1,3,2)
      gradq(1,3,3)= gradqr(1,3,3)+ gradqt(1,3,3)+ gradqct(1,3,3)
      
      gradq(2,3,1)= gradqr(2,3,1)+ gradqt(2,3,1)+ gradqct(2,3,1)
      gradq(2,3,2)= gradqr(2,3,2)+ gradqt(2,3,2)+ gradqct(2,3,2)
      gradq(2,3,3)= gradqr(2,3,3)+ gradqt(2,3,3)+ gradqct(2,3,3)

      gradq(3,3,1)= gradqr(3,3,1)+ gradqt(3,3,1)+ gradqct(3,3,1)
      gradq(3,3,2)= gradqr(3,3,2)+ gradqt(3,3,2)+ gradqct(3,3,2)
      gradq(3,3,3)= gradqr(3,3,3)+ gradqt(3,3,3)+ gradqct(3,3,3)
     
      c3fach1=c3*(r1mag-r2mag)
      c3fach2=-c3fach1
      
      p1=0.42d0+c1*(r1mag+r2mag-2.*r0)+c2*(theta-theta0)
     $     +c3fach1
      p2=0.42d0+c1*(r1mag+r2mag-2.*r0)+c2*(theta-theta0)
     $     +c3fach2
      
      p3=-(p1+p2)
      
c      write(*,*)p1,p2,p3, (p1+p2+p3)
c      write(*,*)r1mag,r2mag
      
      end


       subroutine dipforce
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
           include 'parameter.i'      
               do i=1,nmol

             if(iand.eq.1)then
         potisp1=potcc(i,2)+potcd(i,2)
         potisp2=potcc(i,3)+potcd(i,3)
         potisp3=potcc(i,1)+potcd(i,1)
            else  
         potisp1=potcd(i,2)
         potisp2=potcd(i,3)
         potisp3=potcd(i,1)
            endif


         do ivec=1,3
       FORCE(I,2,ivec)=FORCE(I,2,ivec)-(fq(i,1,1,ivec)*potisp1
     $                                   +fq(i,1,2,ivec)*potisp2
     $                                   +fq(i,1,3,ivec)*potisp3)

       FORCE(i,3,ivec)=FORCE(i,3,ivec)-(fq(i,2,1,ivec)*potisp1
     $                                   +fq(i,2,2,ivec)*potisp2
     $                                   +fq(i,2,3,ivec)*potisp3)
 
      
       FORCE(i,1,ivec)=FORCE(i,1,ivec)-(fq(i,3,1,ivec)*potisp1
     $                                 +fq(i,3,2,ivec)*potisp2
     $                                 +fq(i,3,3,ivec)*potisp3)
              
        enddo
        enddo
 
       
      end



      subroutine subtract(utot,upair,ucc,ucd,udd,uspring,uintra
     $     ,ucx,udx,iton)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         include 'parameter.i'
c         include 'parameter2.i'
              pi=acos(-1.d0)      
                  epso=8.854187817D-12
            esconvert=1.d0/(4.0d0*pi*epso)
         fac(4)=echarge**2*esconvert
           fac(4)=fac(4)*1d+10*1.439325215D+20
        debfac=echarge*1d-10/3.33d-30

              !alph(1,1)=0.1d0
            do i=2,nmol
            alph(i,1)=1.365d0
            alph(i,2)=0.471d0
            alph(i,3)=0.471d0
             enddo

      do  i=1,nmol
           do ia=1,3
         do  j=1,nmol
          do ja=1,3
      diptens1(i,j,ia,ja)=0.d0
      diptens2(i,j,ia,ja)=0.d0
      diptens3(i,j,ia,ja)=0.d0
      diptens4(i,j,ia,ja)=0.d0
      diptens5(i,j,ia,ja)=0.d0
      diptens6(i,j,ia,ja)=0.d0
        enddo
        enddo
        enddo
        enddo
          do i=1,nmol
             do i1=1,3
             do i2=1,3
             do i3=1,3
             do i4=1,3
             dpdr(i,i1,i2,i3,i4)=0.d0
         enddo     
         enddo     
         enddo     
         enddo     
         enddo     

      iand=1!2
      ipdc=1!2
      call loopy(utot,upair,ucc,ucd,udd,uspring,uintra)
           sumdx=debfac*dipole(1,1,1)
           sumdy=debfac*dipole(1,1,2)
           sumdz=debfac*dipole(1,1,3)
         do i=2,nmol
          do iat=1,3
          do idim=1,3
          xd=rat(i,iat,idim)-rat(i,1,idim)
            xd=xd-anint(xd/cvec(1,1))*cvec(1,1)
           if(idim.eq.1)sumdx=sumdx+
     :qat(iat,i)*xd*debfac+
     :(dipole_p(i,iat,idim)+dipole(i,iat,idim))*debfac
           if(idim.eq.2)sumdy=sumdy+
     :qat(iat,i)*xd*debfac+
     :(dipole_p(i,iat,idim)+dipole(i,iat,idim))*debfac
           if(idim.eq.3)sumdz=sumdz+
     :qat(iat,i)*xd*debfac+
     :(dipole_p(i,iat,idim)+dipole(i,iat,idim))*debfac
         enddo
         enddo
         enddo


      utot1=utot
      upair1=upair  
      ucc1=ucc
 
      ucd1=ucd
      udd1=udd
      uspring1=uspring
      uintra1=uintra
      ucx1=ucx
      udx1=udx
       !print*,utot,upair,ucc,ucd,udd,uspring,uintra
        !print*,force(1,1,1)
         !stop
      do i=1,nmol
           do ia=1,3
         do  ivec=1,3
            force1(i,ia,ivec)=force(i,ia,ivec)
            dforce1(i,ia,ivec)=dforce(i,ia,ivec)
            diporig(i,ia,ivec)=dipole(i,ia,ivec)
       enddo
       enddo
       enddo
        !print*,'ok2'
      do  i=1,nmol
           do ia=1,3
         do  j=1,nmol
          do ja=1,3
      diptens1(i,j,ia,ja)=0.d0
      diptens2(i,j,ia,ja)=0.d0
      diptens3(i,j,ia,ja)=0.d0
      diptens4(i,j,ia,ja)=0.d0
      diptens5(i,j,ia,ja)=0.d0
      diptens6(i,j,ia,ja)=0.d0
        enddo
        enddo
        enddo
        enddo
          do i=1,nmol
             do i1=1,3
             do i2=1,3
             do i3=1,3
             do i4=1,3
             dpdr(i,i1,i2,i3,i4)=0.d0
         enddo     
         enddo     
         enddo     
         enddo     
         enddo     

      iand=2
      iton_orig=iton
      iton=1
      ipdc=2
      call loopy(utot,upair,ucc,ucd,udd,uspring,uintra)

      iton=iton_orig
      utot2=utot
      upair2=upair
      ucc2=ucc
      ucd2=ucd
      udd2=udd
      uspring2=uspring
      uintra2=uintra
      ucx2=ucx
      udx2=udx
      
      do i=1,nmol
         do ia=1,3
         do  ivec=1,3
            force2(i,ia,ivec)=force(i,ia,ivec)
            dforce2(i,ia,ivec)=dforce(i,ia,ivec)
       enddo
       enddo
       enddo
      
      do  i=1,nmol
      do  ia=1,3
         do  ivec=1,3
            force(i,ia,ivec)=force1(i,ia,ivec)-force2(i,ia,ivec)
            dforce(i,ia,ivec)=dforce1(i,ia,ivec)-dforce2(i,ia,ivec)
            dipole(i,ia,ivec)=diporig(i,ia,ivec)
        enddo
        enddo
        enddo

      utot=utot1-utot2   
      upair=upair1-upair2
      ucc=ucc1-ucc2
      ucd=ucd1-ucd2
      udd=udd1-udd2
      uspring=uspring1-uspring2
      uintra=uintra1-uintra2
      ucx=ucx1-ucx2
      udx=udx1-udx2
     
c     utot=utot-float(nmolw)*2.94074012287069/pefac
c        print*,utot!,force(128,1,3)
c        stop
        print*,'utot,ucc,ucd,udd,uspring,upair'!,force(128,1,3)
        print*,utot,ucc,ucd,udd,uspring,upair!,force(128,1,3)
      end


      subroutine ddfield
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         include 'parameter.i'      
      if(ipdc.eq.1) then
      
         if(iflag.ne.1) then
            call ewald
            call self
            endif
         endif
           isp=0
      do i=1,nmol
          if(i.eq.1)i_a=1
          if(i.gt.1)i_a=3

          do ia=1,i_a
            isp=isp+1
           jsp=0

      do j=1,nmol
          if(j.eq.1)j_a=1
          if(j.gt.1)j_a=3

          do ja=1,j_a
            jsp=jsp+1
         
            if(jsp.gt.isp) then
               dtens1=diptens1(i,j,ia,ja)
               dtens2=diptens2(i,j,ia,ja) 
               dtens3=diptens3(i,j,ia,ja)
               dtens4=diptens4(i,j,ia,ja)
               dtens5=diptens5(i,j,ia,ja)
               dtens6=diptens6(i,j,ia,ja)
               di1=dipole(i,ia,1)+dipole_p(i,ia,1)
               di2=dipole(i,ia,2)+dipole_p(i,ia,2)
               di3=dipole(i,ia,3)+dipole_p(i,ia,3)
               dj1=dipole(j,ja,1)+dipole_p(j,ja,1)
               dj2=dipole(j,ja,2)+dipole_p(j,ja,2)
               dj3=dipole(j,ja,3)+dipole_p(j,ja,3)
      EFDD(i,ia,1)=EFDD(i,ia,1)+dtens1*dj1
     $                           +dtens2*dj2
     $                           +dtens3*dj3
      EFDD(i,ia,2)=EFDD(i,ia,2)+dtens2*dj1
     $                           +dtens4*dj2
     $                           +dtens5*dj3
      EFDD(i,ia,3)=EFDD(i,ia,3)+dtens3*dj1
     $                           +dtens5*dj2
     $                           +dtens6*dj3
      EFDD(j,ja,1)=EFDD(j,ja,1)+dtens1*di1
     $                           +dtens2*di2
     $                           +dtens3*di3
      EFDD(j,ja,2)=EFDD(j,ja,2)+dtens2*di1
     $                           +dtens4*di2
     $                           +dtens5*di3
      EFDD(j,ja,3)=EFDD(j,ja,3)+dtens3*di1 
     $                           +dtens5*di2
     $                           +dtens6*di3



            endif
          enddo      
          enddo      
          enddo      
          enddo   
            isp=0
            do i=1,nmol
          if(i.eq.1)i_a=1
          if(i.gt.1)i_a=3

               do ia=1,i_a
                  isp=isp+1
         dforce(i,ia,1)=-(efdc(i,ia,1)+efdd(i,ia,1)+efext(1))
     $                  +fac(4)*dipole(i,ia,1)/alph(i,ia)
         dforce(i,ia,2)=-(efdc(i,ia,2)+efdd(i,ia,2)+efext(2))
     $                  +fac(4)*dipole(i,ia,2)/alph(i,ia)
         dforce(i,ia,3)=-(efdc(i,ia,3)+efdd(i,ia,3)+efext(3))
     $                  +fac(4)*dipole(i,ia,3)/alph(i,ia)
          enddo
          enddo
            end

         subroutine smear(dipdip,chgchg,chgdip,rr,rsq,i,j,
     $            iia,jja,ia)  

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         real*8 phi(0:3),bbb(0:3)
         include 'parameter.i'

c     tabulation of the 'Smith' multipole functions
c     for smeared dipole-dipole, charge-dipole and charge-charge
c     interactions [7]
      rsqi=1./rsq
      rr4=rsq*rsq
      rr3i=rsqi/rr   
      rr4i=rsqi*rsqi
      rri=rsqi*rr
         
      if(ipdc.eq.1) then
c     calculate the 'Smith' B (error) functions
c     for real space Ewald sum.  (I've modified
c     Smith's recursion to give -erf, not erfc functions.)
        
         a=ewfac*rr
         bbb(0)=(erfcc(a)-1.d0)*rri
         exp2a=exp(-a*a)

         do  m=1,3   
            fm=float(m)
            bbb(m)=rsqi*((2.d0*fm-1.d0)*bbb(m-1)
     :           +((2.0*ewfac**2.0)**fm)*sqpii*exp2a/ewfac)
         enddo
       endif
         
             
         alphai=alph(i,iia)
         alphaj=alph(j,jja)
         
            
            
c     calculate the charge-charge interactions
      if(ia.ne.0) then
c     i.e. if this is an intermolecular interaction calculate the
c     charge-charge interactions
            afac=.50d0
            call tholesmear(phi, 
     $           afac,alphai,alphaj,rr,rri,rsqi,rr3i,rr4i,4)

        chgchg(0)=phi(0)
        chgchg(1)=-phi(1)*rri
        chgchg(2)=(phi(2)-phi(1)*rri)*rsqi
        chgchg(3)=(3.0*(phi(1)*rri
     $       -phi(2))+phi(3)*rr)*(-rr4i)
      else
c     zero intramolecular charge-charge interactions
         chgchg(0)=0.d0
         chgchg(1)=0.d0
         chgchg(2)=0.d0
         chgchg(3)=0.d0
      endif
     
c     Calculate the charge-dipole interactions
      
c     i.e. if this is an intermolecular interaction or if this is POLIR's model (inter or intra) 
calculate the
c     charge-dipole interactions 
            afac=0.15d0
           if(i.eq.1)afac=0.1d0
            if(ia.eq.0) then
               if(iia.eq.1.and.jja.gt.1) afac=0.65d0000
               if(iia.gt.1.and.jja.eq.1) afac=0.65d0000
               if(iia.gt.1.and.jja.gt.1) afac=0.05d0
            endif   
            call tholesmear(phi,
     $           afac,alphai,alphaj,rr,rri,rsqi,rr3i,rr4i,4)

         chgdip(0)=phi(0)
         chgdip(1)=-phi(1)*rri
         chgdip(2)=(phi(2)-phi(1)*rri)*rsqi
         chgdip(3)=(3.0*(phi(1)*rri
     $        -phi(2))+phi(3)*rr)*(-rr4i)
            

c     Calculate the dipole-dipole interactions
               afac=0.3d0
           !if(i.eq.1)afac=0.001d0
               if(ia.eq.0) then
                  if(iia.eq.1.and.jja.gt.1) afac=0.69000d0
                  if(iia.gt.1.and.jja.eq.1) afac=0.69000d0
                  if(iia.gt.1.and.jja.gt.1) afac=0.05d0
               endif  
               call tholesmear(phi,
     $              afac,alphai,alphaj,rr,rri,rsqi,rr3i,rr4i,4)
            dipdip(0)=phi(0)
            dipdip(1)=-phi(1)*rri
            dipdip(2)=(phi(2)-phi(1)*rri)*rsqi
            dipdip(3)=(3.0*(phi(1)*rri
     $           -phi(2))+phi(3)*rr)*(-rr4i)
            
         
            
      if(ipdc.eq.1) then
         do  ii=0,3
            chgchg(ii)=chgchg(ii)+bbb(ii)
            chgdip(ii)=chgdip(ii)+bbb(ii)
            dipdip(ii)=dipdip(ii)+bbb(ii)
        enddo 
       endif   

              
      end

      subroutine tholesmear(phi,
     $     afac,alphai,alphaj,rr,rri,rsqi,rr3i,rr4i,m)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      real*8 phi(0:3)
         include 'parameter.i'
c     Evaluation of derivatives of the Thole damping term
            
      xm=float(m)
      cc=float(m-1)/float(m)   
      g34=dexp((gammln(cc)))
      aiaj6=(alphai*alphaj)**(1.d0/6.d0)
      sss=afac*aiaj6  
      ra=rr/aiaj6
      ra4=ra**xm
      afra4=afac*(ra**xm)
      efac=exp(-afra4)
       !print*,cc,afra4,rr

      phi(0)=(1.0-efac
     $     +(afra4**(1.0/xm))*gammq(cc,afra4)*g34)*rri
      phi(1)=-rsqi*(1.0-efac)
      phi(2)=(-xm*ra4*afac*efac
     $     +2.0*(1.0-efac))*rr3i
      phi(3)=(-6.0*(1.0-efac)
     $     +(5*xm-xm*xm)*afac*ra4*efac
     $     +xm*xm*afac*afac*ra4*ra4*efac)*rr4i
      end
            
      function gammq(a,x)
      IMPLICIT real*8 (A-H,O-Z)
      real*8 a,gammq,x
      real*8 gammcf,gamser,gln
      if(x.lt.0.d0.or.a.le.0.d0)pause 'bad arguments in gammq'
      if(x.lt.a+1.)then
         call gser(gamser,a,x,gln)
         gammq=1.-gamser
      else
         call gcf(gammcf,a,x,gln)
         gammq=gammcf
      endif
      return
      END
         
      SUBROUTINE gser(gamser,a,x,gln)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER ITMAX
      REAL*8 a,gamser,gln,x,EPS
      PARAMETER (ITMAX=800,EPS=3.e-7)
      INTEGER n   
      REAL*8 ap,del,sum,gammln
      gln=gammln(a)   
      if(x.le.0.) then
         if(x.lt.0.)pause 'x < 0 in gser'
         gamser=0.
         return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
         ap=ap+1.   
         del=del*x/ap
         sum=sum+del
         if(abs(del).lt.abs(sum)*EPS)goto 1
 11   continue
      pause 'a too large, ITMAX too small in gser'
 1    gamser=sum*exp(-x+a*log(x)-gln)
      return
      END
      
      SUBROUTINE gcf(gammcf,a,x,gln)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER ITMAX
      REAL*8 a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
      INTEGER i
      REAL*8 an,b,c,d,del,h,gammln
      gln=gammln(a)   
      b=x+1.-a   
      c=1.d0/FPMIN
      d=1.d0/b
      h=d
      do 11 i=1,ITMAX
         an=-i*(i-a)  
         b=b+2
         d=an*d+b
         if(abs(d).lt.FPMIN)d=FPMIN
         c=b+an/c
         if(abs(c).lt.FPMIN)c=FPMIN
         d=1./d
         del=d*c
         h=h*del
         if(abs(del-1.).lt.EPS)goto 1
 11   continue
      pause 'a too large, ITMAX too small in gcf'
 1    gammcf=exp(-x+a*log(x)-gln)*h
      return
      end
      
      FUNCTION gammln(xx)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp   
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     * 24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     * -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
         y=y+1.d0 
         ser=ser+cof(j)/y   
 11      continue     
         gammln=tmp+log(stp*ser/x)
         return
         END
         
      
      
      subroutine field(dipdip,chgchg,chgdip,i,j,iia,jja,ia,
     $     xdif,ydif,zdif,rr,gcfac,dfacx,dfacy,dfacz,
     $     ecc,ecd,edd)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

          include 'parameter.i'      
      
      
c     For periodic boundary conditions, all interactions all counted
      if(ipdc.eq.1) ia=1
      
         di1=dipole(i,iia,1)+dipole_p(i,iia,1)
         di2=dipole(i,iia,2)+dipole_p(i,iia,2)
         di3=dipole(i,iia,3)+dipole_p(i,iia,3)

         dj1=dipole(j,jja,1)+dipole_p(j,jja,1)
         dj2=dipole(j,jja,2)+dipole_p(j,jja,2)
         dj3=dipole(j,jja,3)+dipole_p(j,jja,3)
      
         dir=di1*xdif+di2*ydif+di3*zdif
         djr=dj1*xdif+dj2*ydif+dj3*zdif
      
         
         didj=di1*dj1+di2*dj2+di3*dj3
         dirdjr=dir*djr


      cidjr=qat(iia,i)*djr
      cjdir=qat(jja,j)*dir
      cicj=qat(iia,i)*qat(jja,j)
     
      ecc=chgchg(0)*cicj*fac(4)
      ecd=chgdip(1)*(cidjr-cjdir)*fac(4)
      edd=(dipdip(1)*didj-dipdip(2)*dirdjr)*fac(4)
      gcfac=(chgchg(1)*cicj+dipdip(2)*didj+chgdip(2)*(cidjr-cjdir)
     $     -dipdip(3)*dirdjr)*fac(4)
      if(iloop.eq.1) then
c     for intermolecular interactions- or for POLIR's model- evaluate dipole-charge field
     
         fdcfacj=chgdip(1)*qat(iia,i)*fac(4)
         fdcfaci=chgdip(1)*qat(jja,j)*fac(4)
         EFDC(j,jja,1)=EFDC(j,jja,1)-(FDCFACJ)*XDIF
         EFDC(j,jja,2)=EFDC(j,jja,2)-(FDCFACJ)*YDIF
         EFDC(j,jja,3)=EFDC(j,jja,3)-(FDCFACJ)*ZDIF
         EFDC(i,iia,1)=EFDC(i,iia,1)+(FDCFACI)*XDIF
         EFDC(i,iia,2)=EFDC(i,iia,2)+(FDCFACI)*YDIF
         EFDC(i,iia,3)=EFDC(i,iia,3)+(FDCFACI)*ZDIF
      endif

c     evaluate the potential of the charge from the other
c     charges and dipoles
      
      pot(i,iia)=pot(i,iia)+(qat(jja,j)*chgchg(0)+djr*chgdip(1))*fac(4)
      pot(j,jja)=pot(j,jja)+(qat(iia,i)*chgchg(0)-dir*chgdip(1))*fac(4)
         
      potcc(i,iia)=potcc(i,iia)+qat(jja,j)*chgchg(0)*fac(4)
      potcc(j,jja)=potcc(j,jja)+qat(iia,i)*chgchg(0)*fac(4)
         
      potcd(i,iia)=potcd(i,iia)+djr*chgdip(1)*fac(4)
      potcd(j,jja)=potcd(j,jja)-dir*chgdip(1)*fac(4)
      
      dfacx=(
     $     (qat(jja,j)*di1-qat(iia,i)*dj1)*chgdip(1)
     $     +(djr*di1+dir*dj1)*dipdip(2))*fac(4)
      
      dfacy=(   
     $     (qat(jja,j)*di2-qat(iia,i)*dj2)*chgdip(1)
     $     +(djr*di2+dir*dj2)*dipdip(2))*fac(4)
      
      dfacz=(
     $     (qat(jja,j)*di3-qat(iia,i)*dj3)*chgdip(1)
     $     +(djr*di3+dir*dj3)*dipdip(2))*fac(4)
      
         
         rirj=xdif*xdif
         diptens1(i,j,iia,jja)=fac(4)*(rirj*dipdip(2)-dipdip(1))
         rirj=xdif*ydif
         diptens2(i,j,iia,jja)=fac(4)*(rirj*dipdip(2))
         rirj=xdif*zdif
         diptens3(i,j,iia,jja)=fac(4)*(rirj*dipdip(2))
         rirj=ydif*ydif
         diptens4(i,j,iia,jja)=fac(4)*(rirj*dipdip(2)-dipdip(1))
         rirj=ydif*zdif
         diptens5(i,j,iia,jja)=fac(4)*(rirj*dipdip(2))
         rirj=zdif*zdif
         diptens6(i,j,iia,jja)=fac(4)*(rirj*dipdip(2)-dipdip(1))
      end

      subroutine coulomb(ucc,ucd,udd,uspring)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

             include 'parameter.i'
c             include 'parameter2.i'
c         if(iflag.eq.1)then
c          print*,dpdr(1,3,3,1,1)
c          stop
c          endif
      ucc=0.d0
      ucd=0.d0
      udc=0.d0
      udd=0.d0
      uspring=0.d0
            
      isp=0
      do i=1,nmol
          if(i.eq.1)i_a=1
          if(i.gt.1)i_a=3

            do ia=1,i_a
            isp=isp+1
               
               
               udc=udc-0.5d0*((dipole(i,ia,1)+dipole_p(i,ia,1))
     :   *efdc(i,ia,1)
     $              +(dipole(i,ia,2)+dipole_p(i,ia,2))*efdc(i,ia,2)
     $              +(dipole(i,ia,3)+dipole_p(i,ia,3))*efdc(i,ia,3))
               udd=udd-0.5d0*((dipole(i,ia,1)+dipole_p(i,ia,1))
     :    *efdd(i,ia,1)
     $              +(dipole(i,ia,2)+dipole_p(i,ia,2))*efdd(i,ia,2)
     $        +(dipole(i,ia,3)+dipole_p(i,ia,3))*efdd(i,ia,3))
               d2=dipole(i,ia,1)*dipole(i,ia,1)
     $              +dipole(i,ia,2)*dipole(i,ia,2)
     $              +dipole(i,ia,3)*dipole(i,ia,3)
            
               uspring=uspring+fac(4)*d2/(2.d0*alph(i,ia))

            ucc=ucc+0.5d0*qat(ia,i)*potcc(i,ia) 
            ucd=ucd+0.5d0*qat(ia,i)*potcd(i,ia)
           enddo
           enddo

c     ucd and udc should be exactly the same (which can be checked).
c     The first is the energy of the charges in the field of the dipoles
c     and the second is the energy of the dipoles in the field of the
c     charges.  They sum to give the total charge-dipole energy.

      ucd=udc+ucd
      do i=1,nmol
          if(i.eq.1)i_a=1
          if(i.gt.1)i_a=3

         do ia=1,i_a
              if(ia.eq.1)iia=3
              if(ia.eq.2)iia=1
              if(ia.eq.3)iia=2

                  do idim=1,3 
                  do ja=1,3                
              if(ja.eq.1)jja=3
              if(ja.eq.2)jja=1
              if(ja.eq.3)jja=2

                  do idim2=1,3
                     force(i,ia,idim)=force(i,ia,idim)
     $                    +(efdc(i,ja,idim2)+efdd(i,ja,idim2))
     $                    *(dpdr(i,jja,iia,idim2,idim))
c            print*,(efdc(i,ja,idim2)+efdd(i,ja,idim2))
c     $                    *(dpdr(i,jja,iia,idim2,idim))
              !print*,dpdr(i,jja,iia,idim2,idim)
              ! print*,i,jja,iia,idim2,idim
            !stop
                  enddo
               enddo
            enddo
         enddo
      enddo
                  
      end

      function erfcc(x)
      real*8 erfcc,x
      real*8 t,z
      z=dabs(x)
      t=1.d0/(1.+0.5d0*z)  
      erfcc=t*exp(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+t*
     $ (.09678418d0+t*(-.18628806d0+t*(.27886807d0+t*(-1.13520398d0+t*
     $ (1.48851587d0+t*(-0.82215223d0+t*.17087277d0)))))))))
      if(x.lt.0.d0) erfcc=2.d0-erfcc
      return
      end
      
         
      function erf(x)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 erf,x
      real*8 gammq
      if(x.lt.0.) then
         erf=1.+gammq(.5,x**2)
      else
         erf=1.-gammq(.5,x**2)  
      endif
      return
      end
         

