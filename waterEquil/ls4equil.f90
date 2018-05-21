      program main
        implicit real*8(a-h,o-z)
        character*150 inputfile           
        character*1 ch1
        include 'parameter.i' 
        real*8 fmass1(3),vel(nmol,3,3)   

!         print*, '** INPUT FILES= '
        read*, xyzFILE,cvecFILE
!         print*, '** OUTPUT FILES= '
        read*, xyZOUT,velOUT,tempOUT,peOUT,keOUT,totOUT,md_trajOUT,md_velOUT

        open(unit=2,file=xyzFILE)
        open(unit=3,file=cvecFILE)

        open(unit=22,file=tempOUT)
        open(unit=23,file=peOUT)
        open(unit=24,file=keOUT)
        open(unit=25,file=totOUT)
        open(unit=300,file=md_trajOUT)
        open(unit=400,file=md_velOUT)

!         print*, '** INPUT PARAMETERS= '
        read*, nstep,firsttimestep,alph(1,1),c_16,c_14,c_12,c_6,fmassmg
        nstep=nstep+firsttimestep      
        fmassmg=fmassmg*1.6605402D-27      

        read(3,*)cvec(1,1),cvec(1,2),cvec(1,3)
        read(3,*)cvec(2,1),cvec(2,2),cvec(2,3)
        read(3,*)cvec(3,1),cvec(3,2),cvec(3,3)
        box=cvec(1,1)
!         nmol=natoms/3
        !CHANGE
        !remove two "READ(2)"commands here

        do i=1,nmol
        !CHANGE
        !fix loop conditions
        !if(i.eq.1)i_a=1
        !if(i.gt.1)i_a=3
        i_a = 3
        do j=1,i_a 
        read(2,*)ch1,rat(i,j,1),rat(i,j,2),rat(i,j,3)
        enddo
        if(i.eq.1)then
!                     rat(i,1,1)=rat(i,2,1)
!                     rat(i,1,2)=rat(i,2,2)
!                     rat(i,1,3)=rat(i,2,3)
        endif
        enddo
        close(2)

        iton=1
!         rat(10,1,2)=rat(10,1,2)+0.0005d0
        pi=acos(-1.d0)
        c16=7.61978e-09*1d+16
        c14=-4.07504e-07*1d+14
        c12=5.80846e-06*1d+12
        c6=-0.00292171*1d+6
        epso=8.854187817D-12
        esconvert=1.d0/(4.0d0*pi*epso)
        temp=300.d0
        temp_t=300.d0
        dt=0.0002d-12
        fkb= 1.3806503d-23
        fmasso=15.9994d0*1.6605402D-27      
        fmassh=1.d0*1.6605402D-27      
        fmass1(1)=fmasso
        fmass1(2)=fmassh
        fmass1(3)=fmassh
        vo=dsqrt(fkb*temp/fmasso)
        vh=dsqrt(fkb*temp/fmassh)

        xmax=10.d0
        ake=0.d0

        do i=1,nmol
        if(i.eq.1)i_a=1
        if(i.gt.1)i_a=3
        do ia=1,i_a
        do idim=1,3
 34               xtry=(2.d0*ran0(idum)-1.d0)*xmax
        ytry=ran0(idum)
        if(ytry.ge.dexp(-0.5d0*xtry*xtry))goto 34
        if(i.eq.1)then
          vel(i,ia,idim)=xtry*sqrt(fkb*temp_t/fmassmg)
          ake=ake+0.5d0*fmassmg*vel(i,ia,idim)**2
        else
          vel(i,ia,idim)=xtry*sqrt(fkb*temp_t/fmass1(ia))
          ake=ake+0.5d0*fmass1(ia)*vel(i,ia,idim)**2
        endif
        enddo
        enddo
        enddo
        sum_x=0.d0
        sum_z=0.d0
        sum_y=0.d0
        temp1=ake/(0.5d0*float((nmol-1)*9+3)*fkb)

        temp_t=300.d0
        do i=1,nmol
        if(i.eq.1)i_a=1
        if(i.gt.1)i_a=3
        do ia=1,i_a
        if(i.eq.1)then
          sum_x=sum_x+fmassmg*vel(i,ia,1)
          sum_y=sum_y+fmassmg*vel(i,ia,2)
          sum_z=sum_z+fmassmg*vel(i,ia,3)
        else
          sum_x=sum_x+fmass1(ia)*vel(i,ia,1)
          sum_y=sum_y+fmass1(ia)*vel(i,ia,2)
          sum_z=sum_z+fmass1(ia)*vel(i,ia,3)
        endif
        enddo
        enddo
        print*,'1',sum_x,sum_y,sum_z
        sum_x=sum_x/(float((nmol-1)*3+1))
        sum_y=sum_y/(float((nmol-1)*3+1))
        sum_z=sum_z/(float((nmol-1)*3+1))
        ake=0.d0

        do i=1,nmol
        if(i.eq.1)i_a=1
        if(i.gt.1)i_a=3
        do ia=1,i_a
        if(i.eq.1)then
          vel(i,ia,1)=vel(i,ia,1)-sum_x/fmassmg
          vel(i,ia,2)=vel(i,ia,2)-sum_y/fmassmg
          vel(i,ia,3)=vel(i,ia,3)-sum_z/fmassmg
          ake=ake+0.5d0*fmassmg*vel(i,ia,1)**2
          ake=ake+0.5d0*fmassmg*vel(i,ia,2)**2
          ake=ake+0.5d0*fmassmg*vel(i,ia,3)**2
        else
          vel(i,ia,1)=vel(i,ia,1)-sum_x/fmass1(ia)
          vel(i,ia,2)=vel(i,ia,2)-sum_y/fmass1(ia)
          vel(i,ia,3)=vel(i,ia,3)-sum_z/fmass1(ia)
          ake=ake+0.5d0*fmass1(ia)*vel(i,ia,1)**2
          ake=ake+0.5d0*fmass1(ia)*vel(i,ia,2)**2
          ake=ake+0.5d0*fmass1(ia)*vel(i,ia,3)**2
        endif

        enddo
        enddo

        do i=1,nmol
        if(i.eq.1)i_a=1
        if(i.gt.1)i_a=3
        do ia=1,i_a
        if(i.eq.1)then
          sum_x=sum_x+fmassmg*vel(i,ia,1)
          sum_y=sum_y+fmassmg*vel(i,ia,2)
          sum_z=sum_z+fmassmg*vel(i,ia,3)
        else
          sum_x=sum_x+fmass1(ia)*vel(i,ia,1)
          sum_y=sum_y+fmass1(ia)*vel(i,ia,2)
          sum_z=sum_z+fmass1(ia)*vel(i,ia,3)
        endif
        enddo
        enddo
        print*,'2',sum_x,sum_y,sum_z

        temp_t=300.d0
        temp1=ake/(0.5d0*float((nmol-1)*9+3)*fkb)
        call scale_temp(temp1,temp_t,vel)

        sum_x=0.d0
        sum_y=0.d0
        sum_z=0.d0

        do i=1,nmol
        if(i.eq.1)i_a=1
        if(i.gt.1)i_a=3
        do ia=1,i_a
        if(i.eq.1)then
          sum_x=sum_x+fmassmg*vel(i,ia,1)
          sum_y=sum_y+fmassmg*vel(i,ia,2)
          sum_z=sum_z+fmassmg*vel(i,ia,3)
        else
          sum_x=sum_x+fmass1(ia)*vel(i,ia,1)
          sum_y=sum_y+fmass1(ia)*vel(i,ia,2)
          sum_z=sum_z+fmass1(ia)*vel(i,ia,3)
        endif
        enddo
        enddo
        print*,'3',sum_x,sum_y,sum_z

        call  recip(vol,voli,rcut)
        avsno=6.0221367D+23

        do iblk=1,1!1
        qat(1,1)=-1.d0!float(iblk-1)*0.2d0
        call subtract(utot,upair,ucc,ucd,udd,uspring,uintra,ucx,udx,iton)
        print*,utot,force(10,1,2)
!            stop
        do n=firsttimestep,nstep
        print*,n!,pott
!               if(n.eq.8)goto 98
        ake=0.d0 
        do i=1,nmol
        if(i.eq.1)i_a=1
        if(i.gt.1)i_a=3

        do ia=1,i_a
        do idim=1,3
        if(i.eq.1) then 
          fmass(ia)=fmassmg
        else
          fmass(ia)=fmass1(ia)
        endif
        acc(i,ia,idim)=force(i,ia,idim)/(1.439325215D+20*fmass(ia))
        acc(i,ia,idim)=acc(i,ia,idim)/1d-10
        rat(i,ia,idim)=rat(i,ia,idim)+dt*vel(i,ia,idim)*1d+10 &
          +0.5d0*dt*dt*acc(i,ia,idim)*1d+10
        enddo
        enddo
        enddo
        call subtract(utot,upair,ucc,ucd,udd,uspring,uintra,ucx,udx,iton)

        do i=1,nmol
        if(i.eq.1)i_a=1
        if(i.gt.1)i_a=3
        do ia=1,i_a
        vv=0.d0
        if(i.eq.1) then 
          fmass(ia)=fmassmg
        else
          fmass(ia)=fmass1(ia)
        endif

        do idim=1,3
        acc1(i,ia,idim)=force(i,ia,idim)/(1.439325215D+20*fmass(ia))
        acc1(i,ia,idim)=acc1(i,ia,idim)/1d-10
        vold=vel(i,ia,idim)
        vel(i,ia,idim)=vold+0.5d0*dt*(acc1(i,ia,idim)+acc(i,ia,idim))
        vv=vv+vel(i,ia,idim)**2
        enddo
        ake=ake+0.5d0*fmass(ia)*vv
        enddo
        enddo

        temp1=ake/(0.5d0*float((nmol-1)*9+3)*fkb)
        !print*,utot+ake*avsno/4184.d0,temp1,utot+0.5d0*temp1*1.985877534d-3*float((nmol-1)*9+3)
        print *, "***************************************************"
        print *, "UTOT: ", utot
        print *, "UCC: ", ucc
        print *, "UCD: ", ucd
        print *, "USPRING: ", uspring
        print *, "UPAIR: ", upair
        call scale_temp(temp1,temp_t,vel)
        print *, "***************************************************"
!               if ((n.ge.1).and.(n.le.1000)) call stat_collect(n,rat,box,nstep)
        write(22,*)n,temp1!utot+ake*avsno/4184.d0
        write(23,*)n,utot!+ake*avsno/4184.d0
        write(24,*)n,ake*avsno/4184.d0
        write(25,*)n,utot+ake*avsno/4184.d0

        !CHANGE
        !remove step check to output traj at every step
        !CHANGE
        !add atom count to md output
        !and space between atom count and output
        write(300,'(a,3f10.4)')'6'
        write(300,*)
        do i=1,nmol
        !if(i.eq.1)i_a=1
        !if(i.gt.1)i_a=3
        i_a=3
        !CHANGE i_a should always be three here

        do ia=1,i_a
        write(400,*)vel(i,ia,1),vel(i,ia,2),vel(i,ia,3)
        if(ia.eq.1)then
          write(300,'(a,3f10.4)')'O',rat(i,ia,1),rat(i,ia,2),rat(i,ia,3)
        else
          write(300,'(a,3f10.4)')'H',rat(i,ia,1),rat(i,ia,2),rat(i,ia,3)
        endif
        enddo
        enddo
        enddo

 98         open(unit=100,file=xyZOUT)
        open(unit=200,file=velOUT)
        write(100,*)
        write(100,*)
        do i=1,nmol
        if(i.eq.1)i_a=1
        if(i.gt.1)i_a=3
        do ia=1,i_a
        write(200,*)vel(i,ia,1),vel(i,ia,2),vel(i,ia,3)
        if(ia.eq.1)then
          write(100,'(a,3f10.4)')'O',rat(i,ia,1),rat(i,ia,2),rat(i,ia,3)
        else
          write(100,'(a,3f10.4)')'H',rat(i,ia,1),rat(i,ia,2),rat(i,ia,3)
        endif
        enddo
        enddo
        close(200)
        close(100)
        enddo
        end         


        subroutine stat_collect(istep,rat,box,nstep)
          implicit real*8(a-h,o-z)
          parameter(nmol=128)
          parameter(dgr=0.05d0)
          dimension rat(nmol,3,3),grhistm(5000),grhist(5000)
          save grhist,grhistm
          character*30, grooOUT,gromOUT
          grooOUT="grWater.out"
          gromOUT="grIon.out"

          do i=2, nmol-1
          do j=i+1, nmol
          x=rat(i, 1, 1)-rat(j, 1, 1)
          y=rat(i, 1, 2)-rat(j, 1, 2)
          z=rat(i, 1, 3)-rat(j, 1, 3)
          x=x-idnint(x/box)*box
          y=y-idnint(y/box)*box
          z=z-idnint(z/box)*box
          rij=sqrt(x**2+y**2+z**2)
          ibin=int(rij/dgr)+1
          grhist(ibin)=grhist(ibin)+2.d0
          enddo
          enddo

          do i=1,1
          do j=i+1, nmol
          x=rat(i, 1, 1)-rat(j, 1, 1)
          y=rat(i, 1, 2)-rat(j, 1, 2)
          z=rat(i, 1, 3)-rat(j, 1, 3)
          x=x-idnint(x/box)*box
          y=y-idnint(y/box)*box
          z=z-idnint(z/box)*box
          rij=sqrt(x**2+y**2+z**2)
          ibin=int(rij/dgr)+1
          grhistm(ibin)=grhistm(ibin)+1.d0
          enddo
          enddo

          if(mod(istep,1000).eq.0)then
            open(unit=20,file=grooOUT)
            open(unit=21,file=gromOUT)
            pi=acos(-1.d0)
            dens=float(nmol)/box**3
            c=4.d0*pi*(dens)/3.d0
            ibin_max=int(0.5d0*box/dgr)
            do ibin=1, ibin_max
            rlower=float(ibin-1)*dgr
            rupper=rlower+dgr
            xideal=c*(rupper**3-rlower**3)
            xliq=grhist(ibin)/float(istep)
            gr=xliq/xideal/float(nmol)
            write (20, *) rlower+0.5d0*dgr, gr

            xliq=grhistm(ibin)/float(istep)
            gr=xliq/xideal/float(nmol)
            write (21, *) rlower+0.5d0*dgr, gr
            enddo
            close (20) 
            close (21) 
          endif
          return
          end

          subroutine scale_temp(temp1,temp_t,vel)
            implicit real*8(a-h,o-z)
            parameter(nmol=128)
            real*8 vel(nmol,3,3)
            scale=dsqrt(temp_t/temp1)
            do i=1,nmol
            if(i.eq.1)i_a=1
            if(i.gt.1)i_a=3
            do ia=1,i_a
            do idim=1,3
            vel(i,ia,idim)=vel(i,ia,idim)*scale
            enddo
            enddo       
            enddo
            end

            FUNCTION ran0(idum)
              IMPLICIT DOUBLE PRECISION (A-H,O-Z)
              INTEGER idum,IA,IM,IQ,IR,MASK
              REAL*8 ran0,AM
              PARAMETER(IA=16807,IM=2147483647,AM=1.d0/IM,IQ=127773,IR=2836,MASK=123459876)
              INTEGER k
              idum=ieor(idum,MASK)
              k=idum/IQ
              idum=IA*(idum-k*IQ)-IR*k
              if(idum.lt.0d0) idum=idum+IM
              ran0=AM*idum
              idum=ieor(idum,MASK)
              return
              end

