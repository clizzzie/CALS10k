      program fielduncert3k
      implicit none
c     Modified version of fieldpred.f to give time series of declination
c     inclination and intensity with uncertainty estimates 
c     from models ARCH3k.1MAST, SED3k.1MAST or CALS3k.3MAST
c     CALS3k.3b, CALS3k.4b
c     The desired model file has to be in the same directory as the program.
c     Monika Korte, November 2010
c
c     uses pieces of code by Jeremy Bloxham, Andrew Jackson, David Gubbins
c     Kathy Whaler, David Barraclough, Rick 0'Connel and Carl de Boor
c---------------------------------------------------------------------
c     for details of B-splines see: 
c            Carl de Boor "A Practical Guide to Splines"
c            Springer-Verlag, 1978.
c--------------------------------------------------------------------

      integer lmax,nspl,n,np,nl,jord,nsplt
      real*8 gt,spl,tknts,g,gd,p,dp,dx,dy,dz,fac
      real*8 dg,dgt,ex,ey,ez,eh,ef,ed,ei
      character*30 outfile,modfile

      parameter (lmax=10)
      parameter (nsplt=402)

      parameter (n=lmax*(lmax+2))
      parameter (np=n*nsplt)
      parameter (nl=(lmax+1)*(lmax+2)/2)

      dimension gt(n,nsplt),dgt(n,nsplt)
      dimension spl(nsplt),tknts(nsplt+4)
      
      dimension g(n),gd(n),dg(n)
      
      dimension p(nl),dp(nl)
      dimension dx(lmax*(lmax+2)),dy(lmax*(lmax+2)),dz(lmax*(lmax+2))

      integer it1,it2,lm,nm,k,j,i,nleft,flag
      real*8 alat,alon,alt,time,theta,phi,rad,sinth,costh,sd,cd
      real*8 x,y,z,h,f,ainc,d,tstartin,tendin

      data jord/4/
      data fac/1.74532925e-2/
      
      write(*,*) '            -Program Fielduncert3k- '
      write(*,*) 'Produces time series of magnetic field components'
      write(*,*) 'declination, inclination and intensity'
      write(*,*) 'for a location of latitude and longitude'
      write(*,*) 'with uncertainty estimates as given by'
      write(*,*) 'MAST estimation procedure for'
      write(*,*) 'the model coefficients.'
      write(*,*) 'Results are written to a plain text output file.'
      write(*,*)
      write(*,*) 'Choose model: 1 - CALS3k.3MAST'
      write(*,*) '              2 - ARCH3k.1MAST'
      write(*,*) '              3 - SED3k.1MAST'
      write(*,*) '              4 - CALS3k.3b'
      write(*,*) '              5 - CALS3k.4b'
      read(*,*) flag
      write(*,*) 'Give output file name:'
      read(*,*) outfile
      write(*,*) 'Give latitude in decimal degrees (S negative):'
      read(*,*) alat
      write(*,*) 'Give longitude in decimal degrees (W negative):'
      read(*,*) alon

      if (flag.eq.1) then
         modfile='CALS3k.3MAST'
      else if (flag.eq.2) then
         modfile='ARCH3k.1MAST'
      else if (flag.eq.3) then
         modfile='SED3k.1MAST'
      else if (flag.eq.4) then
         modfile='CALS3k.3b'
      else if (flag.eq.5) then
         modfile='CALS3k.4b'
      else 
         write(*,*) 'ERROR: invalid model choice'
         stop        
      end if

c********************************************************************
c     read model, block of uncertainties at the end
      open(7,file=modfile)

      read(7,*) tstartin,tendin
      read(7,*) lm,nm,nspl,(tknts(i),i=1,nspl+4)
      read(7,*) gt
      read(7,*) dgt
      close(7)

      it1=-1000
      it2=1990

      open(11,file=outfile)
      write(11,*) 'Model ',modfile,'Location lat=',alat,' long=',alon
      write(11,*) 'Year D(deg) I(deg) F(microT) deltaD deltaI deltaF'
      do i=it1,it2,1
      time = float(i)
      alt=0.0
      theta = (90.0-alat)*fac
      phi   = alon*fac
      call coords(alt,theta,rad,sd,cd)      

      sinth=sin(theta)
      costh=cos(theta)
      call plmbar(p,dp,costh,lmax)

c-----
c     calculate main field coefficients and uncertainties at time time
10    call interv(tknts,time,nspl,nleft)
      call bspline(tknts,time,nspl,jord,nleft,spl(nleft-3))
      
      do  k=1,n
       g(k)=0.0
       dg(k)=0.0
       do j=1,4
        g(k) = g(k) + spl(j+nleft-4)*gt(k,j+nleft-4)
        dg(k)=dg(k) + spl(j+nleft-4)*dgt(k,j+nleft-4)
       enddo 
      enddo 

      call magfdz(p,dp,theta,phi,rad,lmax,g,dx,dy,dz,x,y,z,h,f,
     >ainc,d,sd,cd)

      call errfdz(p,dp,theta,phi,rad,lmax,dg,dx,dy,dz,ex,ey,ez,sd,cd)

      eh=sqrt((1/h**2)*((x*ex)**2+(y*ey)**2))
      ef=sqrt((1/f**2)*((x*ex)**2+(y*ey)**2+(z*ez)**2))
      ei=sqrt((1/(1+(z/h)**2))**2*((ez/h)**2+((z*eh)/h**2)**2))
      ed=sqrt((1/(1+(y/x)**2))**2*((ey/x)**2+((y*ex)/x**2)**2))

      ainc=ainc/fac
      ei=ei/fac
      d=d/fac
      ed=ed/fac
      f=f/1000.
      ef=ef/1000.

      write(11,6200) i,d,ainc,f,ed,ei,ef
      end do
99    continue
 6100 format(i6,6f10.1)
6200  format(i6,2f8.2,f6.1,3f8.2)
      close(11)

      stop
      end      



      subroutine errfdz(p,dp,theta,phi,r,lmax,g,dx,dy,dz,
     >x,y,z,sd,cd)
c     calculate uncertainties for field predictions from uncertainties
c     in coefficients g
c     
c     use error propagation rules

      implicit real*8 (a-h,o-z)
      dimension g(lmax*(lmax+2))
      dimension dx(lmax*(lmax+2)),dy(lmax*(lmax+2)),dz(lmax*(lmax+2))
      dimension p((lmax+1)*(lmax+2)/2),dp((lmax+1)*(lmax+2)/2)
      real*8 i


      b=6371.2/r
      x=0.
      y=0.
      z=0.
      sinth=sin(theta)
      if(abs(sinth).lt.1.e-10) sinth=1.e-10

      do 20 l=1,lmax

      l1=l+1
      bb=b**(l+2)
      k=l*l
      k1=(l*l1)/2+1

      dx(k)=dp(k1)*bb
      dy(k)=0.
      dz(k)=-p(k1)*l1*bb
      x=x+(g(k)*dx(k))**2
      z=z+(g(k)*dz(k))**2

      do 20 m=1,l

      t=float(m)*phi
      k=l*l+2*m-1
      k1=(l*l1)/2+m+1
      sint=sin(t)
      cost=cos(t)

      dxd = dp(k1)*bb
      dx(k) = dxd*cost
      dx(k+1) = dxd*sint
      x = x + (g(k)*dx(k))**2 + (g(k+1)*dx(k+1))**2

      dxy = m*p(k1)*bb/sinth
      dy(k) = dxy*sint
      dy(k+1) = -dxy*cost
      y = y + (g(k)*dy(k))**2 + (g(k+1)*dy(k+1))**2

      dzd = -l1*p(k1)*bb
      dz(k) = dzd*cost
      dz(k+1) = dzd*sint
      z = z + (g(k)*dz(k))**2 + (g(k+1)*dz(k+1))**2

20    continue
      
      xs = x
      x = x*(cd**2) + z*(sd**2)
      z = z*(cd**2) - xs*(sd**2)
   
      x=sqrt(x)
      y=sqrt(y)
      z=sqrt(z)

      do 50 k=1,lmax*(lmax+2)
      dxk = dx(k)
      dzk = dz(k)
      dx(k) = dxk*cd + dzk*sd
      dz(k) = dzk*cd - dxk*sd
50    continue
      
      return
      end

c--------------------------------------------------------------------------      
      
      subroutine interv(tknts,time,nspl,nleft)
      implicit real*8 (a-h,o-z)
      
      dimension tknts(nspl+4)
     
      if(time.lt.tknts(4).or.time.gt.tknts(nspl+1)) return
      
      do 200 n=5,nspl+1
       if(time.le.tknts(n)) then
        nleft=n-1
        goto 210
       endif
200   continue
210   continue


      return
      end

c-------------------------------------------------------------------

       subroutine bspline(tknts,t,nspl,jorder,nleft,spl)
 
c calculate splines of order jorder where 1 <= jorder <= 4
       implicit real*8 (a-h,o-z)
       dimension tknts(nspl+4)
       dimension spl(4)
       
       dimension deltal(4),deltar(4)
       
       spl(1)=1.0
      
       do 200 j=1,jorder-1
       
       deltar(j) = tknts(nleft+j) - t
       deltal(j) = t - tknts(nleft+1-j)
       saved=0.0
       
       do 100 i=1,j
        term = spl(i)/(deltar(i)+deltal(j+1-i))
        spl(i) = saved + deltar(i)*term
        saved = deltal(j+1-i)*term
100    continue

       spl(j+1) = saved
       
200    continue

       
 
       return
       end
        
c-----------------------------------------------------------------

      subroutine coords(h,theta,r,sd,cd)
      implicit real*8 (a-h,o-z)

      pi=3.14159265
      b1=40680925.
      b2=40408585.
      theta=pi/2-theta
      clat=cos(theta)
      slat=sin(theta)
      one=b1*clat*clat
      two=b2*slat*slat
      three=one+two
      four=sqrt(three)
      r=sqrt(h*(h+2.*four)+(b1*one+b2*two)/three)
      cd=(h+four)/r
      sd=(b1-b2)/four*slat*clat/r
      sinth=slat*cd-clat*sd
      costh=clat*cd+slat*sd
      theta=pi/2.-atan2(sinth,costh)
      return
      end

c-----------------------------------------------------------------

      subroutine magfdz(p,dp,theta,phi,r,lmax,g,dx,dy,dz,
     >x,y,z,h,f,i,d,
     >sd,cd)

c
c***************************************************************
c
c     j bloxham  8 nov 1982 & 11 oct 1983
c
c     modified version of dg13.sv.fort:magfd & jb62.sv.progs:magfds
c
c     gives field components at radius r
c
c***************************************************************
c
c
c     this version 16 jan 87
c
c     saves dx dy dz in computation
c
cc======================================================================

      implicit real*8 (a-h,o-z)
      dimension g(lmax*(lmax+2))
      dimension dx(lmax*(lmax+2)),dy(lmax*(lmax+2)),dz(lmax*(lmax+2))
      dimension p((lmax+1)*(lmax+2)/2),dp((lmax+1)*(lmax+2)/2)
      real*8 i


      b=6371.2/r
      x=0.
      y=0.
      z=0.
      sinth=sin(theta)
      if(abs(sinth).lt.1.e-10) sinth=1.e-10

      do 20 l=1,lmax

      l1=l+1
      bb=b**(l+2)
      k=l*l
      k1=(l*l1)/2+1

      dx(k)=dp(k1)*bb
      dy(k)=0.
      dz(k)=-p(k1)*l1*bb
      x=x+g(k)*dx(k)
      z=z+g(k)*dz(k)

      do 20 m=1,l

      t=float(m)*phi
      k=l*l+2*m-1
      k1=(l*l1)/2+m+1
      sint=sin(t)
      cost=cos(t)

      dxd = dp(k1)*bb
      dx(k) = dxd*cost
      dx(k+1) = dxd*sint
      x = x + (g(k)*dx(k)) + (g(k+1)*dx(k+1))

      dxy = m*p(k1)*bb/sinth
      dy(k) = dxy*sint
      dy(k+1) = -dxy*cost
      y = y + (g(k)*dy(k)) + (g(k+1)*dy(k+1))

      dzd = -l1*p(k1)*bb
      dz(k) = dzd*cost
      dz(k+1) = dzd*sint
      z = z + (g(k)*dz(k)) + (g(k+1)*dz(k+1))

20    continue
      
      xs = x
      x = x*cd + z*sd
      z = z*cd - xs*sd
   
      do 50 k=1,lmax*(lmax+2)
      dxk = dx(k)
      dzk = dz(k)
      dx(k) = dxk*cd + dzk*sd
      dz(k) = dzk*cd - dxk*sd
50    continue
    
      
      h=sqrt(x*x+y*y)
      f=sqrt(h*h+z*z)
      i=asin(z/f)
      d=atan2(y,x)

      return
      end
       
c---------------------------------------------------------------------

      subroutine plmbar(p,dp,z,lmax)
c
c  evaluates normalized associated legendre function p(l,m) as function of
c   z=cos(colatitude) using recurrence relation starting with p(l,l) 
c   and then increasing l keeping m fixed.  normalization is: 
c   integral(y(l,m)*y(l,m))=4.*pi, where y(l,m) = p(l,m)*exp(i*m*longitude),
c   which is incorporated into the recurrence relation. p(k) contains p(l,m)
c   with k=(l+1)*l/2+m+1; i.e. m increments through range 0 to l before 
c   incrementing l. routine is stable in single and double precision to
c   l,m = 511 at least; timing proportional to lmax**2
c   r.j.o'connell 7 sept. 1989

c   a.jackson 19 october 1989  code added at end:
c   (2) derivatives added and stored in dp(k)
c       using same arrangement as for p(k)
c
      implicit real*8(a-h,o-z)
      dimension p(*),dp(*)
c     --dimension of p, dp must be (lmax+1)*(lmax+2)/2 in calling program
      if (lmax.lt.0.or.abs(z).gt.1.d0) stop 'bad arguments'
c       --case for p(l,0) 
        pm2=1.d0
        p(1)=1.d0
        dp(1)=0.d0
        if (lmax .eq. 0) return
        pm1=z
        p(2)=dsqrt(3.d0)*pm1
        k=2
        do 4 l=2,lmax
          k=k+l
          plm=(dfloat(2*l-1)*z*pm1-dfloat(l-1)*pm2)/dfloat(l)
          p(k)=dsqrt(dfloat(2*l+1))*plm
          pm2=pm1
4         pm1=plm
c       --case for m > 0
        pmm = 1.d0
        sintsq = (1.d0-z)*(1.d0+z)
        fnum = -1.0d0
        fden = 0.0d0
        kstart = 1
        do 20 m =1 ,lmax
c         --case for p(m,m) 
          kstart = kstart+m+1
          fnum = fnum+2.0d0
          fden = fden+2.0d0
          pmm = pmm*sintsq*fnum/fden
          pm2 = dsqrt(dfloat(4*m+2)*pmm)
          p(kstart) = pm2
          if (m .eq. lmax) goto 100
c         --case for p(m+1,m)
          pm1=z*dsqrt(dfloat(2*m+3))*pm2
          k = kstart+m+1
          p(k) = pm1
c         --case for p(l,m) with l > m+1
          if (m .lt. (lmax-1)) then
           do 10 l = m+2,lmax
            k = k+l
            f1=dsqrt(dfloat((2*l+1)*(2*l-1))/dfloat((l+m)*(l-m)))
            f2=dsqrt(dfloat((2*l+1)*(l-m-1)*(l+m-1))
     &              /dfloat((2*l-3)*(l+m)*(l-m)))
            plm=z*f1*pm1-f2*pm2
            p(k) = plm
            pm2 = pm1
10          pm1 = plm
          endif
20        continue

100     continue

c       Gauss-Schmidt normalisation:
        k=1
        do 30 l=1,lmax
        fac=1.d0/dsqrt(dfloat(2*l+1))
        do 30 m=0,l
        k=k+1
        p(k)=p(k)*fac
30      continue

c       now find derivatives of p(z) wrt theta, where z=cos(theta)
        dp(2)=-p(3)
        dp(3)=p(2)
        k=3
        do 200 l=2,lmax
        
          k=k+1
c         treat m=0 and m=l separately
          dp(k)=-dsqrt(dfloat(l*(l+1))/2.d0)*p(k+1)
          dp(k+l)=dsqrt(dfloat(l)/2.d0)*p(k+l-1)
          do 300 m=1,l-1
            k=k+1
            fac1=dsqrt( dfloat( (l-m)*(l+m+1) ) )
            fac2=dsqrt( dfloat( (l+m)*(l-m+1) ) )
            if(m.eq.1)fac2=fac2*dsqrt(2.d0)
            dp(k)=0.5d0*( fac2*p(k-1) - fac1*p(k+1) )
300       continue
          k=k+1

200     continue
        return
        end

c------------------------------------------------------------------

      subroutine bspline1(tknts,t,nspl,nl,spl)
      
      implicit real*8(a-h,o-z)
      parameter(nsplmx=164)
      dimension tknts(nspl+4)
      dimension spl(nspl),spl1(nsplmx)
     
      data jord/3/
      
      if(nspl.gt.nsplmx)then
      write(6,*)' increase dimensions of spl1 in bspline1'
      stop
      endif
      
      do 100 is=1,nspl
       spl(is)=0.0
100   continue
      
      call bspline(tknts,t,nspl,jord,nl,spl1(nl-2))
      
      spl(nl) = b0(tknts,nl)*spl1(nl)
      
      spl(nl-1) = b0(tknts,nl-1)*spl1(nl-1) + b1(tknts,nl-1)*spl1(nl)
      spl(nl-2) = b0(tknts,nl-2)*spl1(nl-2) + b1(tknts,nl-2)*spl1(nl-1)
      spl(nl-3) = b1(tknts,nl-3)*spl1(nl-2)
      
      return
      end

      function b0(tknts,i)
      implicit real*8(a-h,o-z)
      dimension tknts(*)
      
      b0 = 3.0/( tknts(i+3) - tknts(i) )
 
      return
      end      
      
      function b1(tknts,i)
      implicit real*8(a-h,o-z)
      dimension tknts(*)
      
      b1 = -3.0/( tknts(i+4) - tknts(i+1) )
 
      return
      end
      
