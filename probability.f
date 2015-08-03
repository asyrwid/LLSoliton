      PROGRAM LIEB_SOLITONS  
c     Many-particle probability by Metropolis algorithm. 

      implicit double precision (a-h,o-z)
c      parameter (Natoms=11,nperm=39916800)
      parameter (Natoms=8,nperm=40320)
c      parameter (Natoms=5,nperm=120)
c      parameter (Natoms=3,nperm=6)
      parameter (gam=0.01d0)
      logical nextp,nextpsign,ltemp
      integer permutation(Natoms),permsign
      real*8 k(Natoms)
      real*8 x(Natoms),xt(Natoms),xs(Natoms)
      complex*16 ampl(nperm)  !! size=Factorial(Natoms)
      complex*16 wf,imag,ctemp


      pi=4.d0*datan(1.d0)
      pi2=2.d0*pi
      imag=dcmplx(0.d0,1.d0)

      c=gam*Natoms
      write(*,*) 'c=',c      

c      open(1,file='k1.LL',status='old')
      open(1,file='k0.LL',status='old')
      do j=1,Natoms
        read(1,*) k(j)
      enddo
      close(1)

      do ind=1,Natoms
        permutation(ind)=ind
      enddo
      permsign=1

      do indper=1,nperm  !! amplitudes of many-particle wf
        ctemp=1.d0
        do i=1,Natoms-1
        do j=i+1,Natoms
          ctemp=ctemp*(k(permutation(j))-k(permutation(i))-imag*c)
        enddo
        enddo
        ampl(indper)=permsign*ctemp
        ltemp=nextpsign(Natoms,permutation,permsign)
      enddo

c ---------------------------------------------------------------------
c ------------ initial point for Metropolis algorithm -----------------

      idum=-1
      do j=1,Natoms
        x(j)=ran1(idum)
        xs(j)=x(j)
      enddo

      do j=Natoms,1,-1
        do i=1,j
          if(xs(i-1) .gt. xs(i)) then
            temp=xs(i)
            xs(i)=xs(i-1)
            xs(i-1)=temp
          endif
        enddo
      enddo  

      do ind=1,Natoms
        permutation(ind)=ind
      enddo
      wf=0.d0
      do indper=1,nperm  !! many-particle wave function
        temp=0.d0
        do j=1,Natoms
          temp=temp+k(permutation(j))*xs(j)
        enddo
        wf=wf+ampl(indper)*cdexp(imag*temp)
        ltemp=nextp(Natoms,permutation)
      enddo
      prob=cdabs(wf)**2.d0

c ---------------------------------------------------------------------
c ---------------------- Metropolis algorithm -------------------------

      open(1,file='x.out',status='unknown')
      open(2,file='prob.out',status='unknown')

      delta=0.1d0
      do iMetropolis=1,30000   !! beginning of Metropolis algorithm

      do j=1,Natoms
        xt(j)=x(j)+delta*(2.d0*ran1(idum)-1.d0)
        xt(j)=1.d0+xt(j)-int(1.d0+xt(j))
        xs(j)=xt(j)
      enddo

      do j=Natoms,1,-1
        do i=1,j
          if(xs(i-1) .gt. xs(i)) then
            temp=xs(i)
            xs(i)=xs(i-1)
            xs(i-1)=temp
          endif
        enddo
      enddo  

      do ind=1,Natoms
        permutation(ind)=ind
      enddo
      wf=0.d0
      do indper=1,nperm  !! many-particle wave function
        temp=0.d0
        do j=1,Natoms
          temp=temp+k(permutation(j))*xs(j)
        enddo
        wf=wf+ampl(indper)*cdexp(imag*temp)
        ltemp=nextp(Natoms,permutation)
      enddo
      probnext=cdabs(wf)**2.d0
      eta=ran1(idum)
      if(probnext/prob .gt. eta) then
        write(1,*) (xt(j),j=1,Natoms)
        write(2,*) probnext
        do j=1,Natoms
          x(j)=xt(j)
        enddo
        prob=probnext
        else
        write(1,*) (x(j),j=1,Natoms)
        write(2,*) prob
      endif

      enddo   !! end of Metropolis algorithm

      close(1)
      close(2)

      end



      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*8 ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1.d0/IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2d-7,RNMX=1.d0-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END



      FUNCTION nextp(n,a)
      integer n,a,i,j,k,t
      logical nextp
      dimension a(n)

      i=n-1
   10 if(a(i).lt.a(i+1)) go to 20
      i=i-1
      if(i.eq.0) go to 20
      go to 10
   20 j=i+1
      k=n
   30 t=a(j)
      a(j)=a(k)
      a(k)=t
      j=j+1
      k=k-1
      if(j.lt.k) go to 30
      j=i
      if(j.ne.0) go to 40
      nextp=.false.
      return
   40 j=j+1
      if(a(j).lt.a(i)) go to 40
      t=a(i)
      a(i)=a(j)
      a(j)=t
      nextp=.true.

      END



      FUNCTION nextpsign(n,a,permsign)
      integer n,a,i,j,k,t,permsign
      logical nextpsign
      dimension a(n)

      i=n-1
   10 if(a(i).lt.a(i+1)) go to 20
      i=i-1
      if(i.eq.0) go to 20
      go to 10
   20 j=i+1
      k=n
   30 t=a(j)
      a(j)=a(k)
      a(k)=t
        if(j.ne.k) permsign=-permsign
      j=j+1
      k=k-1
      if(j.lt.k) go to 30
      j=i
      if(j.ne.0) go to 40
      nextpsign=.false.
      return
   40 j=j+1
      if(a(j).lt.a(i)) go to 40
      t=a(i)
      a(i)=a(j)
      a(j)=t
        if(j.ne.i) permsign=-permsign
      nextpsign=.true.

      END

