      subroutine afun(xobs,xgrid,n,m,amat,loclin,dloclin,h,d)

c 
c computation of constraint coefficient matrix (amat) for kernel 
c regression sharpening -- regression function is assumed to be 
c increasing
c 
c local constant (degree=0); local linear (degree=1)
c
c normal kernel only
c
c n = number of observations
c m = number of gridpoints
c

      implicit none
      integer d, i, j, m, n
      double precision amat(m*n), xi, xobs(n), dkernel, dloclin(n), 
     &h, loclin(n), s1, s2, s0, s0p, temp, gridpt, 
     &xgrid(m), kernel, y, zero, d2kernel
      parameter(zero=0.0d0)

      do 30 j=1,m

      s0=zero
      s1=zero
      s2=zero
      s0p=zero

      gridpt = xgrid(j)

      do 10 i=1,n
        xi = xobs(i)
        call kern(gridpt, xi, h, kernel, dkernel, d2kernel)
        loclin(i) = kernel
        dloclin(i) = dkernel
        s0 = s0+kernel 
        s1 = s1+kernel*(xi - gridpt)
        s2 = s2+kernel*(xi - gridpt)**2
        temp = dkernel
        s0p = s0p + temp
10     continue
        temp=s1**2 - s0*s2
      do 20 i=1,n    
        if (d .eq. 0) then 
          amat((j-1)*n+i) = -(dloclin(i)*s0 - loclin(i)*s0p)/s0**2
c           amat((j-1)*n+i) = loclin(i)/s0
        else
           amat((j-1)*n+i) = (s1*(xobs(i) - gridpt) - s2)*
     &loclin(i)/temp
        endif
20    continue
30    continue

      if (d .eq. 1) then
      do 50 j = 1, m-1
      do 40 i = 1, n
        amat((j-1)*n+i) = amat(j*n+i) - amat((j-1)*n + i)
40    continue
50    continue
      endif 

      return
      end

      subroutine kern(x, xi, h, kernel, dkernel, d2kernel)

      implicit none
      double precision h, x, xi, z, tmp, kernel, dkernel, d2kernel
      double precision  pi
      parameter(pi=3.14159265358979d0)

      z = (xi - x)/h

c normal 

      kernel=exp(-z**2/2.0)/(h*sqrt(2.0*pi))
      dkernel = -z*kernel/h
      d2kernel = -kernel/h**2+z**2*kernel/h**2
      return
      end

 
