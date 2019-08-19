      subroutine afun(xobs, xgrid, n, m, amat, choice, h)

c 
c computation of constraint coefficient matrix (amat) for kernel 
c regression sharpening -- regression function is assumed to be 
c increasing
c 
c local constant (degree=0)
c
c kernel choices: 1 - biweight
c                 2 - normal
c                 3 - epanechnikov
c
c n = number of observations
c m = number of gridpoints
c

      implicit none
      integer i, j, choice, m, n
      double precision amat(m*n), xgrid(m), xobs(n), x, xo, dy, d2y, 
     &y,  h, sumkjx, sumdkjx
        
      do 30 j=1, m
        sumkjx = 0.0
        sumdkjx = 0.0
        do 10 i=1, n
          x = xgrid(j)
          xo = xobs(i)

c            y = kernel evaluated at x-xo
c           dy = kernel derivated evaluated at x-xo, etc.

          call kern(x, choice, xo, h, y, dy, d2y)
          sumkjx = sumkjx + y
          sumdkjx = sumdkjx + dy
10      continue
        do 20 i=1, n
          x = xgrid(j)
          xo = xobs(i)
          call kern(x, choice, xo, h, y, dy, d2y)
          amat((j-1)*n+i) = (dy*sumkjx - y*sumdkjx)/sumkjx**2
20      continue
30    continue

      return
      end


      subroutine kern(x, choice, xo, h, kernel, dkernel, d2kernel)

      implicit none
      integer choice   
      double precision h, x, xo, z, tmp, kernel, dkernel, d2kernel
      double precision  pi
      parameter(pi=3.14159265358979d0)

      z = (x-xo)/h

c biweight (choice 1):
      if (choice .eq. 1) then
      if (abs(z) .lt. 1.0d0) then
            kernel = 15*(1-z**2)**2/16
            dkernel =  15*z*(z**2-1)/4
            d2kernel = 15*(3*z**2-1)/4
      else
            kernel = 0.0d0
            dkernel = 0.0d0
            d2kernel = 0.0d0
      end if
        
      else

c normal (choice 2):
      if (choice .eq. 2) then
           kernel=exp(-z**2/2.0)/sqrt(2.0*pi)
           dkernel = -z*kernel
           d2kernel = -kernel-z*dkernel
      else

c epanechnikov (choice 3):
          if (choice .eq. 3) then
          if (abs(z) .lt. 1.0) then
                kernel=0.75*(1-z**2)
                dkernel = -1.5*z
                d2kernel = -1.5
          else
                kernel = 0.0
                dkernel = 0.0
                d2kernel = 0.0
          end if

          end if
          end if
      end if

        if (h .ne. 1.0) then
        kernel = kernel/h
        dkernel = dkernel/h**2
        d2kernel = d2kernel/h**3
        end if

      return
      end

 
      subroutine dllc(xobs,xgrid,n,m,amat,loclin,dloclin,choice,h,work)

c local linear nonparametric regression weights (and derivatives)
c evaluated at a point gridpt; xobs = horizontal coordinates of original data

      implicit none
      integer i, j, choice, m, n
      double precision xo, xobs(n), dkernel, dloclin(n), d2kernel,   
     &h, work(n), loclin(n), s1, s2, s0, s1p, s2p, s0p, temp, gridpt, 
     &xgrid(m), kernel, amat(m*n), zero
      parameter(zero=0.0d0)

      do 40 j=1,m

      s0=zero
      s1=zero
      s2=zero
      s1p=zero
      s2p=zero
      s0p=zero

      gridpt = xgrid(j)

      do 10 i=1,n
        xo = xobs(i)
        call kern(gridpt, choice, xo, h, kernel, dkernel, d2kernel)
        loclin(i) = kernel
        work(i)= kernel
        dloclin(i) = dkernel
        s0 = s0+kernel 
        s1 = s1+kernel*(gridpt - xobs(i))
        s2 = s2+kernel*(gridpt - xobs(i))**2
        temp = dkernel
        s0p = s0p + temp
        temp=temp*(gridpt-xobs(i))
        s1p = s1p + temp
        temp=temp*(gridpt-xobs(i))
        s2p = s2p + temp
10     continue
        s1p=s0+s1p
        s2p=2*s1+s2p
        temp=s2*s0-s1**2
      do 20 i=1,n    
        loclin(i)=(s2 - s1*(gridpt - xobs(i)))*loclin(i)/temp
       dloclin(i)=((s2p-s1p*(gridpt-xobs(i))-s1)*work(i)   
     &+(s2-s1*(gridpt-xobs(i)))*dloclin(i)-loclin(i)*
     &(s2p*s0+s0p*s2-2*s1*s1p))/temp
20    continue

      do 30 i=1,n
      amat((j-1)*n+i) = dloclin(i)
30    continue      

40    continue

      return
      end
      subroutine ll(xobs, yobs, xgrid, ygrid, n, m, choice, h, work)
      
      implicit none
      integer choice, i, j, m, n
      double precision yhat, xobs(n), yobs(n), xgrid(m), ygrid(m), h, 
     &work(m), gridpt
      
      do 20 j=1,m
      gridpt = xgrid(j)
      call llwts(xobs, gridpt, n, choice, work, h)
      yhat = 0.0
      do 10 i=1,n
      yhat = yhat + work(i)*yobs(i)
10    continue     
      ygrid(j) = yhat
20    continue 
                   
      return
      end
      
      subroutine llwts(xobs, gridpt, n, choice, llwt, h)

! local linear nonparametric regression weights
! evaluated at a point x=gridpt; data.x = horizontal coordinates of original data

      implicit none
      integer choice, i, n
      double precision xobs(n), kernel, h, llwt(n),
     &s1, s2, s0, temp, gridpt, zero, dkernel, d2kernel, xo
      parameter(zero=0.0d0)

      s0=zero
      s1=zero
      s2=zero
      do 10 i=1,n
        xo = xobs(i)
        call kern(gridpt, choice, xo, h, kernel, dkernel, d2kernel)
        llwt(i) = kernel
        
        s0 = s0+llwt(i) 
        s1 = s1+llwt(i)*(gridpt - xobs(i))
        s2 = s2+llwt(i)*(gridpt - xobs(i))**2
10    continue
        temp=s2*s0-s1**2
      do 20 i=1,n    
        llwt(i)=(s2 - s1*(gridpt - xobs(i)))*llwt(i)/temp
20    continue

      return
      end

