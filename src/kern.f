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

