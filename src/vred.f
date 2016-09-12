c# For compiling the Fortran file :
c# R CMD SHLIB vred.f


      subroutine vredpdf(ci,ip1,ip2,ipb,pdf1,pdf2,pdfb,ipar,
     1             par,model)
c     Fortran rotuine for finding the aggregated semivariance between two observations
c     Uses the pdf of distances between the two observations


      implicit double precision (a-h,o-z)
      double precision pdf1(ip1,2), pdf2(ip2,2), pdfb(ipb,2),par(ipar)
      integer ip1, ip2, ipb
      
c      write(*,*) ci,ip1,ip2,ipb,par
      w = 0
      c1 = 0
      open(unit = 23,file = "log.txt")
      do i = 1,ip1
c        write(*,*) i, (pdf1(i,j),j=1,2)
        xd = pdf1(i,1)
        gamma = vario(xd,par,model)
        w1 = pdf1(i,2)
        c1 = c1+gamma*w1
        w = w+w1
c        write(23,'(2i5,2f10.2,f13.2,2f9.0,f10.4)') 
c     1                        1,i,xd,gamma,c1,w1,w,c1/w
      enddo
      c1 = c1/w

      w = 0.
      c2 = 0.
      do i = 1,ip2
        xd = pdf2(i,1)
        gamma = vario(xd,par,model)
        w2 = pdf2(i,2)
        c2 = c2+gamma*w2
        w = w+w2
c        write(23,'(2i5,2f10.2,f13.2,2f9.0,f10.4)')  
c     1             2,i,xd,gamma,c2,w2,w,c2/w
      enddo
      c2 = c2/w


      w = 0.
      cb = 0.
      do i = 1,ipb
        xd = pdfb(i,1)
        gamma = vario(xd,par,model)
        wb = pdfb(i,2)
        cb = cb+gamma*wb
        w = w+wb
c        write(23,'(2i5,2f10.2,f13.2,2f9.0,f10.4)')  
c     1             3,i,xd,gamma,cb,wb,w,cb/w
      enddo
      cb = cb/w

      ci = cb-(c1+c2)/2.
c      write(*,*) ci, w, cb, c1, c2
      close(23)
      end




      subroutine vredhyp(ci,a1,a2,dist,ipar,par,resol,model)
c     Fortran rotuine for finding the aggregated semivariance between two observations
c     The observations are given by their areas and the distances between them.
      implicit double precision (a-h,o-z)
      double precision par(ipar)
      integer resol
      
      ar1 = sqrt(a1)
      ar2 = sqrt(a2)
      xd1 = ar1/resol
      xd2 = ar2/resol
      x01 = 0.-floor(resol/2.)*xd1
      y01 = 0.-floor(resol/2.)*xd1
      x02 = dist-floor(resol/2.)*xd2
      y02 = 0.-floor(resol/2.)*xd2
      x1max = ar1
      x2min = dist

      c1 = 0.
      c2 = 0.
      cb = 0.
c      write(*,*) ci,a1,a2,ar1, ar2, xd1,xd2,x02,dist,ipar
      do i = 1,resol
        x11 = x01 + i*xd1
        x12 = x02 + i*xd2
        do j = 1,resol
          x21 = x01 + j*xd1
          x22 = x02 + j*xd2
c          write(*,'(2i3,15f10.2)') 
c     1          i,j,x11,x12,y11,y12,x01,x02,y01,y02,xd1,xd2
          xs1 = x21 -x11
          xs2 = x22 -x12
          xsb = x22 -x11
          xdd1 = xs1*xs1
          xdd2 = xs2*xs2
          xddb = xsb*xsb
          do k = 1,resol
            y11 = y01 + k*xd1
            y12 = y02 + k*xd2
            do l = 1,resol
              y21 = y01 + j*xd1
              y22 = y02 + j*xd2
              ys1 = y21 -y11
              ys2 = y22 -y12
              ysb = y22 -y11
              ydd1 = ys1*ys1
              ydd2 = ys2*ys2
              yddb = ysb*ysb

              hd1 = sqrt(xdd1 + ydd1)
              hd2 = sqrt(xdd2 + ydd2)
              hdb = sqrt(xddb + yddb)

              g1 = vario(hd1,par,model)
              g2 = vario(hd2,par,model)
              gb = vario(hdb,par,model)

              c1 = c1+g1
              c2 = c2+g2
              cb = cb+gb
c              if(i .eq. 2 .and. j .eq. 2) then
c                write(*,*) model
c                write(*,*)  i,j,k,l,x11,x12,y11,y12
c                write(*,*) x21,x22,y21,y22
c                write(*,*) hd1,hd2,hdb,g1,g2,gb,c1,c2,cb
c              endif
            enddo
          enddo
        enddo
      enddo

      ires = resol*resol*resol*resol
      c1 = c1/ires
      c2 = c2/ires
      cb = cb/ires
      ci = cb-(c1+c2)/2.

c      write(*,*) ires,c1,c2,cb,ci
      end



c   vreda = .Fortran("vredind",ci,ip1,ip2,a1,a2,length(pars),pars,model)

      subroutine vredind(ci,ip1,ip2,a1,a2,ipar,pars,model)
      implicit double precision (a-h,o-z)
      double precision a1(ip1,2), a2(ip2,2),pars(ipar),par(20)
c      write(*,*) model,model,model
c      write(*,*) par
      ct = 0
      do i = 1,ipar
        par(i) = pars(i)
      enddo
      do i = 1,ip1
        x1 = a1(i,1)
        y1 = a1(i,2)
        do j = 1,ip2
          x2 = a2(j,1)
          y2 = a2(j,2)
          xs = x2-x1
          ys = y2-y1
          xss = xs*xs+ys*ys
          if (xss .le. 0) then
            xd = 0
          else
            xd = sqrt(xss)
          endif
c          write(*,*) i,j,x1,x2,par(1)
          xgamma = vario(xd,par,model)
c          write(*,*) i,j,x1,x2,xd,par(1),xgamma,model
          ct = ct+xgamma
        enddo
      enddo
      ci = ct/(ip1*ip2*1.)
      end




      subroutine varioex(res,skor,ip,pa,model)
c     Function for external requests for variogram value
      implicit double precision (a-h,o-z)
      dimension pa(ip),par(20)
      do i = 1,ip
        par(i) = pa(i)
      enddo
c      write(58,*) skor,ip,(par(i),i=1,ip),model,res
      res = vario(skor,par,model)
      end








      double precision function vario(skor,par,model)
c     Function choosing between the different variogram models available
c     The numbers should match the numbers of the R-function rtop:::imodel
      implicit double precision (a-h,o-z)
      dimension par(20)
c      write(*,*) model,model
c      write(*,*) skor,(par(i),i=1,5),model

      if (model .eq. 1) then
        vario = exponential(skor,par)
        return
      else if (model .eq. 2) then
        vario = exp1(skor,par)
        return
      else if (model .eq. 3) then
        vario = gaussian(skor,par)
        return
      else if (model .eq. 4) then
        vario = gau1(skor,par)
        return
      else if (model .eq. 5) then
        vario = skor
        return
      else if (model .eq. 6) then
        vario = spherical(skor,par)
        return
      else if (model .eq. 7) then
        vario = sph1(skor, par)
        return
      else if (model .eq. 8) then
        vario = fractal(skor, par)
        return
      else 
        vario = -999
        return
      endif
      end




      double precision function exponential(skor,par)
      implicit double precision (a-h,o-z)
      dimension par(20)
      xd = skor
  
      c = par(2)
      epar = (xd/c)
      exponential = par(1)*(1-exp(-epar))
      end



      double precision function exp1(skor,par)
      implicit double precision (a-h,o-z)
      dimension par(20)
      xd = skor
      rlim = 1e-7

      c = par(2)
      if (xd > rlim) then
        xfrac = (xd)**par(4)
      else 
        xfrac = 0
      endif
      epar = (xd/c)**par(5)
      exp1 = par(1)*xfrac*(1-exp(-epar)) 
      end


      double precision function gaussian(skor,par)
      implicit double precision (a-h,o-z)
      dimension par(20)
      
      xd = skor
      c = par(2)
      epar = (xd**2)/(c**2)
      gaussian = par(1)*(1-exp(-epar))
      end


      double precision function gau1(skor,par)
      implicit double precision (a-h,o-z)
      dimension par(20)
      rlim = 1e-9      
      xd = skor
      if (xd > rlim) then
        xfrac = (xd)**par(4)
      else 
        xfrac = 0
      endif
      c = par(2)
      epar = (xd**2)/(c**2)
      gau1 = par(1)*xfrac*(1-exp(-epar))
      end


      double precision function spherical(skor,par)
      implicit double precision (a-h,o-z)
      dimension par(20)

      xd = skor
      c = par(2)
      epar = (xd/c)
      spherical = par(1)*(3./2.*epar-.5*(epar**3))
      end
      
      
      
      double precision function sph1(skor,par)
      implicit double precision (a-h,o-z)
      dimension par(20)

      xd = skor
      rlim = 1e-9
      c = par(2)
      if (xd > rlim) then
        xfrac = (xd)**par(4)
      else 
        xfrac = 0
      endif
      epar = (xd/c)
      sph1 = par(1)*xfrac*(3./2.*epar-.5*(epar**3))
      end



      double precision function fractal(skor, par)
      implicit double precision (a-h,o-z)
      dimension par(20)
      xd = skor
      rlim = 1e-7
      if (xd > rlim) then
        xfrac = (xd)**par(2)
      else 
        xfrac = 0
      endif
      fractal = par(1)*xfrac
      end
