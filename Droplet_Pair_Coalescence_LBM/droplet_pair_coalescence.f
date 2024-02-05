*************************************************************
      program main
*************************************************************

      implicit double precision (a-h, o-z)

      open(10, file = 'input', status = 'old')

      read(10,*) nx, ny, nz

      call solve(nx, ny, nz)

      stop
      end

*************************************************************
      subroutine solve(nx,ny,nz)
*************************************************************

      implicit double precision (a-h, o-z)

      dimension rr0(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension rc0(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension cp0(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension dif(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension dif1(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension pp0(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension ge(-1:nx+1,-1:ny+1,-1:nz+1,0:26)
      dimension ff(-1:nx+1,-1:ny+1,-1:nz+1,0:26)
      dimension gg(-1:nx+1,-1:ny+1,-1:nz+1,0:26)
      dimension dox(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension doy(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension doz(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension dcx(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension dcy(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension dcz(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension dpx(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension dpy(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension dpz(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension doxm(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension doym(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension dozm(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension dcxm(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension dcym(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension dczm(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension dpxm(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension dpym(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension dpzm(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension uu0(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension vv0(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension ww0(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension tau(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension hgt(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension wa(0:26), ic(0:26,3), io(26), dst(26)

      character filename*50, filenam1*50, num*3

      read(10,*) Re, Ca, theta
      read(10,*) rad, rmu
      read(10,*) max_iter, delta, r_ratio, interv

      open(28, file = 'radius.dat', status = 'unknown')

      RT = 1.d0/3.d0
      RTi = 3.d0

      Pi = 4.d0*atan(1.d0)

      tau0 = 0.5d0

      ind = 0
      indd = 0
      iter_s = 0

C     INITIALIZE MACROSCOPIC VARIABLES

      call weight(wa,ic,io,dst,dstm)

      do ii = 0, 26
      do k = -1, nz+1
      do j = -1, ny+1
      do i = -1, nx+1

         ff(i,j,k,ii) = 0.d0
         gg(i,j,k,ii) = 0.d0

      enddo
      enddo
      enddo
      enddo

      x0 = nx/2
      y0 = ny/2
      z0 = nz/2

! ======= Incorporating droplet pair ======================

      x1 = (nx/2) - 40
      y1 = (ny/2) + 27.5
      z1 = nz/2

      x2 = (nx/2) + 40
      y2 = ny/2 - 27.5
      z2 = nz/2

!======================================================

      rho_h = 1.0d0

      rho_l = rho_h/r_ratio

      drdc = rho_h-rho_l

      dcdr = 1.d0/drdc

      rl_dcdr = rho_l*dcdr

      ind = 0

!      uc = 0.003d0

      uc = 0.001d0

      amul = rho_l*uc*2.d0*rad/Re

      amuh = amul*rmu

      sigma = amul*uc/Ca

      interv0 = Interv

      print*, 'Interv0 = ', interv0

c      tauh = 0.02d0

      tauh = 3.d0*amuh/rho_h

      taul = 3.d0*amul/rho_l

      !taul = tauh*rmu

      rrm = 0.5d0*(rho_h+rho_l)

      rrh = 0.5d0*(rho_h-rho_l)

      rc_h = rho_h

      rc_l = 0.d0

      rcm = 0.5d0*(rc_h+rc_l)

      rch = 0.5d0*(rc_h-rc_l)

      drdc = rho_h-rho_l

      beta = 12.d0*sigma/(rc_h-rc_l)**4/delta

      akappa = beta*delta**2*(rc_h-rc_l)**2/8.d0

!      aMb = 0.05/4./beta
!      aMb = 0.05/100./beta
!       aMb = 0.00009/beta
       Pec = 500.0
       !aMb = uc*dsqrt(akappa/beta)/Pec/beta
       !aMb = 0.01/4./beta/4
       aMb = 0.8


      print*, 'Tau_h = ', tauh, ', Tau_l = ', taul
      print*, 'Rho_h = ', rho_h, ', Rho_l = ', rho_l
      print*, 'Mu_h = ', amuh, ', Mu_l = ', amul
      print*, 'Sigma = ', sigma
      print*, 'Re = ', uc*2.*rad/amuh*rho_h
      print*, 'We = ', We
      print*, 'uc = ', uc
      print*, 'aMb = ', aMb


      ini_iter = -sqrt(ti)

      ini_iter = 0

      print*, 'Beta = ', beta, ', Kappa = ', akappa

      gry0 = 0.

      Bo = 3.475E-4

      Bo = 0.d0

      gry1 = -Bo/rho_h/(2.*rad)**2*sigma

      print*, 'Gravity = ', -gry1

      do k = -2, nz+2
      do j = -2, ny+2
      do i = -2, nx+2




         if (i.LT.x0) then
            c = (dsqrt((i-x1)**2+(j-y1)**2+(k-z1)**2*0)-rad)/delta*2.d0
         else
            c = (dsqrt((i-x2)**2+(j-y2)**2+(k-z2)**2*0)-rad)/delta*2.d0
         end if



         rc0(i,j,k) = rcm-rch*tanh(c)
         pp0(i,j,k) = 0.d0
         hgt(i,j,k) = j
         cp0(i,j,k) = 0.d0
         dif(i,j,k) = 0.d0
         dif1(i,j,k) = 0.d0
         uu0(i,j,k) = 0.d0
         vv0(i,j,k) = 0.d0
         ww0(i,j,k) = 0.d0
         dox(i,j,k) = 0.d0
         doy(i,j,k) = 0.d0
         doz(i,j,k) = 0.d0
         dcx(i,j,k) = 0.d0
         dcy(i,j,k) = 0.d0
         dcz(i,j,k) = 0.d0
         dpx(i,j,k) = 0.d0
         dpy(i,j,k) = 0.d0
         dpz(i,j,k) = 0.d0
         doxm(i,j,k) = 0.d0
         doym(i,j,k) = 0.d0
         dozm(i,j,k) = 0.d0
         dcxm(i,j,k) = 0.d0
         dcym(i,j,k) = 0.d0
         dczm(i,j,k) = 0.d0
         dpxm(i,j,k) = 0.d0
         dpym(i,j,k) = 0.d0
         dpzm(i,j,k) = 0.d0

      enddo
      enddo
      enddo

      do k = 0, nz
      do j = 0, ny
      do i = 0, nx

         rr0(i,j,k) = rc0(i,j,k)*rho_h+(1.d0-rc0(i,j,k))*rho_l

         tau(i,j,k) = rc0(i,j,k)/tauh+(1.d0-rc0(i,j,k))/taul

         tau(i,j,k) = 1.d0/tau(i,j,k)

         if (j.eq.0) tau(i,j,k) = 0.5d0
         if (j.eq.ny) tau(i,j,k) = 0.5d0

      enddo
      enddo
      enddo

      call gradient(rc0,doxm,doym,dozm,dox,doy,doz,nx,ny,nz)

      call gradient(pp0,dpxm,dpym,dpzm,dpx,dpy,dpz,nx,ny,nz)

      call diffusion(rc0,dif,nx,ny,nz)

      do k = 0, nz
      do j = 0, ny
      do i = 0, nx

         cp0(i,j,k) = 4.d0*beta*rc0(i,j,k)
     *               *(rc0(i,j,k)-0.5d0)
     *               *(rc0(i,j,k)-1.0d0)
     *               -akappa*dif(i,j,k)

      enddo
      enddo
      enddo

      call diffusion(cp0,dif,nx,ny,nz)

      call gradient(cp0,dcxm,dcym,dczm,dcx,dcy,dcz,nx,ny,nz)

c     Time marching

      do ii = 0, 26

      do k = 0, nz
      do j = 0, ny
      do i = 0, nx

         gamm = 3.d0*(ic(ii,1)*uu0(i,j,k)
     *               +ic(ii,2)*vv0(i,j,k)
     *               +ic(ii,3)*ww0(i,j,k))
     *        +4.5d0*(ic(ii,1)*uu0(i,j,k)
     *               +ic(ii,2)*vv0(i,j,k)
     *               +ic(ii,3)*ww0(i,j,k))**2
     *        -1.5d0*(uu0(i,j,k)**2
     *               +vv0(i,j,k)**2
     *               +ww0(i,j,k)**2)

         gamma = wa(ii)*(1.d0+gamm)

         ge(i,j,k,ii) = gamma

         o1 = rc0(i+ic(ii,1),j+ic(ii,2),k+ic(ii,3))-rc0(i,j,k)
         o0 = rc0(i,j,k)-rc0(i-ic(ii,1),j-ic(ii,2),k-ic(ii,3))

         c1 = cp0(i+ic(ii,1),j+ic(ii,2),k+ic(ii,3))-cp0(i,j,k)
         c0 = cp0(i,j,k)-cp0(i-ic(ii,1),j-ic(ii,2),k-ic(ii,3))

         g1 = hgt(i+ic(ii,1),j+ic(ii,2),k+ic(ii,3))-hgt(i,j,k)
         g0 = hgt(i,j,k)-hgt(i-ic(ii,1),j-ic(ii,2),k-ic(ii,3))

         udo = uu0(i,j,k)*dox(i,j,k)
     *        +vv0(i,j,k)*doy(i,j,k)
     *        +ww0(i,j,k)*doz(i,j,k)

         ofc = o1-0.5d0*(o1-o0)-udo

         udc = uu0(i,j,k)*dcx(i,j,k)
     *        +vv0(i,j,k)*dcy(i,j,k)
     *        +ww0(i,j,k)*dcz(i,j,k)

         cfc = rc0(i,j,k)
     *            *(c1-0.5d0*(c1-c0)-udc)

         coef = rr0(i,j,k)-rho_l

         gfc = coef*gry0
     *        *(g1-0.5d0*(g1-g0)-vv0(i,j,k))

        ffeq = rc0(i,j,k)*ge(i,j,k,ii)
     *            -(ofc-(cfc-gfc)*RTi*rc0(i,j,k)/rr0(i,j,k))
     *            *ge(i,j,k,ii)*0.5d0

         ggeq = wa(ii)*(pp0(i,j,k)+rr0(i,j,k)*RT*gamm)
     *         -ofc*drdc*RT*(ge(i,j,k,ii)-wa(ii))*0.5d0
     *         +(cfc-gfc)*ge(i,j,k,ii)*0.5d0

         ff(i,j,k,ii) = ffeq

         gg(i,j,k,ii) = ggeq

      enddo
      enddo
      enddo

      enddo

      do 1000 iter = ini_iter, max_iter

C     COLLISION

      do ii = 0, 26

      do k = 0, nz
      do j = 0, ny
      do i = 0, nx

         gamm = 3.d0*(ic(ii,1)*uu0(i,j,k)
     *               +ic(ii,2)*vv0(i,j,k)
     *               +ic(ii,3)*ww0(i,j,k))
     *        +4.5d0*(ic(ii,1)*uu0(i,j,k)
     *               +ic(ii,2)*vv0(i,j,k)
     *               +ic(ii,3)*ww0(i,j,k))**2
     *        -1.5d0*(uu0(i,j,k)**2
     *               +vv0(i,j,k)**2
     *               +ww0(i,j,k)**2)

         gamma = wa(ii)*(1.d0+gamm)

         ge(i,j,k,ii) = gamma

         o2 = rc0(i+2*ic(ii,1),j+2*ic(ii,2),k+2*ic(ii,3))
     *       -rc0(i+ic(ii,1),j+ic(ii,2),k+ic(ii,3))
         o1 = rc0(i+ic(ii,1),j+ic(ii,2),k+ic(ii,3))-rc0(i,j,k)
         o0 = rc0(i,j,k)-rc0(i-ic(ii,1),j-ic(ii,2),k-ic(ii,3))

         c2 = cp0(i+2*ic(ii,1),j+2*ic(ii,2),k+2*ic(ii,3))
     *       -cp0(i+ic(ii,1),j+ic(ii,2),k+ic(ii,3))
         c1 = cp0(i+ic(ii,1),j+ic(ii,2),k+ic(ii,3))-cp0(i,j,k)
         c0 = cp0(i,j,k)-cp0(i-ic(ii,1),j-ic(ii,2),k-ic(ii,3))

         p2 = pp0(i+2*ic(ii,1),j+2*ic(ii,2),k+2*ic(ii,3))
     *       -pp0(i+ic(ii,1),j+ic(ii,2),k+ic(ii,3))
         p1 = pp0(i+ic(ii,1),j+ic(ii,2),k+ic(ii,3))-pp0(i,j,k)
         p0 = pp0(i,j,k)-pp0(i-ic(ii,1),j-ic(ii,2),k-ic(ii,3))

         g2 = hgt(i+2*ic(ii,1),j+2*ic(ii,2),k+2*ic(ii,3))
     *       -hgt(i+ic(ii,1),j+ic(ii,2),k+ic(ii,3))
         g1 = hgt(i+ic(ii,1),j+ic(ii,2),k+ic(ii,3))-hgt(i,j,k)
         g0 = hgt(i,j,k)-hgt(i-ic(ii,1),j-ic(ii,2),k-ic(ii,3))

         udom = uu0(i,j,k)*doxm(i,j,k)
     *         +vv0(i,j,k)*doym(i,j,k)
     *         +ww0(i,j,k)*dozm(i,j,k)

         udo = uu0(i,j,k)*dox(i,j,k)
     *        +vv0(i,j,k)*doy(i,j,k)
     *        +ww0(i,j,k)*doz(i,j,k)

         ofu = o1-0.25d0*(o2-o0)-udom

         ofc = o1-0.5d0*(o1-o0)-udo

         udcm = uu0(i,j,k)*dcxm(i,j,k)
     *         +vv0(i,j,k)*dcym(i,j,k)
     *         +ww0(i,j,k)*dczm(i,j,k)

         udc = uu0(i,j,k)*dcx(i,j,k)
     *        +vv0(i,j,k)*dcy(i,j,k)
     *        +ww0(i,j,k)*dcz(i,j,k)

         cfu = rc0(i,j,k)
     *        *(c1-0.25d0*(c2-c0)-udcm)

        cfc = rc0(i,j,k)
     *        *(c1-0.5d0*(c1-c0)-udc)

         udpm = uu0(i,j,k)*dpxm(i,j,k)
     *         +vv0(i,j,k)*dpym(i,j,k)
     *         +ww0(i,j,k)*dpzm(i,j,k)

         udp = uu0(i,j,k)*dpx(i,j,k)
     *        +vv0(i,j,k)*dpy(i,j,k)
     *        +ww0(i,j,k)*dpz(i,j,k)

         pfu = p1-0.25d0*(p2-p0)-udpm

         pfc = p1-0.5d0*(p1-p0)-udp

         coef = rr0(i,j,k)-rho_l

         gfu = (coef*gry0+rr0(i,j,k)*gry1)
     *        *(g1-0.25d0*(g2-g0)-vv0(i,j,k))

         gfc = (coef*gry0+rr0(i,j,k)*gry1)
     *        *(g1-0.5d0*(g1-g0)-vv0(i,j,k))

         ffeq = rc0(i,j,k)*ge(i,j,k,ii)
     *         -(ofc
     *         -(cfc-gfc)*RTi*rc0(i,j,k)/rr0(i,j,k))
     *          *ge(i,j,k,ii)*0.5d0

         ggeq = wa(ii)*(pp0(i,j,k)+rr0(i,j,k)*RT*gamm)
     *         -(ofc*drdc*RT)*(ge(i,j,k,ii)-wa(ii))*0.5d0
     *         +(cfc-gfc)*ge(i,j,k,ii)*0.5d0

         ff(i,j,k,ii) = ffeq
     *                +(ofu
     *                -(cfu-gfu)*RTi*rc0(i,j,k)/rr0(i,j,k))
     *                *ge(i,j,k,ii)
     *                +0.5d0*aMb*dif(i,j,k)*ge(i,j,k,ii)

         gg(i,j,k,ii) = gg(i,j,k,ii)
     *                -(gg(i,j,k,ii)-ggeq)/(tau(i,j,k)+0.5d0)
     *                +(ofu*drdc*RT)*(ge(i,j,k,ii)-wa(ii))
     *                -(cfu-gfu)*ge(i,j,k,ii)

      enddo
      enddo
      enddo

      enddo

C     STREAMING

      call PS(ff,ic,nx,ny,nz)  ! Perfect shift step
      call PS(gg,ic,nx,ny,nz)

      do ii = 0, 26

      do k = 0, nz
      do j = 0, ny
      do i = 0, nx

         ff(i,j,k,ii) = ff(i,j,k,ii)
     *                 +0.5d0*aMb*dif(i,j,k)*ge(i,j,k,ii)

      enddo
      enddo
      enddo

      enddo

C     BOUNCE-BACK

      do ii = 1, 26
      do k = 0, nz
      do j = 0, ny, ny
      do i = 0, nx

         if ((j-ic(ii,2)).lt.0) then

         ff(i,j,k,ii) = ff(i,j,k,io(ii))

         gg(i,j,k,ii) = gg(i,j,k,io(ii))

         endif

      enddo
      enddo
      enddo
      enddo

C     RECOVER MACROSCOPIC VARIABLES

C     RHO

      do k = 0, nz
      do j = 0, ny
      do i = 0, nx

         rc0(i,j,k) = 0.d0

      enddo
      enddo
      enddo

      do ii = 0, 26
      do k = 0, nz
      do j = 0, ny
      do i = 0, nx

         rc0(i,j,k) = rc0(i,j,k)+ff(i,j,k,ii)

      enddo
      enddo
      enddo
      enddo

      call diffusion(rc0,dif,nx,ny,nz)

      do k = 0, nz
      do j = 0, ny
      do i = 0, nx

         cp0(i,j,k) = 4.d0*beta*rc0(i,j,k)
     *               *(rc0(i,j,k)-0.5d0)
     *               *(rc0(i,j,k)-1.0d0)
     *               -akappa*dif(i,j,k)

      enddo
      enddo
      enddo

      call diffusion(cp0,dif,nx,ny,nz)

      call gradient(cp0,dcxm,dcym,dczm,dcx,dcy,dcz,nx,ny,nz)

      call gradient(rc0,doxm,doym,dozm,dox,doy,doz,nx,ny,nz)

      do k = 0, nz
      do j = 0, ny
      do i = 0, nx

         rr0(i,j,k) = rc0(i,j,k)*rho_h+(1.d0-rc0(i,j,k))*rho_l

         tau(i,j,k) = rc0(i,j,k)/tauh+(1.d0-rc0(i,j,k))/taul

         tau(i,j,k) = 1.d0/tau(i,j,k)

      enddo
      enddo
      enddo

      do k = 0, nz
      do j = 0, ny, ny
      do i = 0, nx

         tau(i,j,k) = 0.5d0

      enddo
      enddo
      enddo

C     VELOCITIES

      do k = 0, nz
      do j = 0, ny
      do i = 0, nx

         uu0(i,j,k) = 0.d0
         vv0(i,j,k) = 0.d0
         ww0(i,j,k) = 0.d0

      enddo
      enddo
      enddo

      do ii = 1, 26

      do k = 0, nz
      do j = 0, ny
      do i = 0, nx

         uu0(i,j,k) = uu0(i,j,k)+ic(ii,1)*gg(i,j,k,ii)
         vv0(i,j,k) = vv0(i,j,k)+ic(ii,2)*gg(i,j,k,ii)
         ww0(i,j,k) = ww0(i,j,k)+ic(ii,3)*gg(i,j,k,ii)

      enddo
      enddo
      enddo

      enddo

      do k = 0, nz
      do j = 0, ny
      do i = 0, nx

         coef = rr0(i,j,k)-rho_l

         uu0(i,j,k) = uu0(i,j,k)/rr0(i,j,k)/RT
     *               -0.5d0*rc0(i,j,k)*dcx(i,j,k)
     *                     /rr0(i,j,k)

         vv0(i,j,k) = vv0(i,j,k)/rr0(i,j,k)/RT
     *               -0.5d0*rc0(i,j,k)*dcy(i,j,k)
     *                     /rr0(i,j,k)
     *               +0.5d0*(coef*gry0+rr0(i,j,k)*gry1)
     *                     /rr0(i,j,k)

         ww0(i,j,k) = ww0(i,j,k)/rr0(i,j,k)/RT
     *               -0.5d0*rc0(i,j,k)*dcz(i,j,k)
     *                     /rr0(i,j,k)

      enddo
      enddo
      enddo

      do k = 0, nz   ! the velocities we apply to top and bottom wall
      do i = 0, nx

         uu0(i,0,k) = -uc
         uu0(i,ny,k) = uc
         vv0(i,0,k) = 0.d0
         vv0(i,ny,k) = 0.d0
         ww0(i,0,k) = 0.d0
         ww0(i,ny,k) = 0.d0

      enddo
      enddo

C     PRESSURE

      do k = 0, nz
      do j = 0, ny
      do i = 0, nx

         pp0(i,j,k) = 0.d0

      enddo
      enddo
      enddo

      do ii = 0, 26

      do k = 0, nz
      do j = 0, ny
      do i = 0, nx

         pp0(i,j,k) = pp0(i,j,k)+gg(i,j,k,ii)

      enddo
      enddo
      enddo

      enddo

      do k = 0, nz
      do j = 0, ny
      do i = 0, nx

         pp0(i,j,k) = pp0(i,j,k)
     *               +0.5d0*RT*(uu0(i,j,k)*dox(i,j,k)
     *                         +vv0(i,j,k)*doy(i,j,k)
     *                         +ww0(i,j,k)*doz(i,j,k))
     *                        *drdc

      enddo
      enddo
      enddo

      call gradient(pp0,dpxm,dpym,dpzm,dpx,dpy,dpz,nx,ny,nz)

C     POSTPROCESSING

      if (mod(iter,10).eq.0) then
      print*, 'Iter = ', iter
      endif

      if (mod(iter,interv).eq.0) then

c     Write

      md = iter/interv

      write(num,'(i3.3)') md

!      filename = 'uns'//num(1:3)//'.dat'

!      open(11, file = filename, status = 'unknown')


!      nx0 = nx
!      ny0 = ny
!      nz0 = nz

!      call diffusion(rc0,dif1,nx,ny,nz)

!      write(11,*) 'variables="x","y","z","r","u","v","w","ph","pt"'
!      write(11,*) 'zone I=',nx0+1, ', J=',ny0+1, ', K=',2*nz0+1

!      do k = -nz0, nz0
!      do j = 0, ny0
!      do i = 0, nx0

!        il = abs(i)
!         jl = abs(j)
!         kl = abs(k)

!         im = sign(1,i)
!         jm = sign(1,j)
!         km = sign(1,k)

!         pph = pp0(il,jl,kl)

!         ppt = 4.d0*beta*rc0(il,jl,kl)**2
!     *                  *(rc0(il,jl,kl)-0.5d0)
!     *                  *(rc0(il,jl,kl)-1.d0)
!     *        -beta*rc0(il,jl,kl)**2
!     *             *(rc0(il,jl,kl)-1.d0)**2
!     *        -akappa*rc0(il,jl,kl)*dif1(il,jl,kl)
!     *        +0.5d0*akappa*(dox(il,jl,kl)**2
!     *                      +doy(il,jl,kl)**2
!     *                      +doz(il,jl,kl)**2)

!         write(11,203) i/rad,j/rad,k/rad,rr0(il,jl,kl),
!     *                 im*uu0(il,jl,kl),jm*vv0(il,jl,kl),
!     *                 km*ww0(il,jl,kl),pph,ppt

!      enddo
!      enddo
!      enddo

      filenam1 = 'u2d'//num(1:3)//'.dat'

      open(16, file = filenam1, status = 'unknown')

      nx0 = nx
      ny0 = ny

      write(16,*) 'variables="x","y","r","u","v","ph","pt"'
      write(16,*) 'zone I=',nx0+1, ', J=',ny0+1

      do j = 0, ny0
      do i = 0, nx0

         il = abs(i)
         jl = abs(j)

         im = sign(1,i)
         jm = sign(1,j)

         pph = pp0(il,jl,0)

         ppt = 4.d0*beta*rc0(il,jl,0)**2
     *                  *(rc0(il,jl,0)-0.5d0)
     *                  *(rc0(il,jl,0)-1.d0)
     *        -beta*rc0(il,jl,0)**2
     *             *(rc0(il,jl,0)-1.d0)**2
     *        -akappa*rc0(il,jl,0)*dif1(il,jl,0)
     *        +0.5d0*akappa*(dox(il,jl,0)**2
     *                      +doy(il,jl,0)**2
     *                      +doz(il,jl,0)**2)

         write(16,204) i/rad,j/rad,rc0(il,jl,0),
     *                 im*uu0(il,jl,0),jm*vv0(il,jl,0),
     *                 pph,ppt

      enddo
      enddo

203   format(9e12.4)
204   format(7e12.4)

!      close(11)
      close(16)
      print*, 'Results Written!'
      print*

      endif
1000  continue

      close(28)

      return
      end



**********************************************************
      subroutine bc(p,nx,ny,nz)
**********************************************************

      implicit double precision (a-h, o-z)
      dimension p(-2:nx+2,-2:ny+2,-2:nz+2)

C     SYMMETRIC BC

C     Y = 0 & Y = NY PLANES

      do k = 0, nz
      do i = 0, nx

         p(i,-1,k) = p(i,1,k)
         p(i,-2,k) = p(i,2,k)
         p(i,ny+1,k) = p(i,ny-1,k)
         p(i,ny+2,k) = p(i,ny-2,k)

      enddo
      enddo

C     Z = 0 & Z = NZ PLANES

      do j = -2, ny+2
      do i = 0, nx

         p(i,j,-1) = p(i,j,1)
         p(i,j,-2) = p(i,j,2)
         p(i,j,nz+1) = p(i,j,nz-1)
         p(i,j,nz+2) = p(i,j,nz-2)

      enddo
      enddo

C     X = 0 & X = NX PLANES

      do k = -2, nz+2
      do j = -2, ny+2

         p(-1,j,k) = p(nx-1,j,k)
         p(-2,j,k) = p(nx-2,j,k)
         p(nx+1,j,k) = p(1,j,k)
         p(nx+2,j,k) = p(2,j,k)

      enddo
      enddo

      return
      end



**********************************************************
      subroutine gradient(p,dxm,dym,dzm,dx,dy,dz,nx,ny,nz)
**********************************************************

      implicit double precision (a-h, o-z)
      dimension p(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension dxm(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension dym(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension dzm(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension dx(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension dy(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension dz(-2:nx+2,-2:ny+2,-2:nz+2)

C     SYMMETRIC BC

C     Y = 0 & Y = NY PLANES

      do k = 0, nz
      do i = 0, nx

         p(i,-1,k) = p(i,1,k)
         p(i,-2,k) = p(i,2,k)
         p(i,ny+1,k) = p(i,ny-1,k)
         p(i,ny+2,k) = p(i,ny-2,k)
c         p(i,ny+1,k) = p(i,ny,k)
c         p(i,ny+2,k) = p(i,ny,k)

      enddo
      enddo

C     Z = 0 & Z = NZ PLANES

      do j = -2, ny+2
      do i = 0, nx

         p(i,j,-1) = p(i,j,1)
         p(i,j,-2) = p(i,j,2)
         p(i,j,nz+1) = p(i,j,nz-1)
         p(i,j,nz+2) = p(i,j,nz-2)
c         p(i,j,nz+1) = p(i,j,nz)
c         p(i,j,nz+2) = p(i,j,nz)

      enddo
      enddo

C     X = 0 & X = NX PLANES

      do k = -2, nz+2
      do j = -2, ny+2

         p(-1,j,k) = p(nx-1,j,k)
         p(-2,j,k) = p(nx-2,j,k)
         p(nx+1,j,k) = p(1,j,k)
         p(nx+2,j,k) = p(2,j,k)
c         p(nx+1,j,k) = p(nx,j,k)
c         p(nx+2,j,k) = p(nx,j,k)

      enddo
      enddo

      c1 = 2.d0/9.d0
      c2 = 1.d0/18.d0
      c3 = 1.d0/72.d0

      c1m = 1.d0/9.d0
      c2m = 1.d0/36.d0
      c3m = 1.d0/144.d0

      do k = 0, nz
      do j = 0, ny
      do i = 0, nx

         dxm(i,j,k) = (-p(i+2,j,k)+4.d0*p(i+1,j,k)
     *                 -4.d0*p(i-1,j,k)+p(i-2,j,k))*c1m
     *               +(-p(i+2,j+2,k)+4.d0*p(i+1,j+1,k)
     *                 -4.d0*p(i-1,j-1,k)+p(i-2,j-2,k)
     *                 -p(i+2,j-2,k)+4.d0*p(i+1,j-1,k)
     *                 -4.d0*p(i-1,j+1,k)+p(i-2,j+2,k)
     *                 -p(i+2,j,k+2)+4.d0*p(i+1,j,k+1)
     *                 -4.d0*p(i-1,j,k-1)+p(i-2,j,k-2)
     *                 -p(i+2,j,k-2)+4.d0*p(i+1,j,k-1)
     *                 -4.d0*p(i-1,j,k+1)+p(i-2,j,k+2))*c2m
     *               +(-p(i+2,j+2,k+2)+4.d0*p(i+1,j+1,k+1)
     *                 -4.d0*p(i-1,j-1,k-1)+p(i-2,j-2,k-2)
     *                 -p(i+2,j-2,k+2)+4.d0*p(i+1,j-1,k+1)
     *                 -4.d0*p(i-1,j+1,k-1)+p(i-2,j+2,k-2)
     *                 -p(i+2,j+2,k-2)+4.d0*p(i+1,j+1,k-1)
     *                 -4.d0*p(i-1,j-1,k+1)+p(i-2,j-2,k+2)
     *                 -p(i+2,j-2,k-2)+4.d0*p(i+1,j-1,k-1)
     *                 -4.d0*p(i-1,j+1,k+1)+p(i-2,j+2,k+2))*c3m

         dym(i,j,k) = (-p(i,j+2,k)+4.d0*p(i,j+1,k)
     *                 -4.d0*p(i,j-1,k)+p(i,j-2,k))*c1m
     *               +(-p(i+2,j+2,k)+4.d0*p(i+1,j+1,k)
     *                 -4.d0*p(i-1,j-1,k)+p(i-2,j-2,k)
     *                 -p(i-2,j+2,k)+4.d0*p(i-1,j+1,k)
     *                 -4.d0*p(i+1,j-1,k)+p(i+2,j-2,k)
     *                 -p(i,j+2,k+2)+4.d0*p(i,j+1,k+1)
     *                 -4.d0*p(i,j-1,k-1)+p(i,j-2,k-2)
     *                 -p(i,j+2,k-2)+4.d0*p(i,j+1,k-1)
     *                 -4.d0*p(i,j-1,k+1)+p(i,j-2,k+2))*c2m
     *               +(-p(i+2,j+2,k+2)+4.d0*p(i+1,j+1,k+1)
     *                 -4.d0*p(i-1,j-1,k-1)+p(i-2,j-2,k-2)
     *                 -p(i+2,j+2,k-2)+4.d0*p(i+1,j+1,k-1)
     *                 -4.d0*p(i-1,j-1,k+1)+p(i-2,j-2,k+2)
     *                 -p(i-2,j+2,k+2)+4.d0*p(i-1,j+1,k+1)
     *                 -4.d0*p(i+1,j-1,k-1)+p(i+2,j-2,k-2)
     *                 -p(i-2,j+2,k-2)+4.d0*p(i-1,j+1,k-1)
     *                 -4.d0*p(i+1,j-1,k+1)+p(i+2,j-2,k+2))*c3m

         dzm(i,j,k) = (-p(i,j,k+2)+4.d0*p(i,j,k+1)
     *                 -4.d0*p(i,j,k-1)+p(i,j,k-2))*c1m
     *               +(-p(i+2,j,k+2)+4.d0*p(i+1,j,k+1)
     *                 -4.d0*p(i-1,j,k-1)+p(i-2,j,k-2)
     *                 -p(i-2,j,k+2)+4.d0*p(i-1,j,k+1)
     *                 -4.d0*p(i+1,j,k-1)+p(i+2,j,k-2)
     *                 -p(i,j+2,k+2)+4.d0*p(i,j+1,k+1)
     *                 -4.d0*p(i,j-1,k-1)+p(i,j-2,k-2)
     *                 -p(i,j-2,k+2)+4.d0*p(i,j-1,k+1)
     *                 -4.d0*p(i,j+1,k-1)+p(i,j+2,k-2))*c2m
     *               +(-p(i+2,j+2,k+2)+4.d0*p(i+1,j+1,k+1)
     *                 -4.d0*p(i-1,j-1,k-1)+p(i-2,j-2,k-2)
     *                 -p(i+2,j-2,k+2)+4.d0*p(i+1,j-1,k+1)
     *                 -4.d0*p(i-1,j+1,k-1)+p(i-2,j+2,k-2)
     *                 -p(i-2,j+2,k+2)+4.d0*p(i-1,j+1,k+1)
     *                 -4.d0*p(i+1,j-1,k-1)+p(i+2,j-2,k-2)
     *                 -p(i-2,j-2,k+2)+4.d0*p(i-1,j-1,k+1)
     *                 -4.d0*p(i+1,j+1,k-1)+p(i+2,j+2,k-2))*c3m

         dx(i,j,k) = (p(i+1,j,k)-p(i-1,j,k))*c1
     *              +(p(i+1,j+1,k)-p(i-1,j-1,k)
     *               +p(i+1,j-1,k)-p(i-1,j+1,k)
     *               +p(i+1,j,k+1)-p(i-1,j,k-1)
     *               +p(i+1,j,k-1)-p(i-1,j,k+1))*c2
     *              +(p(i+1,j+1,k+1)-p(i-1,j-1,k-1)
     *               +p(i+1,j-1,k+1)-p(i-1,j+1,k-1)
     *               +p(i+1,j+1,k-1)-p(i-1,j-1,k+1)
     *               +p(i+1,j-1,k-1)-p(i-1,j+1,k+1))*c3

         dy(i,j,k) = (p(i,j+1,k)-p(i,j-1,k))*c1
     *              +(p(i+1,j+1,k)-p(i-1,j-1,k)
     *               +p(i-1,j+1,k)-p(i+1,j-1,k)
     *               +p(i,j+1,k+1)-p(i,j-1,k-1)
     *               +p(i,j+1,k-1)-p(i,j-1,k+1))*c2
     *              +(p(i+1,j+1,k+1)-p(i-1,j-1,k-1)
     *               +p(i+1,j+1,k-1)-p(i-1,j-1,k+1)
     *               +p(i-1,j+1,k+1)-p(i+1,j-1,k-1)
     *               +p(i-1,j+1,k-1)-p(i+1,j-1,k+1))*c3

         dz(i,j,k) = (p(i,j,k+1)-p(i,j,k-1))*c1
     *              +(p(i+1,j,k+1)-p(i-1,j,k-1)
     *               +p(i-1,j,k+1)-p(i+1,j,k-1)
     *               +p(i,j+1,k+1)-p(i,j-1,k-1)
     *               +p(i,j-1,k+1)-p(i,j+1,k-1))*c2
     *              +(p(i+1,j+1,k+1)-p(i-1,j-1,k-1)
     *               +p(i+1,j-1,k+1)-p(i-1,j+1,k-1)
     *               +p(i-1,j+1,k+1)-p(i+1,j-1,k-1)
     *               +p(i-1,j-1,k+1)-p(i+1,j+1,k-1))*c3

      enddo
      enddo
      enddo

      do k = 0, nz
      do j = 0, ny
      do i = 0, nx

         dxm(i,j,k) = 0.5d0*(dxm(i,j,k)+dx(i,j,k))

         dym(i,j,k) = 0.5d0*(dym(i,j,k)+dy(i,j,k))

         dzm(i,j,k) = 0.5d0*(dzm(i,j,k)+dz(i,j,k))

      enddo
      enddo
      enddo

      return
      end




******************************************************************
      subroutine diffusion(p,dif,nx,ny,nz)
******************************************************************

      implicit double precision (a-h, o-z)
      dimension p(-2:nx+2,-2:ny+2,-2:nz+2)
      dimension dif(-2:nx+2,-2:ny+2,-2:nz+2)

C     SYMMETRIC BC

C     Y = 0 & Y = NY PLANES

      do k = 0, nz
      do i = 0, nx

         p(i,-1,k) = p(i,1,k)
         p(i,-2,k) = p(i,2,k)
         p(i,ny+1,k) = p(i,ny-1,k)
         p(i,ny+2,k) = p(i,ny-2,k)
c         p(i,ny+1,k) = p(i,ny,k)
c         p(i,ny+2,k) = p(i,ny,k)

      enddo
      enddo

C     Z = 0 & Z = NZ PLANES

      do j = -2, ny+2
      do i = 0, nx

         p(i,j,-1) = p(i,j,1)
         p(i,j,-2) = p(i,j,2)
         p(i,j,nz+1) = p(i,j,nz-1)
         p(i,j,nz+2) = p(i,j,nz-2)
c         p(i,j,nz+1) = p(i,j,nz)
c         p(i,j,nz+2) = p(i,j,nz)

      enddo
      enddo

C     X = 0 & X = NX PLANES

      do k = -2, nz+2
      do j = -2, ny+2

         p(-1,j,k) = p(nx-1,j,k)
         p(-2,j,k) = p(nx-2,j,k)
         p(nx+1,j,k) = p(1,j,k)
         p(nx+2,j,k) = p(2,j,k)
c         p(nx+1,j,k) = p(nx,j,k)
c         p(nx+2,j,k) = p(nx,j,k)

      enddo
      enddo

      c1 = 4.d0/9.d0
      c2 = 1.d0/9.d0
      c3 = 1.d0/36.d0
      c4 = 38d0/9.d0

C     DIFFUSION

      do k = 0, nz
      do j = 0, ny
      do i = 0, nx

         dif(i,j,k) = (p(i+1,j,k)+p(i-1,j,k)
     *                +p(i,j+1,k)+p(i,j-1,k)
     *                +p(i,j,k+1)+p(i,j,k-1))*c1
     *               +(p(i+1,j+1,k)+p(i-1,j-1,k)
     *                +p(i+1,j-1,k)+p(i-1,j+1,k)
     *                +p(i+1,j,k+1)+p(i-1,j,k-1)
     *                +p(i+1,j,k-1)+p(i-1,j,k+1)
     *                +p(i,j+1,k+1)+p(i,j-1,k-1)
     *                +p(i,j+1,k-1)+p(i,j-1,k+1))*c2
     *               +(p(i+1,j+1,k+1)+p(i+1,j+1,k-1)
     *                +p(i+1,j-1,k+1)+p(i+1,j-1,k-1)
     *                +p(i-1,j+1,k+1)+p(i-1,j+1,k-1)
     *                +p(i-1,j-1,k+1)+p(i-1,j-1,k-1))*c3
     *                -p(i,j,k)*c4

      enddo
      enddo
      enddo

      return
      end




********************************************************
      subroutine PS(ff,ic,nx,ny,nz)
********************************************************

      implicit double precision (a-h, o-z)

      dimension ff(-1:nx+1,-1:ny+1,-1:nz+1,0:26)
      dimension ft(0:nx,0:ny,0:nz)
      dimension ic(0:26,3)

C     SYMMETRIC BC

C     Y = 0 & Y = NY PLANES

      do k = 0, nz
      do i = 0, nx

         ff(i,-1,k,0) = ff(i,1,k,0)
         ff(i,-1,k,1) = ff(i,1,k,1)
         ff(i,-1,k,2) = ff(i,1,k,2)
         ff(i,-1,k,3) = ff(i,1,k,4)
         ff(i,-1,k,4) = ff(i,1,k,3)
         ff(i,-1,k,5) = ff(i,1,k,5)
         ff(i,-1,k,6) = ff(i,1,k,6)
         ff(i,-1,k,7) = ff(i,1,k,9)
         ff(i,-1,k,8) = ff(i,1,k,10)
         ff(i,-1,k,9) = ff(i,1,k,7)
         ff(i,-1,k,10) = ff(i,1,k,8)
         ff(i,-1,k,11) = ff(i,1,k,11)
         ff(i,-1,k,12) = ff(i,1,k,12)
         ff(i,-1,k,13) = ff(i,1,k,13)
         ff(i,-1,k,14) = ff(i,1,k,14)
         ff(i,-1,k,15) = ff(i,1,k,18)
         ff(i,-1,k,16) = ff(i,1,k,17)
         ff(i,-1,k,17) = ff(i,1,k,16)
         ff(i,-1,k,18) = ff(i,1,k,15)
         ff(i,-1,k,19) = ff(i,1,k,21)
         ff(i,-1,k,20) = ff(i,1,k,22)
         ff(i,-1,k,21) = ff(i,1,k,19)
         ff(i,-1,k,22) = ff(i,1,k,20)
         ff(i,-1,k,23) = ff(i,1,k,25)
         ff(i,-1,k,24) = ff(i,1,k,26)
         ff(i,-1,k,25) = ff(i,1,k,23)
         ff(i,-1,k,26) = ff(i,1,k,24)

         ff(i,ny+1,k,0) = ff(i,ny-1,k,0)
         ff(i,ny+1,k,1) = ff(i,ny-1,k,1)
         ff(i,ny+1,k,2) = ff(i,ny-1,k,2)
         ff(i,ny+1,k,3) = ff(i,ny-1,k,4)
         ff(i,ny+1,k,4) = ff(i,ny-1,k,3)
         ff(i,ny+1,k,5) = ff(i,ny-1,k,5)
         ff(i,ny+1,k,6) = ff(i,ny-1,k,6)
         ff(i,ny+1,k,7) = ff(i,ny-1,k,9)
         ff(i,ny+1,k,8) = ff(i,ny-1,k,10)
         ff(i,ny+1,k,9) = ff(i,ny-1,k,7)
         ff(i,ny+1,k,10) = ff(i,ny-1,k,8)
         ff(i,ny+1,k,11) = ff(i,ny-1,k,11)
         ff(i,ny+1,k,12) = ff(i,ny-1,k,12)
         ff(i,ny+1,k,13) = ff(i,ny-1,k,13)
         ff(i,ny+1,k,14) = ff(i,ny-1,k,14)
         ff(i,ny+1,k,15) = ff(i,ny-1,k,18)
         ff(i,ny+1,k,16) = ff(i,ny-1,k,17)
         ff(i,ny+1,k,17) = ff(i,ny-1,k,16)
         ff(i,ny+1,k,18) = ff(i,ny-1,k,15)
         ff(i,ny+1,k,19) = ff(i,ny-1,k,21)
         ff(i,ny+1,k,20) = ff(i,ny-1,k,22)
         ff(i,ny+1,k,21) = ff(i,ny-1,k,19)
         ff(i,ny+1,k,22) = ff(i,ny-1,k,20)
         ff(i,ny+1,k,23) = ff(i,ny-1,k,25)
         ff(i,ny+1,k,24) = ff(i,ny-1,k,26)
         ff(i,ny+1,k,25) = ff(i,ny-1,k,23)
         ff(i,ny+1,k,26) = ff(i,ny-1,k,24)

      enddo
      enddo

C     Z = 0 & Z = NZ PLANES

      do j = -1, ny+1
      do i = 0, nx

         ff(i,j,-1,0) = ff(i,j,1,0)
         ff(i,j,-1,1) = ff(i,j,1,1)
         ff(i,j,-1,2) = ff(i,j,1,2)
         ff(i,j,-1,3) = ff(i,j,1,3)
         ff(i,j,-1,4) = ff(i,j,1,4)
         ff(i,j,-1,5) = ff(i,j,1,6)
         ff(i,j,-1,6) = ff(i,j,1,5)
         ff(i,j,-1,7) = ff(i,j,1,7)
         ff(i,j,-1,8) = ff(i,j,1,8)
         ff(i,j,-1,9) = ff(i,j,1,9)
         ff(i,j,-1,10) = ff(i,j,1,10)
         ff(i,j,-1,11) = ff(i,j,1,13)
         ff(i,j,-1,12) = ff(i,j,1,14)
         ff(i,j,-1,13) = ff(i,j,1,11)
         ff(i,j,-1,14) = ff(i,j,1,12)
         ff(i,j,-1,15) = ff(i,j,1,17)
         ff(i,j,-1,16) = ff(i,j,1,18)
         ff(i,j,-1,17) = ff(i,j,1,15)
         ff(i,j,-1,18) = ff(i,j,1,16)
         ff(i,j,-1,19) = ff(i,j,1,23)
         ff(i,j,-1,20) = ff(i,j,1,24)
         ff(i,j,-1,21) = ff(i,j,1,25)
         ff(i,j,-1,22) = ff(i,j,1,26)
         ff(i,j,-1,23) = ff(i,j,1,19)
         ff(i,j,-1,24) = ff(i,j,1,20)
         ff(i,j,-1,25) = ff(i,j,1,21)
         ff(i,j,-1,26) = ff(i,j,1,22)

         ff(i,j,nz+1,0) = ff(i,j,nz-1,0)
         ff(i,j,nz+1,1) = ff(i,j,nz-1,1)
         ff(i,j,nz+1,2) = ff(i,j,nz-1,2)
         ff(i,j,nz+1,3) = ff(i,j,nz-1,3)
         ff(i,j,nz+1,4) = ff(i,j,nz-1,4)
         ff(i,j,nz+1,5) = ff(i,j,nz-1,6)
         ff(i,j,nz+1,6) = ff(i,j,nz-1,5)
         ff(i,j,nz+1,7) = ff(i,j,nz-1,7)
         ff(i,j,nz+1,8) = ff(i,j,nz-1,8)
         ff(i,j,nz+1,9) = ff(i,j,nz-1,9)
         ff(i,j,nz+1,10) = ff(i,j,nz-1,10)
         ff(i,j,nz+1,11) = ff(i,j,nz-1,13)
         ff(i,j,nz+1,12) = ff(i,j,nz-1,14)
         ff(i,j,nz+1,13) = ff(i,j,nz-1,11)
         ff(i,j,nz+1,14) = ff(i,j,nz-1,12)
         ff(i,j,nz+1,15) = ff(i,j,nz-1,17)
         ff(i,j,nz+1,16) = ff(i,j,nz-1,18)
         ff(i,j,nz+1,17) = ff(i,j,nz-1,15)
         ff(i,j,nz+1,18) = ff(i,j,nz-1,16)
         ff(i,j,nz+1,19) = ff(i,j,nz-1,23)
         ff(i,j,nz+1,20) = ff(i,j,nz-1,24)
         ff(i,j,nz+1,21) = ff(i,j,nz-1,25)
         ff(i,j,nz+1,22) = ff(i,j,nz-1,26)
         ff(i,j,nz+1,23) = ff(i,j,nz-1,19)
         ff(i,j,nz+1,24) = ff(i,j,nz-1,20)
         ff(i,j,nz+1,25) = ff(i,j,nz-1,21)
         ff(i,j,nz+1,26) = ff(i,j,nz-1,22)

      enddo
      enddo

C     X = 0 & X = NX PLANES

      do k = -1, nz+1
      do j = -1, ny+1

         ff(-1,j,k,0) = ff(1,j,k,0)
         ff(-1,j,k,1) = ff(1,j,k,2)
         ff(-1,j,k,2) = ff(1,j,k,1)
         ff(-1,j,k,3) = ff(1,j,k,3)
         ff(-1,j,k,4) = ff(1,j,k,4)
         ff(-1,j,k,5) = ff(1,j,k,5)
         ff(-1,j,k,6) = ff(1,j,k,6)
         ff(-1,j,k,7) = ff(1,j,k,10)
         ff(-1,j,k,8) = ff(1,j,k,9)
         ff(-1,j,k,9) = ff(1,j,k,8)
         ff(-1,j,k,10) = ff(1,j,k,7)
         ff(-1,j,k,11) = ff(1,j,k,14)
         ff(-1,j,k,12) = ff(1,j,k,13)
         ff(-1,j,k,13) = ff(1,j,k,12)
         ff(-1,j,k,14) = ff(1,j,k,11)
         ff(-1,j,k,15) = ff(1,j,k,15)
         ff(-1,j,k,16) = ff(1,j,k,16)
         ff(-1,j,k,17) = ff(1,j,k,17)
         ff(-1,j,k,18) = ff(1,j,k,18)
         ff(-1,j,k,19) = ff(1,j,k,26)
         ff(-1,j,k,20) = ff(1,j,k,25)
         ff(-1,j,k,21) = ff(1,j,k,24)
         ff(-1,j,k,22) = ff(1,j,k,23)
         ff(-1,j,k,23) = ff(1,j,k,22)
         ff(-1,j,k,24) = ff(1,j,k,21)
         ff(-1,j,k,25) = ff(1,j,k,20)
         ff(-1,j,k,26) = ff(1,j,k,19)

         ff(nx+1,j,k,0) = ff(nx-1,j,k,0)
         ff(nx+1,j,k,1) = ff(nx-1,j,k,2)
         ff(nx+1,j,k,2) = ff(nx-1,j,k,1)
         ff(nx+1,j,k,3) = ff(nx-1,j,k,3)
         ff(nx+1,j,k,4) = ff(nx-1,j,k,4)
         ff(nx+1,j,k,5) = ff(nx-1,j,k,5)
         ff(nx+1,j,k,6) = ff(nx-1,j,k,6)
         ff(nx+1,j,k,7) = ff(nx-1,j,k,10)
         ff(nx+1,j,k,8) = ff(nx-1,j,k,9)
         ff(nx+1,j,k,9) = ff(nx-1,j,k,8)
         ff(nx+1,j,k,10) = ff(nx-1,j,k,7)
         ff(nx+1,j,k,11) = ff(nx-1,j,k,14)
         ff(nx+1,j,k,12) = ff(nx-1,j,k,13)
         ff(nx+1,j,k,13) = ff(nx-1,j,k,12)
         ff(nx+1,j,k,14) = ff(nx-1,j,k,11)
         ff(nx+1,j,k,15) = ff(nx-1,j,k,15)
         ff(nx+1,j,k,16) = ff(nx-1,j,k,16)
         ff(nx+1,j,k,17) = ff(nx-1,j,k,17)
         ff(nx+1,j,k,18) = ff(nx-1,j,k,18)
         ff(nx+1,j,k,19) = ff(nx-1,j,k,26)
         ff(nx+1,j,k,20) = ff(nx-1,j,k,25)
         ff(nx+1,j,k,21) = ff(nx-1,j,k,24)
         ff(nx+1,j,k,22) = ff(nx-1,j,k,23)
         ff(nx+1,j,k,23) = ff(nx-1,j,k,22)
         ff(nx+1,j,k,24) = ff(nx-1,j,k,21)
         ff(nx+1,j,k,25) = ff(nx-1,j,k,20)
         ff(nx+1,j,k,26) = ff(nx-1,j,k,19)

         do ii = 0, 26

         ff(-1,j,k,ii) = ff(nx-1,j,k,ii)
         ff(nx+1,j,k,ii) = ff(1,j,k,ii)

         enddo

      enddo
      enddo

      do ii = 1, 26

      do k = 0, nz
      do j = 0, ny
      do i = 0, nx

         ft(i,j,k) = ff(i-ic(ii,1),j-ic(ii,2),k-ic(ii,3),ii)

      enddo
      enddo
      enddo

      do k = 0, nz
      do j = 0, ny
      do i = 0, nx

         ff(i,j,k,ii) = ft(i,j,k)

      enddo
      enddo
      enddo

      enddo

      return
      end







**************************************************************
      subroutine weight(wa,ic,io,dst,dstm)
**************************************************************

      implicit double precision (a-h, o-z)
      dimension wa(0:26), ic(0:26,3), io(26), dst(26)

C     WEIGHTS

      wa(0) = 8.d0/27.d0
      wa(1) = 2.d0/27.d0
      wa(2) = 2.d0/27.d0
      wa(3) = 2.d0/27.d0
      wa(4) = 2.d0/27.d0
      wa(5) = 2.d0/27.d0
      wa(6) = 2.d0/27.d0
      wa(7) = 1.d0/54.d0
      wa(8) = 1.d0/54.d0
      wa(9) = 1.d0/54.d0
      wa(10) = 1.d0/54.d0
      wa(11) = 1.d0/54.d0
      wa(12) = 1.d0/54.d0
      wa(13) = 1.d0/54.d0
      wa(14) = 1.d0/54.d0
      wa(15) = 1.d0/54.d0
      wa(16) = 1.d0/54.d0
      wa(17) = 1.d0/54.d0
      wa(18) = 1.d0/54.d0
      wa(19) = 1.d0/216.d0
      wa(20) = 1.d0/216.d0
      wa(21) = 1.d0/216.d0
      wa(22) = 1.d0/216.d0
      wa(23) = 1.d0/216.d0
      wa(24) = 1.d0/216.d0
      wa(25) = 1.d0/216.d0
      wa(26) = 1.d0/216.d0

C     PARTICLE VELOCITY

      ic(0,1) = 0
      ic(0,2) = 0
      ic(0,3) = 0

      ic(1,1) = 1
      ic(1,2) = 0
      ic(1,3) = 0

      ic(2,1) = -1
      ic(2,2) = 0
      ic(2,3) = 0

      ic(3,1) = 0
      ic(3,2) = 1
      ic(3,3) = 0

      ic(4,1) = 0
      ic(4,2) = -1
      ic(4,3) = 0

      ic(5,1) = 0
      ic(5,2) = 0
      ic(5,3) = 1

      ic(6,1) = 0
      ic(6,2) = 0
      ic(6,3) = -1

      ic(7,1) = 1
      ic(7,2) = 1
      ic(7,3) = 0

      ic(8,1) = -1
      ic(8,2) = -1
      ic(8,3) = 0

      ic(9,1) = 1
      ic(9,2) = -1
      ic(9,3) = 0

      ic(10,1) = -1
      ic(10,2) = 1
      ic(10,3) = 0

      ic(11,1) = 1
      ic(11,2) = 0
      ic(11,3) = 1

      ic(12,1) = -1
      ic(12,2) = 0
      ic(12,3) = -1

      ic(13,1) = 1
      ic(13,2) = 0
      ic(13,3) = -1

      ic(14,1) = -1
      ic(14,2) = 0
      ic(14,3) = 1

      ic(15,1) = 0
      ic(15,2) = 1
      ic(15,3) = 1

      ic(16,1) = 0
      ic(16,2) = -1
      ic(16,3) = -1

      ic(17,1) = 0
      ic(17,2) = 1
      ic(17,3) = -1

      ic(18,1) = 0
      ic(18,2) = -1
      ic(18,3) = 1

      ic(19,1) = 1
      ic(19,2) = 1
      ic(19,3) = 1

      ic(20,1) = -1
      ic(20,2) = -1
      ic(20,3) = -1

      ic(21,1) = 1
      ic(21,2) = -1
      ic(21,3) = 1

      ic(22,1) = -1
      ic(22,2) = 1
      ic(22,3) = -1

      ic(23,1) = 1
      ic(23,2) = 1
      ic(23,3) = -1

      ic(24,1) = -1
      ic(24,2) = -1
      ic(24,3) = 1

      ic(25,1) = 1
      ic(25,2) = -1
      ic(25,3) = -1

      ic(26,1) = -1
      ic(26,2) = 1
      ic(26,3) = 1

C     BOUNCE BACK DIRECTIONS

      io(1) = 2
      io(2) = 1
      io(3) = 4
      io(4) = 3
      io(5) = 6
      io(6) = 5
      io(7) = 8
      io(8) = 7
      io(9) = 10
      io(10) = 9
      io(11) = 12
      io(12) = 11
      io(13) = 14
      io(14) = 13
      io(15) = 16
      io(16) = 15
      io(17) = 18
      io(18) = 17
      io(19) = 20
      io(20) = 19
      io(21) = 22
      io(22) = 21
      io(23) = 24
      io(24) = 23
      io(25) = 26
      io(26) = 25

C     WEIGHT BASED ON DISTANCE

      dst0 = 0.d0

      do ii = 1, 26

         dst0 = dst0
     *         +1.d0/dsqrt(ic(ii,1)**2+ic(ii,2)**2+ic(ii,3)**2+0.d0)

      enddo

      dstm = 0.d0

      do ii = 1, 26

         dst(ii) = 1.d0/dst0
     *            /dsqrt(ic(ii,1)**2+ic(ii,2)**2+ic(ii,3)**2+0.d0)

         if (dstm.lt.dst(ii)) dstm = dst(ii)

      enddo

      return
      end
