

#include "macrodef.fh"

  function compute_rhs_bssn(ex, T,X, Y, Z,                                     &
               chi    ,   trK    ,                                             &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,     &
               Axx    ,   Axy    ,   Axz    ,   Ayy    ,   Ayz    ,   Azz,     &
               Gamx   ,  Gamy    ,  Gamz    ,                                  &
               Lap    ,  betax   ,  betay   ,  betaz   ,                       &
               dtSfx  ,  dtSfy   ,  dtSfz   ,                                  &
               chi_rhs,   trK_rhs,                                             &
               gxx_rhs,   gxy_rhs,   gxz_rhs,   gyy_rhs,   gyz_rhs,   gzz_rhs, &
               Axx_rhs,   Axy_rhs,   Axz_rhs,   Ayy_rhs,   Ayz_rhs,   Azz_rhs, &
               Gamx_rhs,  Gamy_rhs,  Gamz_rhs,                                 &
               Lap_rhs,  betax_rhs,  betay_rhs,  betaz_rhs,                    &
               dtSfx_rhs,  dtSfy_rhs,  dtSfz_rhs,                              &
               rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz,                           &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,                      &
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,                      &
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,                      &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,                                        &
               ham_Res, movx_Res, movy_Res, movz_Res,                          &
                        Gmx_Res, Gmy_Res, Gmz_Res,                             &
               Symmetry,Lev,eps,co)  result(gont)
! calculate constraint violation when co=0               
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3), Symmetry,Lev,co
  real*8, intent(in ):: T
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: chi,dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: trK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ), target :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamx,Gamy,Gamz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Lap, betax, betay, betaz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dtSfx,  dtSfy,  dtSfz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: chi_rhs,trK_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gxx_rhs,gxy_rhs,gxz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gyy_rhs,gyz_rhs,gzz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Axx_rhs,Axy_rhs,Axz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Ayy_rhs,Ayz_rhs,Azz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamx_rhs,Gamy_rhs,Gamz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Lap_rhs, betax_rhs, betay_rhs, betaz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: dtSfx_rhs,dtSfy_rhs,dtSfz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: rho,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Sxx,Sxy,Sxz,Syy,Syz,Szz
! when out, physical second kind of connection  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out), target :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out), target :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out), target :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out), target :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out), target :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out), target :: Gamzyy, Gamzyz, Gamzzz
! when out, physical Ricci tensor  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out), target :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
  real*8,intent(in) :: eps
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: ham_Res, movx_Res, movy_Res, movz_Res
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gmx_Res, Gmy_Res, Gmz_Res
!  gont = 0: success; gont = 1: something wrong
  integer::gont

!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)), target :: gxx,gyy,gzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: chix,chiy,chiz
  real*8, dimension(ex(1),ex(2),ex(3)), target :: gxxx,gxyx,gxzx,gyyx,gyzx,gzzx
  real*8, dimension(ex(1),ex(2),ex(3)), target :: gxxy,gxyy,gxzy,gyyy,gyzy,gzzy
  real*8, dimension(ex(1),ex(2),ex(3)), target :: gxxz,gxyz,gxzz,gyyz,gyzz,gzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Lapx,Lapy,Lapz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betaxx,betaxy,betaxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betayx,betayy,betayz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betazx,betazy,betazz
  real*8, dimension(ex(1),ex(2),ex(3)), target :: Gamxx,Gamxy,Gamxz
  real*8, dimension(ex(1),ex(2),ex(3)), target :: Gamyx,Gamyy,Gamyz
  real*8, dimension(ex(1),ex(2),ex(3)), target :: Gamzx,Gamzy,Gamzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Kx,Ky,Kz,div_beta,S
  real*8, dimension(ex(1),ex(2),ex(3)) :: f,fxx,fxy,fxz,fyy,fyz,fzz
  real*8, dimension(ex(1),ex(2),ex(3)), target :: Gamxa,Gamya,Gamza,alpn1,chin1
  real*8, dimension(ex(1),ex(2),ex(3)), target :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)), target :: gupyy,gupyz,gupzz

  real*8,dimension(3) ::SSS,AAS,ASA,SAA,ASS,SAS,SSA
  real*8            :: dX, dY, dZ, PI
  real*8, parameter :: ZEO = 0.d0,ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0
  real*8, parameter :: EIGHT = 8.D0, HALF = 0.5D0, THR = 3.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  double precision,parameter::FF = 0.75d0,eta=2.d0
  real*8, parameter :: F1o3 = 1.D0/3.D0, F2o3 = 2.D0/3.D0,F3o2=1.5d0, F1o6 = 1.D0/6.D0
  real*8, parameter :: F16=1.6d1,F8=8.d0

#if (GAUGE == 2 || GAUGE == 3 || GAUGE == 4 || GAUGE == 5)
  real*8, dimension(ex(1),ex(2),ex(3)) :: reta
#endif

#if (GAUGE == 6 || GAUGE == 7)
  integer :: BHN,i,j,k
  real*8, dimension(9) :: Porg
  real*8, dimension(3) :: Mass
  real*8 :: r1,r2,M,A,w1,w2,C1,C2
  real*8, dimension(ex(1),ex(2),ex(3)) :: reta

  call getpbh(BHN,Porg,Mass)
#endif

!!! sanity check
  dX = sum(chi)+sum(trK)+sum(dxx)+sum(gxy)+sum(gxz)+sum(dyy)+sum(gyz)+sum(dzz) &
      +sum(Axx)+sum(Axy)+sum(Axz)+sum(Ayy)+sum(Ayz)+sum(Azz)                   &
      +sum(Gamx)+sum(Gamy)+sum(Gamz)                                           &
      +sum(Lap)+sum(betax)+sum(betay)+sum(betaz)
  if(dX.ne.dX) then
     if(sum(chi).ne.sum(chi))write(*,*)"bssn.f90: find NaN in chi"
     if(sum(trK).ne.sum(trK))write(*,*)"bssn.f90: find NaN in trk"
     if(sum(dxx).ne.sum(dxx))write(*,*)"bssn.f90: find NaN in dxx"
     if(sum(gxy).ne.sum(gxy))write(*,*)"bssn.f90: find NaN in gxy"
     if(sum(gxz).ne.sum(gxz))write(*,*)"bssn.f90: find NaN in gxz"
     if(sum(dyy).ne.sum(dyy))write(*,*)"bssn.f90: find NaN in dyy"
     if(sum(gyz).ne.sum(gyz))write(*,*)"bssn.f90: find NaN in gyz"
     if(sum(dzz).ne.sum(dzz))write(*,*)"bssn.f90: find NaN in dzz"
     if(sum(Axx).ne.sum(Axx))write(*,*)"bssn.f90: find NaN in Axx"
     if(sum(Axy).ne.sum(Axy))write(*,*)"bssn.f90: find NaN in Axy"
     if(sum(Axz).ne.sum(Axz))write(*,*)"bssn.f90: find NaN in Axz"
     if(sum(Ayy).ne.sum(Ayy))write(*,*)"bssn.f90: find NaN in Ayy"
     if(sum(Ayz).ne.sum(Ayz))write(*,*)"bssn.f90: find NaN in Ayz"
     if(sum(Azz).ne.sum(Azz))write(*,*)"bssn.f90: find NaN in Azz"
     if(sum(Gamx).ne.sum(Gamx))write(*,*)"bssn.f90: find NaN in Gamx"
     if(sum(Gamy).ne.sum(Gamy))write(*,*)"bssn.f90: find NaN in Gamy"
     if(sum(Gamz).ne.sum(Gamz))write(*,*)"bssn.f90: find NaN in Gamz"
     if(sum(Lap).ne.sum(Lap))write(*,*)"bssn.f90: find NaN in Lap"
     if(sum(betax).ne.sum(betax))write(*,*)"bssn.f90: find NaN in betax"
     if(sum(betay).ne.sum(betay))write(*,*)"bssn.f90: find NaN in betay"
     if(sum(betaz).ne.sum(betaz))write(*,*)"bssn.f90: find NaN in betaz"
     gont = 1
     return
  endif

  PI = dacos(-ONE)

  dX = X(2) - X(1)
  dY = Y(2) - Y(1)
  dZ = Z(2) - Z(1)

  alpn1 = Lap + ONE
  chin1 = chi + ONE
  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE

  call fderivs(ex,betax,betaxx,betaxy,betaxz,X,Y,Z,ANTI, SYM, SYM,Symmetry,Lev)
  call fderivs(ex,betay,betayx,betayy,betayz,X,Y,Z, SYM,ANTI, SYM,Symmetry,Lev)
  call fderivs(ex,betaz,betazx,betazy,betazz,X,Y,Z, SYM, SYM,ANTI,Symmetry,Lev)
  call fderivs(ex,chi,chix,chiy,chiz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)
  call fderivs(ex,dxx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,gxy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,Lev)
  call fderivs(ex,gxz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,Lev)
  call fderivs(ex,dyy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,gyz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,Lev)
  call fderivs(ex,dzz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)

  div_beta = betaxx + betayy + betazz

  chi_rhs = F2o3 *chin1*( alpn1 * trK - div_beta ) !rhs for chi

  gxx_rhs = - TWO * alpn1 * Axx    -  F2o3 * gxx * div_beta          + &
              TWO *(  gxx * betaxx +   gxy * betayx +   gxz * betazx)

  gyy_rhs = - TWO * alpn1 * Ayy    -  F2o3 * gyy * div_beta          + &
              TWO *(  gxy * betaxy +   gyy * betayy +   gyz * betazy)

  gzz_rhs = - TWO * alpn1 * Azz    -  F2o3 * gzz * div_beta          + &
              TWO *(  gxz * betaxz +   gyz * betayz +   gzz * betazz)

  gxy_rhs = - TWO * alpn1 * Axy    +  F1o3 * gxy    * div_beta       + &
                      gxx * betaxy                  +   gxz * betazy + &
                                       gyy * betayx +   gyz * betazx   &
                                                    -   gxy * betazz

  gyz_rhs = - TWO * alpn1 * Ayz    +  F1o3 * gyz    * div_beta       + &
                      gxy * betaxz +   gyy * betayz                  + &
                      gxz * betaxy                  +   gzz * betazy   &
                                                    -   gyz * betaxx
 
  gxz_rhs = - TWO * alpn1 * Axz    +  F1o3 * gxz    * div_beta       + &
                      gxx * betaxz +   gxy * betayz                  + &
                                       gyz * betayx +   gzz * betazx   &
                                                    -   gxz * betayy     !rhs for gij

! invert tilted metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz

  if(co == 0)then
! Gam^i_Res = Gam^i + gup^ij_,j
  Gmx_Res = Gamx - (gupxx*(gupxx*gxxx+gupxy*gxyx+gupxz*gxzx)&
                   +gupxy*(gupxx*gxyx+gupxy*gyyx+gupxz*gyzx)&
                   +gupxz*(gupxx*gxzx+gupxy*gyzx+gupxz*gzzx)&
                   +gupxx*(gupxy*gxxy+gupyy*gxyy+gupyz*gxzy)&
                   +gupxy*(gupxy*gxyy+gupyy*gyyy+gupyz*gyzy)&
                   +gupxz*(gupxy*gxzy+gupyy*gyzy+gupyz*gzzy)&
                   +gupxx*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)&
                   +gupxy*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)&
                   +gupxz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz))
  Gmy_Res = Gamy - (gupxx*(gupxy*gxxx+gupyy*gxyx+gupyz*gxzx)&
                   +gupxy*(gupxy*gxyx+gupyy*gyyx+gupyz*gyzx)&
                   +gupxz*(gupxy*gxzx+gupyy*gyzx+gupyz*gzzx)&
                   +gupxy*(gupxy*gxxy+gupyy*gxyy+gupyz*gxzy)&
                   +gupyy*(gupxy*gxyy+gupyy*gyyy+gupyz*gyzy)&
                   +gupyz*(gupxy*gxzy+gupyy*gyzy+gupyz*gzzy)&
                   +gupxy*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)&
                   +gupyy*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)&
                   +gupyz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz))
  Gmz_Res = Gamz - (gupxx*(gupxz*gxxx+gupyz*gxyx+gupzz*gxzx)&
                   +gupxy*(gupxz*gxyx+gupyz*gyyx+gupzz*gyzx)&
                   +gupxz*(gupxz*gxzx+gupyz*gyzx+gupzz*gzzx)&
                   +gupxy*(gupxz*gxxy+gupyz*gxyy+gupzz*gxzy)&
                   +gupyy*(gupxz*gxyy+gupyz*gyyy+gupzz*gyzy)&
                   +gupyz*(gupxz*gxzy+gupyz*gyzy+gupzz*gzzy)&
                   +gupxz*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)&
                   +gupyz*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)&
                   +gupzz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz))
  endif

! second kind of connection
  Gamxxx =HALF*( gupxx*gxxx + gupxy*(TWO*gxyx - gxxy ) + gupxz*(TWO*gxzx - gxxz ))
  Gamyxx =HALF*( gupxy*gxxx + gupyy*(TWO*gxyx - gxxy ) + gupyz*(TWO*gxzx - gxxz ))
  Gamzxx =HALF*( gupxz*gxxx + gupyz*(TWO*gxyx - gxxy ) + gupzz*(TWO*gxzx - gxxz ))
 
  Gamxyy =HALF*( gupxx*(TWO*gxyy - gyyx ) + gupxy*gyyy + gupxz*(TWO*gyzy - gyyz ))
  Gamyyy =HALF*( gupxy*(TWO*gxyy - gyyx ) + gupyy*gyyy + gupyz*(TWO*gyzy - gyyz ))
  Gamzyy =HALF*( gupxz*(TWO*gxyy - gyyx ) + gupyz*gyyy + gupzz*(TWO*gyzy - gyyz ))

  Gamxzz =HALF*( gupxx*(TWO*gxzz - gzzx ) + gupxy*(TWO*gyzz - gzzy ) + gupxz*gzzz)
  Gamyzz =HALF*( gupxy*(TWO*gxzz - gzzx ) + gupyy*(TWO*gyzz - gzzy ) + gupyz*gzzz)
  Gamzzz =HALF*( gupxz*(TWO*gxzz - gzzx ) + gupyz*(TWO*gyzz - gzzy ) + gupzz*gzzz)

  Gamxxy =HALF*( gupxx*gxxy + gupxy*gyyx + gupxz*( gxzy + gyzx - gxyz ) )
  Gamyxy =HALF*( gupxy*gxxy + gupyy*gyyx + gupyz*( gxzy + gyzx - gxyz ) )
  Gamzxy =HALF*( gupxz*gxxy + gupyz*gyyx + gupzz*( gxzy + gyzx - gxyz ) )

  Gamxxz =HALF*( gupxx*gxxz + gupxy*( gxyz + gyzx - gxzy ) + gupxz*gzzx )
  Gamyxz =HALF*( gupxy*gxxz + gupyy*( gxyz + gyzx - gxzy ) + gupyz*gzzx )
  Gamzxz =HALF*( gupxz*gxxz + gupyz*( gxyz + gyzx - gxzy ) + gupzz*gzzx )

  Gamxyz =HALF*( gupxx*( gxyz + gxzy - gyzx ) + gupxy*gyyz + gupxz*gzzy )
  Gamyyz =HALF*( gupxy*( gxyz + gxzy - gyzx ) + gupyy*gyyz + gupyz*gzzy )
  Gamzyz =HALF*( gupxz*( gxyz + gxzy - gyzx ) + gupyz*gyyz + gupzz*gzzy )
! Raise indices of \tilde A_{ij} and store in R_ij

  Rxx =    gupxx * gupxx * Axx + gupxy * gupxy * Ayy + gupxz * gupxz * Azz + &
      TWO*(gupxx * gupxy * Axy + gupxx * gupxz * Axz + gupxy * gupxz * Ayz)

  Ryy =    gupxy * gupxy * Axx + gupyy * gupyy * Ayy + gupyz * gupyz * Azz + &
      TWO*(gupxy * gupyy * Axy + gupxy * gupyz * Axz + gupyy * gupyz * Ayz)

  Rzz =    gupxz * gupxz * Axx + gupyz * gupyz * Ayy + gupzz * gupzz * Azz + &
      TWO*(gupxz * gupyz * Axy + gupxz * gupzz * Axz + gupyz * gupzz * Ayz)

  Rxy =    gupxx * gupxy * Axx + gupxy * gupyy * Ayy + gupxz * gupyz * Azz + &
          (gupxx * gupyy       + gupxy * gupxy)* Axy                       + &
          (gupxx * gupyz       + gupxz * gupxy)* Axz                       + &
          (gupxy * gupyz       + gupxz * gupyy)* Ayz

  Rxz =    gupxx * gupxz * Axx + gupxy * gupyz * Ayy + gupxz * gupzz * Azz + &
          (gupxx * gupyz       + gupxy * gupxz)* Axy                       + &
          (gupxx * gupzz       + gupxz * gupxz)* Axz                       + &
          (gupxy * gupzz       + gupxz * gupyz)* Ayz

  Ryz =    gupxy * gupxz * Axx + gupyy * gupyz * Ayy + gupyz * gupzz * Azz + &
          (gupxy * gupyz       + gupyy * gupxz)* Axy                       + &
          (gupxy * gupzz       + gupyz * gupxz)* Axz                       + &
          (gupyy * gupzz       + gupyz * gupyz)* Ayz

! Right hand side for Gam^i without shift terms...
  call fderivs(ex,Lap,Lapx,Lapy,Lapz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  call fderivs(ex,trK,Kx,Ky,Kz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)

   Gamx_rhs = - TWO * (   Lapx * Rxx +   Lapy * Rxy +   Lapz * Rxz ) + &
        TWO * alpn1 * (                                                &
        -F3o2/chin1 * (   chix * Rxx +   chiy * Rxy +   chiz * Rxz ) - &
              gupxx * (   F2o3 * Kx  +  EIGHT * PI * Sx            ) - &
              gupxy * (   F2o3 * Ky  +  EIGHT * PI * Sy            ) - &
              gupxz * (   F2o3 * Kz  +  EIGHT * PI * Sz            ) + &
                        Gamxxx * Rxx + Gamxyy * Ryy + Gamxzz * Rzz   + &
                TWO * ( Gamxxy * Rxy + Gamxxz * Rxz + Gamxyz * Ryz ) )

   Gamy_rhs = - TWO * (   Lapx * Rxy +   Lapy * Ryy +   Lapz * Ryz ) + &
        TWO * alpn1 * (                                                &
        -F3o2/chin1 * (   chix * Rxy +  chiy * Ryy +    chiz * Ryz ) - &
              gupxy * (   F2o3 * Kx  +  EIGHT * PI * Sx            ) - &
              gupyy * (   F2o3 * Ky  +  EIGHT * PI * Sy            ) - &
              gupyz * (   F2o3 * Kz  +  EIGHT * PI * Sz            ) + &
                        Gamyxx * Rxx + Gamyyy * Ryy + Gamyzz * Rzz   + &
                TWO * ( Gamyxy * Rxy + Gamyxz * Rxz + Gamyyz * Ryz ) )

   Gamz_rhs = - TWO * (   Lapx * Rxz +   Lapy * Ryz +   Lapz * Rzz ) + &
        TWO * alpn1 * (                                                &
        -F3o2/chin1 * (   chix * Rxz +  chiy * Ryz +    chiz * Rzz ) - &
              gupxz * (   F2o3 * Kx  +  EIGHT * PI * Sx            ) - &
              gupyz * (   F2o3 * Ky  +  EIGHT * PI * Sy            ) - &
              gupzz * (   F2o3 * Kz  +  EIGHT * PI * Sz            ) + &
                        Gamzxx * Rxx + Gamzyy * Ryy + Gamzzz * Rzz   + &
                TWO * ( Gamzxy * Rxy + Gamzxz * Rxz + Gamzyz * Ryz ) )

  call fdderivs(ex,betax,gxxx,gxyx,gxzx,gyyx,gyzx,gzzx,&
                X,Y,Z,ANTI,SYM, SYM ,Symmetry,Lev)
  call fdderivs(ex,betay,gxxy,gxyy,gxzy,gyyy,gyzy,gzzy,&
                X,Y,Z,SYM ,ANTI,SYM ,Symmetry,Lev)
  call fdderivs(ex,betaz,gxxz,gxyz,gxzz,gyyz,gyzz,gzzz,&
                X,Y,Z,SYM ,SYM, ANTI,Symmetry,Lev)

  fxx = gxxx + gxyy + gxzz
  fxy = gxyx + gyyy + gyzz
  fxz = gxzx + gyzy + gzzz

  Gamxa =       gupxx * Gamxxx + gupyy * Gamxyy + gupzz * Gamxzz + &
          TWO*( gupxy * Gamxxy + gupxz * Gamxxz + gupyz * Gamxyz )
  Gamya =       gupxx * Gamyxx + gupyy * Gamyyy + gupzz * Gamyzz + &
          TWO*( gupxy * Gamyxy + gupxz * Gamyxz + gupyz * Gamyyz )
  Gamza =       gupxx * Gamzxx + gupyy * Gamzyy + gupzz * Gamzzz + &
          TWO*( gupxy * Gamzxy + gupxz * Gamzxz + gupyz * Gamzyz )

  call fderivs(ex,Gamx,Gamxx,Gamxy,Gamxz,X,Y,Z,ANTI,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,Gamy,Gamyx,Gamyy,Gamyz,X,Y,Z,SYM ,ANTI,SYM ,Symmetry,Lev)
  call fderivs(ex,Gamz,Gamzx,Gamzy,Gamzz,X,Y,Z,SYM ,SYM ,ANTI,Symmetry,Lev)

  Gamx_rhs =               Gamx_rhs +  F2o3 *  Gamxa * div_beta        - &
                     Gamxa * betaxx - Gamya * betaxy - Gamza * betaxz  + &
             F1o3 * (gupxx * fxx    + gupxy * fxy    + gupxz * fxz    ) + &
                     gupxx * gxxx   + gupyy * gyyx   + gupzz * gzzx    + &
              TWO * (gupxy * gxyx   + gupxz * gxzx   + gupyz * gyzx  )

  Gamy_rhs =               Gamy_rhs +  F2o3 *  Gamya * div_beta        - &
                     Gamxa * betayx - Gamya * betayy - Gamza * betayz  + &
             F1o3 * (gupxy * fxx    + gupyy * fxy    + gupyz * fxz    ) + &
                     gupxx * gxxy   + gupyy * gyyy   + gupzz * gzzy    + &
              TWO * (gupxy * gxyy   + gupxz * gxzy   + gupyz * gyzy  )

  Gamz_rhs =               Gamz_rhs +  F2o3 *  Gamza * div_beta        - &
                     Gamxa * betazx - Gamya * betazy - Gamza * betazz  + &
             F1o3 * (gupxz * fxx    + gupyz * fxy    + gupzz * fxz    ) + &
                     gupxx * gxxz   + gupyy * gyyz   + gupzz * gzzz    + &
              TWO * (gupxy * gxyz   + gupxz * gxzz   + gupyz * gyzz  )    !rhs for Gam^i

!first kind of connection stored in gij,k
  gxxx = gxx * Gamxxx + gxy * Gamyxx + gxz * Gamzxx
  gxyx = gxx * Gamxxy + gxy * Gamyxy + gxz * Gamzxy
  gxzx = gxx * Gamxxz + gxy * Gamyxz + gxz * Gamzxz
  gyyx = gxx * Gamxyy + gxy * Gamyyy + gxz * Gamzyy
  gyzx = gxx * Gamxyz + gxy * Gamyyz + gxz * Gamzyz
  gzzx = gxx * Gamxzz + gxy * Gamyzz + gxz * Gamzzz

  gxxy = gxy * Gamxxx + gyy * Gamyxx + gyz * Gamzxx
  gxyy = gxy * Gamxxy + gyy * Gamyxy + gyz * Gamzxy
  gxzy = gxy * Gamxxz + gyy * Gamyxz + gyz * Gamzxz
  gyyy = gxy * Gamxyy + gyy * Gamyyy + gyz * Gamzyy
  gyzy = gxy * Gamxyz + gyy * Gamyyz + gyz * Gamzyz
  gzzy = gxy * Gamxzz + gyy * Gamyzz + gyz * Gamzzz

  gxxz = gxz * Gamxxx + gyz * Gamyxx + gzz * Gamzxx
  gxyz = gxz * Gamxxy + gyz * Gamyxy + gzz * Gamzxy
  gxzz = gxz * Gamxxz + gyz * Gamyxz + gzz * Gamzxz
  gyyz = gxz * Gamxyy + gyz * Gamyyy + gzz * Gamzyy
  gyzz = gxz * Gamxyz + gyz * Gamyyz + gzz * Gamzyz
  gzzz = gxz * Gamxzz + gyz * Gamyzz + gzz * Gamzzz

!compute Ricci tensor for tilted metric
   call fdderivs(ex,dxx,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
   Rxx =  - HALF * ( gupxx * fxx + gupyy * fyy + gupzz * fzz) - &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) 

   call fdderivs(ex,dyy,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
   Ryy =  - HALF * ( gupxx * fxx + gupyy * fyy + gupzz * fzz) - &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) 

   call fdderivs(ex,dzz,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
   Rzz =  - HALF * ( gupxx * fxx + gupyy * fyy + gupzz * fzz) - &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz )

   call fdderivs(ex,gxy,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,ANTI, ANTI,SYM ,symmetry,Lev)
   Rxy =  - HALF * ( gupxx * fxx + gupyy * fyy + gupzz * fzz) - &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz )

   call fdderivs(ex,gxz,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,ANTI ,SYM ,ANTI,symmetry,Lev)
   Rxz =  - HALF * ( gupxx * fxx + gupyy * fyy + gupzz * fzz) - &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz )

   call fdderivs(ex,gyz,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,ANTI ,ANTI,symmetry,Lev)
   Ryz =  - HALF * ( gupxx * fxx + gupyy * fyy + gupzz * fzz) - &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz )

  block
    real*8, pointer :: Rxx1d(:), Ryy1d(:), Rzz1d(:), Rxy1d(:), Rxz1d(:), Ryz1d(:)
    real*8, pointer :: Gamxa1d(:), Gamya1d(:), Gamza1d(:)
    real*8, pointer :: Gamxx1d(:), Gamxy1d(:), Gamxz1d(:)
    real*8, pointer :: Gamyx1d(:), Gamyy1d(:), Gamyz1d(:)
    real*8, pointer :: Gamzx1d(:), Gamzy1d(:), Gamzz1d(:)
    real*8, pointer :: Gamxxx1d(:), Gamxxy1d(:), Gamxxz1d(:)
    real*8, pointer :: Gamxyy1d(:), Gamxyz1d(:), Gamxzz1d(:)
    real*8, pointer :: Gamyxx1d(:), Gamyxy1d(:), Gamyxz1d(:)
    real*8, pointer :: Gamyyy1d(:), Gamyyz1d(:), Gamyzz1d(:)
    real*8, pointer :: Gamzxx1d(:), Gamzxy1d(:), Gamzxz1d(:)
    real*8, pointer :: Gamzyy1d(:), Gamzyz1d(:), Gamzzz1d(:)
    real*8, pointer :: gxx1d(:),gyy1d(:),gzz1d(:)
    real*8, pointer :: gxy1d(:),gxz1d(:),gyz1d(:)
    real*8, pointer :: gxxx1d(:),gxyx1d(:),gxzx1d(:),gyyx1d(:),gyzx1d(:),gzzx1d(:)
    real*8, pointer :: gxxy1d(:),gxyy1d(:),gxzy1d(:),gyyy1d(:),gyzy1d(:),gzzy1d(:)
    real*8, pointer :: gxxz1d(:),gxyz1d(:),gxzz1d(:),gyyz1d(:),gyzz1d(:),gzzz1d(:)
    real*8, pointer :: gupxx1d(:),gupxy1d(:),gupxz1d(:),gupyy1d(:),gupyz1d(:),gupzz1d(:)
    integer :: N, i
    N = ex(1)*ex(2)*ex(3)
    Rxx1d(1: N) => Rxx
    Ryy1d(1: N) => Ryy
    Rzz1d(1: N) => Rzz
    Rxy1d(1: N) => Rxy
    Rxz1d(1: N) => Rxz
    Ryz1d(1: N) => Ryz

    Gamxa1d(1: N) => Gamxa
    Gamya1d(1: N) => Gamya
    Gamza1d(1: N) => Gamza

    Gamxx1d(1: N) => Gamxx
    Gamxy1d(1: N) => Gamxy
    Gamxz1d(1: N) => Gamxz
    Gamyx1d(1: N) => Gamyx
    Gamyy1d(1: N) => Gamyy
    Gamyz1d(1: N) => Gamyz
    Gamzx1d(1: N) => Gamzx
    Gamzy1d(1: N) => Gamzy
    Gamzz1d(1: N) => Gamzz

    Gamxxx1d(1: N) => Gamxxx
    Gamxxy1d(1: N) => Gamxxy
    Gamxxz1d(1: N) => Gamxxz
    Gamxyy1d(1: N) => Gamxyy
    Gamxyz1d(1: N) => Gamxyz
    Gamxzz1d(1: N) => Gamxzz
    Gamyxx1d(1: N) => Gamyxx
    Gamyxy1d(1: N) => Gamyxy
    Gamyxz1d(1: N) => Gamyxz
    Gamyyy1d(1: N) => Gamyyy
    Gamyyz1d(1: N) => Gamyyz
    Gamyzz1d(1: N) => Gamyzz
    Gamzxx1d(1: N) => Gamzxx
    Gamzxy1d(1: N) => Gamzxy
    Gamzxz1d(1: N) => Gamzxz
    Gamzyy1d(1: N) => Gamzyy
    Gamzyz1d(1: N) => Gamzyz
    Gamzzz1d(1: N) => Gamzzz

    gxx1d(1: N) => gxx
    gyy1d(1: N) => gyy
    gzz1d(1: N) => gzz

    gxy1d(1: N) => gxy
    gxz1d(1: N) => gxz
    gyz1d(1: N) => gyz

    gxxx1d(1: N) => gxxx
    gxyx1d(1: N) => gxyx
    gxzx1d(1: N) => gxzx
    gyyx1d(1: N) => gyyx
    gyzx1d(1: N) => gyzx
    gzzx1d(1: N) => gzzx
    gxxy1d(1: N) => gxxy
    gxyy1d(1: N) => gxyy
    gxzy1d(1: N) => gxzy
    gyyy1d(1: N) => gyyy
    gyzy1d(1: N) => gyzy
    gzzy1d(1: N) => gzzy
    gxxz1d(1: N) => gxxz
    gxyz1d(1: N) => gxyz
    gxzz1d(1: N) => gxzz
    gyyz1d(1: N) => gyyz
    gyzz1d(1: N) => gyzz
    gzzz1d(1: N) => gzzz

    gupxx1d(1: N) => gupxx
    gupxy1d(1: N) => gupxy
    gupxz1d(1: N) => gupxz
    gupyy1d(1: N) => gupyy
    gupyz1d(1: N) => gupyz
    gupzz1d(1: N) => gupzz

    do i = 1, N
      Rxx1d(i) = Rxx1d(i) + &
               gxx1d(i) * Gamxx1d(i) + gxy1d(i) * Gamyx1d(i) + gxz1d(i) * Gamzx1d(i) + &
               gxxx1d(i) * Gamxa1d(i) + gxyx1d(i) * Gamya1d(i) + gxzx1d(i) * Gamza1d(i) + &
               gupxx1d(i) * ( &
                   TWO*(Gamxxx1d(i) * gxxx1d(i) + Gamyxx1d(i) * gxyx1d(i) + Gamzxx1d(i) * gxzx1d(i)) + &
                        Gamxxx1d(i) * gxxx1d(i) + Gamyxx1d(i) * gxxy1d(i) + Gamzxx1d(i) * gxxz1d(i) )+ &
               gupxy1d(i) * ( &
                   TWO*(Gamxxx1d(i) * gxyx1d(i) + Gamyxx1d(i) * gyyx1d(i) + Gamzxx1d(i) * gyzx1d(i) + &
                        Gamxxy1d(i) * gxxx1d(i) + Gamyxy1d(i) * gxyx1d(i) + Gamzxy1d(i) * gxzx1d(i)) + &
                        Gamxxy1d(i) * gxxx1d(i) + Gamyxy1d(i) * gxxy1d(i) + Gamzxy1d(i) * gxxz1d(i) + &
                        Gamxxx1d(i) * gxyx1d(i) + Gamyxx1d(i) * gxyy1d(i) + Gamzxx1d(i) * gxyz1d(i) )+ &
               gupxz1d(i) * ( &
                   TWO*(Gamxxx1d(i) * gxzx1d(i) + Gamyxx1d(i) * gyzx1d(i) + Gamzxx1d(i) * gzzx1d(i) + &
                        Gamxxz1d(i) * gxxx1d(i) + Gamyxz1d(i) * gxyx1d(i) + Gamzxz1d(i) * gxzx1d(i)) + &
                        Gamxxz1d(i) * gxxx1d(i) + Gamyxz1d(i) * gxxy1d(i) + Gamzxz1d(i) * gxxz1d(i) + &
                        Gamxxx1d(i) * gxzx1d(i) + Gamyxx1d(i) * gxzy1d(i) + Gamzxx1d(i) * gxzz1d(i) )+ &
               gupyy1d(i) * ( &
                   TWO*(Gamxxy1d(i) * gxyx1d(i) + Gamyxy1d(i) * gyyx1d(i) + Gamzxy1d(i) * gyzx1d(i)) + &
                        Gamxxy1d(i) * gxyx1d(i) + Gamyxy1d(i) * gxyy1d(i) + Gamzxy1d(i) * gxyz1d(i) )+ &
               gupyz1d(i) * ( &
                   TWO*(Gamxxy1d(i) * gxzx1d(i) + Gamyxy1d(i) * gyzx1d(i) + Gamzxy1d(i) * gzzx1d(i) + &
                        Gamxxz1d(i) * gxyx1d(i) + Gamyxz1d(i) * gyyx1d(i) + Gamzxz1d(i) * gyzx1d(i)) + &
                        Gamxxz1d(i) * gxyx1d(i) + Gamyxz1d(i) * gxyy1d(i) + Gamzxz1d(i) * gxyz1d(i) + &
                        Gamxxy1d(i) * gxzx1d(i) + Gamyxy1d(i) * gxzy1d(i) + Gamzxy1d(i) * gxzz1d(i) )+ &
               gupzz1d(i) * ( &
                   TWO*(Gamxxz1d(i) * gxzx1d(i) + Gamyxz1d(i) * gyzx1d(i) + Gamzxz1d(i) * gzzx1d(i)) + &
                        Gamxxz1d(i) * gxzx1d(i) + Gamyxz1d(i) * gxzy1d(i) + Gamzxz1d(i) * gxzz1d(i) )
        
      Ryy1d(i) = Ryy1d(i) + &
               gxy1d(i) * Gamxy1d(i) + gyy1d(i) * Gamyy1d(i) + gyz1d(i) * Gamzy1d(i) + &
               gxyy1d(i) * Gamxa1d(i) + gyyy1d(i) * Gamya1d(i) + gyzy1d(i) * Gamza1d(i) + &
               gupxx1d(i) * ( &
                   TWO*(Gamxxy1d(i) * gxxy1d(i) + Gamyxy1d(i) * gxyy1d(i) + Gamzxy1d(i) * gxzy1d(i)) + &
                        Gamxxy1d(i) * gxyx1d(i) + Gamyxy1d(i) * gxyy1d(i) + Gamzxy1d(i) * gxyz1d(i) )+ &
               gupxy1d(i) * ( &
                   TWO*(Gamxxy1d(i) * gxyy1d(i) + Gamyxy1d(i) * gyyy1d(i) + Gamzxy1d(i) * gyzy1d(i) + &
                        Gamxyy1d(i) * gxxy1d(i) + Gamyyy1d(i) * gxyy1d(i) + Gamzyy1d(i) * gxzy1d(i)) + &
                        Gamxyy1d(i) * gxyx1d(i) + Gamyyy1d(i) * gxyy1d(i) + Gamzyy1d(i) * gxyz1d(i) + &
                        Gamxxy1d(i) * gyyx1d(i) + Gamyxy1d(i) * gyyy1d(i) + Gamzxy1d(i) * gyyz1d(i) )+ &
               gupxz1d(i) * ( &
                   TWO*(Gamxxy1d(i) * gxzy1d(i) + Gamyxy1d(i) * gyzy1d(i) + Gamzxy1d(i) * gzzy1d(i) + &
                        Gamxyz1d(i) * gxxy1d(i) + Gamyyz1d(i) * gxyy1d(i) + Gamzyz1d(i) * gxzy1d(i)) + &
                        Gamxyz1d(i) * gxyx1d(i) + Gamyyz1d(i) * gxyy1d(i) + Gamzyz1d(i) * gxyz1d(i) + &
                        Gamxxy1d(i) * gyzx1d(i) + Gamyxy1d(i) * gyzy1d(i) + Gamzxy1d(i) * gyzz1d(i) )+ &
               gupyy1d(i) * ( &
                   TWO*(Gamxyy1d(i) * gxyy1d(i) + Gamyyy1d(i) * gyyy1d(i) + Gamzyy1d(i) * gyzy1d(i)) + &
                        Gamxyy1d(i) * gyyx1d(i) + Gamyyy1d(i) * gyyy1d(i) + Gamzyy1d(i) * gyyz1d(i) )+ &
               gupyz1d(i) * ( &
                   TWO*(Gamxyy1d(i) * gxzy1d(i) + Gamyyy1d(i) * gyzy1d(i) + Gamzyy1d(i) * gzzy1d(i) + &
                        Gamxyz1d(i) * gxyy1d(i) + Gamyyz1d(i) * gyyy1d(i) + Gamzyz1d(i) * gyzy1d(i)) + &
                        Gamxyz1d(i) * gyyx1d(i) + Gamyyz1d(i) * gyyy1d(i) + Gamzyz1d(i) * gyyz1d(i) + &
                        Gamxyy1d(i) * gyzx1d(i) + Gamyyy1d(i) * gyzy1d(i) + Gamzyy1d(i) * gyzz1d(i) )+ &
               gupzz1d(i) * ( &
                   TWO*(Gamxyz1d(i) * gxzy1d(i) + Gamyyz1d(i) * gyzy1d(i) + Gamzyz1d(i) * gzzy1d(i)) + &
                        Gamxyz1d(i) * gyzx1d(i) + Gamyyz1d(i) * gyzy1d(i) + Gamzyz1d(i) * gyzz1d(i) )
      
      Rzz1d(i) = Rzz1d(i) + &
               gxz1d(i) * Gamxz1d(i) + gyz1d(i) * Gamyz1d(i) + gzz1d(i) * Gamzz1d(i) + &
               gxzz1d(i) * Gamxa1d(i) + gyzz1d(i) * Gamya1d(i) + gzzz1d(i) * Gamza1d(i) + &
               gupxx1d(i) * ( &
                   TWO*(Gamxxz1d(i) * gxxz1d(i) + Gamyxz1d(i) * gxyz1d(i) + Gamzxz1d(i) * gxzz1d(i)) + &
                        Gamxxz1d(i) * gxzx1d(i) + Gamyxz1d(i) * gxzy1d(i) + Gamzxz1d(i) * gxzz1d(i) )+ &
               gupxy1d(i) * ( &
                   TWO*(Gamxxz1d(i) * gxyz1d(i) + Gamyxz1d(i) * gyyz1d(i) + Gamzxz1d(i) * gyzz1d(i) + &
                        Gamxyz1d(i) * gxxz1d(i) + Gamyyz1d(i) * gxyz1d(i) + Gamzyz1d(i) * gxzz1d(i)) + &
                        Gamxyz1d(i) * gxzx1d(i) + Gamyyz1d(i) * gxzy1d(i) + Gamzyz1d(i) * gxzz1d(i) + &
                        Gamxxz1d(i) * gyzx1d(i) + Gamyxz1d(i) * gyzy1d(i) + Gamzxz1d(i) * gyzz1d(i) )+ &
               gupxz1d(i) * ( &
                   TWO*(Gamxxz1d(i) * gxzz1d(i) + Gamyxz1d(i) * gyzz1d(i) + Gamzxz1d(i) * gzzz1d(i) + &
                        Gamxzz1d(i) * gxxz1d(i) + Gamyzz1d(i) * gxyz1d(i) + Gamzzz1d(i) * gxzz1d(i)) + &
                        Gamxzz1d(i) * gxzx1d(i) + Gamyzz1d(i) * gxzy1d(i) + Gamzzz1d(i) * gxzz1d(i) + &
                        Gamxxz1d(i) * gzzx1d(i) + Gamyxz1d(i) * gzzy1d(i) + Gamzxz1d(i) * gzzz1d(i) )+ &
               gupyy1d(i) * ( &
                   TWO*(Gamxyz1d(i) * gxyz1d(i) + Gamyyz1d(i) * gyyz1d(i) + Gamzyz1d(i) * gyzz1d(i)) + &
                        Gamxyz1d(i) * gyzx1d(i) + Gamyyz1d(i) * gyzy1d(i) + Gamzyz1d(i) * gyzz1d(i) )+ &
               gupyz1d(i) * ( &
                   TWO*(Gamxyz1d(i) * gxzz1d(i) + Gamyyz1d(i) * gyzz1d(i) + Gamzyz1d(i) * gzzz1d(i) + &
                        Gamxzz1d(i) * gxyz1d(i) + Gamyzz1d(i) * gyyz1d(i) + Gamzzz1d(i) * gyzz1d(i)) + &
                        Gamxzz1d(i) * gyzx1d(i) + Gamyzz1d(i) * gyzy1d(i) + Gamzzz1d(i) * gyzz1d(i) + &
                        Gamxyz1d(i) * gzzx1d(i) + Gamyyz1d(i) * gzzy1d(i) + Gamzyz1d(i) * gzzz1d(i) )+ &
               gupzz1d(i) * ( &
                   TWO*(Gamxzz1d(i) * gxzz1d(i) + Gamyzz1d(i) * gyzz1d(i) + Gamzzz1d(i) * gzzz1d(i)) + &
                        Gamxzz1d(i) * gzzx1d(i) + Gamyzz1d(i) * gzzy1d(i) + Gamzzz1d(i) * gzzz1d(i) )
      Rxy1d(i) = Rxy1d(i) + &
        HALF*( gxx1d(i) * Gamxy1d(i) + gxy1d(i) * Gamyy1d(i) + gxz1d(i) * Gamzy1d(i) + &
               gxy1d(i) * Gamxx1d(i) + gyy1d(i) * Gamyx1d(i) + gyz1d(i) * Gamzx1d(i) + &
               gxyx1d(i) * Gamxa1d(i) + gyyx1d(i) * Gamya1d(i) + gyzx1d(i) * Gamza1d(i) + &
               gxxy1d(i) * Gamxa1d(i) + gxyy1d(i) * Gamya1d(i) + gxzy1d(i) * Gamza1d(i) )+ &
        gupxx1d(i) * ( &
            Gamxxx1d(i) * gxxy1d(i) + Gamyxx1d(i) * gxyy1d(i) + Gamzxx1d(i) * gxzy1d(i) + &
            Gamxxy1d(i) * gxxx1d(i) + Gamyxy1d(i) * gxyx1d(i) + Gamzxy1d(i) * gxzx1d(i) + &
            Gamxxx1d(i) * gxyx1d(i) + Gamyxx1d(i) * gxyy1d(i) + Gamzxx1d(i) * gxyz1d(i) )+ &
        gupxy1d(i) * ( &
            Gamxxx1d(i) * gxyy1d(i) + Gamyxx1d(i) * gyyy1d(i) + Gamzxx1d(i) * gyzy1d(i) + &
            Gamxxy1d(i) * gxyx1d(i) + Gamyxy1d(i) * gyyx1d(i) + Gamzxy1d(i) * gyzx1d(i) + &
            Gamxxy1d(i) * gxyx1d(i) + Gamyxy1d(i) * gxyy1d(i) + Gamzxy1d(i) * gxyz1d(i) + &
            Gamxxy1d(i) * gxxy1d(i) + Gamyxy1d(i) * gxyy1d(i) + Gamzxy1d(i) * gxzy1d(i) + &
            Gamxyy1d(i) * gxxx1d(i) + Gamyyy1d(i) * gxyx1d(i) + Gamzyy1d(i) * gxzx1d(i) + &
            Gamxxx1d(i) * gyyx1d(i) + Gamyxx1d(i) * gyyy1d(i) + Gamzxx1d(i) * gyyz1d(i) )+ &
        gupxz1d(i) * ( &
            Gamxxx1d(i) * gxzy1d(i) + Gamyxx1d(i) * gyzy1d(i) + Gamzxx1d(i) * gzzy1d(i) + &
            Gamxxy1d(i) * gxzx1d(i) + Gamyxy1d(i) * gyzx1d(i) + Gamzxy1d(i) * gzzx1d(i) + &
            Gamxxz1d(i) * gxyx1d(i) + Gamyxz1d(i) * gxyy1d(i) + Gamzxz1d(i) * gxyz1d(i) + &
            Gamxxz1d(i) * gxxy1d(i) + Gamyxz1d(i) * gxyy1d(i) + Gamzxz1d(i) * gxzy1d(i) + &
            Gamxyz1d(i) * gxxx1d(i) + Gamyyz1d(i) * gxyx1d(i) + Gamzyz1d(i) * gxzx1d(i) + &
            Gamxxx1d(i) * gyzx1d(i) + Gamyxx1d(i) * gyzy1d(i) + Gamzxx1d(i) * gyzz1d(i) )+ &
        gupyy1d(i) * ( &
            Gamxxy1d(i) * gxyy1d(i) + Gamyxy1d(i) * gyyy1d(i) + Gamzxy1d(i) * gyzy1d(i) + &
            Gamxyy1d(i) * gxyx1d(i) + Gamyyy1d(i) * gyyx1d(i) + Gamzyy1d(i) * gyzx1d(i) + &
            Gamxxy1d(i) * gyyx1d(i) + Gamyxy1d(i) * gyyy1d(i) + Gamzxy1d(i) * gyyz1d(i) )+ &
        gupyz1d(i) * ( &
            Gamxxy1d(i) * gxzy1d(i) + Gamyxy1d(i) * gyzy1d(i) + Gamzxy1d(i) * gzzy1d(i) + &
            Gamxyy1d(i) * gxzx1d(i) + Gamyyy1d(i) * gyzx1d(i) + Gamzyy1d(i) * gzzx1d(i) + &
            Gamxxz1d(i) * gyyx1d(i) + Gamyxz1d(i) * gyyy1d(i) + Gamzxz1d(i) * gyyz1d(i) + &
            Gamxxz1d(i) * gxyy1d(i) + Gamyxz1d(i) * gyyy1d(i) + Gamzxz1d(i) * gyzy1d(i) + &
            Gamxyz1d(i) * gxyx1d(i) + Gamyyz1d(i) * gyyx1d(i) + Gamzyz1d(i) * gyzx1d(i) + &
            Gamxxy1d(i) * gyzx1d(i) + Gamyxy1d(i) * gyzy1d(i) + Gamzxy1d(i) * gyzz1d(i) )+ &
        gupzz1d(i) * ( &
            Gamxxz1d(i) * gxzy1d(i) + Gamyxz1d(i) * gyzy1d(i) + Gamzxz1d(i) * gzzy1d(i) + &
            Gamxyz1d(i) * gxzx1d(i) + Gamyyz1d(i) * gyzx1d(i) + Gamzyz1d(i) * gzzx1d(i) + &
            Gamxxz1d(i) * gyzx1d(i) + Gamyxz1d(i) * gyzy1d(i) + Gamzxz1d(i) * gyzz1d(i) )
      Rxz1d(i) = Rxz1d(i) + &
        HALF*( gxx1d(i) * Gamxz1d(i) + gxy1d(i) * Gamyz1d(i) + gxz1d(i) * Gamzz1d(i) + &
               gxz1d(i) * Gamxx1d(i) + gyz1d(i) * Gamyx1d(i) + gzz1d(i) * Gamzx1d(i) + &
               gxzx1d(i) * Gamxa1d(i) + gyzx1d(i) * Gamya1d(i) + gzzx1d(i) * Gamza1d(i) + &
               gxxz1d(i) * Gamxa1d(i) + gxyz1d(i) * Gamya1d(i) + gxzz1d(i) * Gamza1d(i) )+ &
        gupxx1d(i) * ( &
            Gamxxx1d(i) * gxxz1d(i) + Gamyxx1d(i) * gxyz1d(i) + Gamzxx1d(i) * gxzz1d(i) + &
            Gamxxz1d(i) * gxxx1d(i) + Gamyxz1d(i) * gxyx1d(i) + Gamzxz1d(i) * gxzx1d(i) + &
            Gamxxx1d(i) * gxzx1d(i) + Gamyxx1d(i) * gxzy1d(i) + Gamzxx1d(i) * gxzz1d(i) )+ &
        gupxy1d(i) * ( &
            Gamxxx1d(i) * gxyz1d(i) + Gamyxx1d(i) * gyyz1d(i) + Gamzxx1d(i) * gyzz1d(i) + &
            Gamxxz1d(i) * gxyx1d(i) + Gamyxz1d(i) * gyyx1d(i) + Gamzxz1d(i) * gyzx1d(i) + &
            Gamxxy1d(i) * gxzx1d(i) + Gamyxy1d(i) * gxzy1d(i) + Gamzxy1d(i) * gxzz1d(i) + &
            Gamxxy1d(i) * gxxz1d(i) + Gamyxy1d(i) * gxyz1d(i) + Gamzxy1d(i) * gxzz1d(i) + &
            Gamxyz1d(i) * gxxx1d(i) + Gamyyz1d(i) * gxyx1d(i) + Gamzyz1d(i) * gxzx1d(i) + &
            Gamxxx1d(i) * gyzx1d(i) + Gamyxx1d(i) * gyzy1d(i) + Gamzxx1d(i) * gyzz1d(i) )+ &
        gupxz1d(i) * ( &
            Gamxxx1d(i) * gxzz1d(i) + Gamyxx1d(i) * gyzz1d(i) + Gamzxx1d(i) * gzzz1d(i) + &
            Gamxxz1d(i) * gxzx1d(i) + Gamyxz1d(i) * gyzx1d(i) + Gamzxz1d(i) * gzzx1d(i) + &
            Gamxxz1d(i) * gxzx1d(i) + Gamyxz1d(i) * gxzy1d(i) + Gamzxz1d(i) * gxzz1d(i) + &
            Gamxxz1d(i) * gxxz1d(i) + Gamyxz1d(i) * gxyz1d(i) + Gamzxz1d(i) * gxzz1d(i) + &
            Gamxzz1d(i) * gxxx1d(i) + Gamyzz1d(i) * gxyx1d(i) + Gamzzz1d(i) * gxzx1d(i) + &
            Gamxxx1d(i) * gzzx1d(i) + Gamyxx1d(i) * gzzy1d(i) + Gamzxx1d(i) * gzzz1d(i) )+ &
        gupyy1d(i) * ( &
            Gamxxy1d(i) * gxyz1d(i) + Gamyxy1d(i) * gyyz1d(i) + Gamzxy1d(i) * gyzz1d(i) + &
            Gamxyz1d(i) * gxyx1d(i) + Gamyyz1d(i) * gyyx1d(i) + Gamzyz1d(i) * gyzx1d(i) + &
            Gamxxy1d(i) * gyzx1d(i) + Gamyxy1d(i) * gyzy1d(i) + Gamzxy1d(i) * gyzz1d(i) )+ &
        gupyz1d(i) * ( &
            Gamxxy1d(i) * gxzz1d(i) + Gamyxy1d(i) * gyzz1d(i) + Gamzxy1d(i) * gzzz1d(i) + &
            Gamxyz1d(i) * gxzx1d(i) + Gamyyz1d(i) * gyzx1d(i) + Gamzyz1d(i) * gzzx1d(i) + &
            Gamxxz1d(i) * gyzx1d(i) + Gamyxz1d(i) * gyzy1d(i) + Gamzxz1d(i) * gyzz1d(i) + &
            Gamxxz1d(i) * gxyz1d(i) + Gamyxz1d(i) * gyyz1d(i) + Gamzxz1d(i) * gyzz1d(i) + &
            Gamxzz1d(i) * gxyx1d(i) + Gamyzz1d(i) * gyyx1d(i) + Gamzzz1d(i) * gyzx1d(i) + &
            Gamxxy1d(i) * gzzx1d(i) + Gamyxy1d(i) * gzzy1d(i) + Gamzxy1d(i) * gzzz1d(i) )+ &
        gupzz1d(i) * ( &
            Gamxxz1d(i) * gxzz1d(i) + Gamyxz1d(i) * gyzz1d(i) + Gamzxz1d(i) * gzzz1d(i) + &
            Gamxzz1d(i) * gxzx1d(i) + Gamyzz1d(i) * gyzx1d(i) + Gamzzz1d(i) * gzzx1d(i) + &
            Gamxxz1d(i) * gzzx1d(i) + Gamyxz1d(i) * gzzy1d(i) + Gamzxz1d(i) * gzzz1d(i) )
      Ryz1d(i) = Ryz1d(i) + &
        HALF*( gxy1d(i) * Gamxz1d(i) + gyy1d(i) * Gamyz1d(i) + gyz1d(i) * Gamzz1d(i) + &
               gxz1d(i) * Gamxy1d(i) + gyz1d(i) * Gamyy1d(i) + gzz1d(i) * Gamzy1d(i) + &
               gxzy1d(i) * Gamxa1d(i) + gyzy1d(i) * Gamya1d(i) + gzzy1d(i) * Gamza1d(i) + &
               gxyz1d(i) * Gamxa1d(i) + gyyz1d(i) * Gamya1d(i) + gyzz1d(i) * Gamza1d(i) )+ &
        gupxx1d(i) * ( &
            Gamxxy1d(i) * gxxz1d(i) + Gamyxy1d(i) * gxyz1d(i) + Gamzxy1d(i) * gxzz1d(i) + &
            Gamxxz1d(i) * gxxy1d(i) + Gamyxz1d(i) * gxyy1d(i) + Gamzxz1d(i) * gxzy1d(i) + &
            Gamxxy1d(i) * gxzx1d(i) + Gamyxy1d(i) * gxzy1d(i) + Gamzxy1d(i) * gxzz1d(i) )+ &
        gupxy1d(i) * ( &
            Gamxxy1d(i) * gxyz1d(i) + Gamyxy1d(i) * gyyz1d(i) + Gamzxy1d(i) * gyzz1d(i) + &
            Gamxxz1d(i) * gxyy1d(i) + Gamyxz1d(i) * gyyy1d(i) + Gamzxz1d(i) * gyzy1d(i) + &
            Gamxyy1d(i) * gxzx1d(i) + Gamyyy1d(i) * gxzy1d(i) + Gamzyy1d(i) * gxzz1d(i) + &
            Gamxyy1d(i) * gxxz1d(i) + Gamyyy1d(i) * gxyz1d(i) + Gamzyy1d(i) * gxzz1d(i) + &
            Gamxyz1d(i) * gxxy1d(i) + Gamyyz1d(i) * gxyy1d(i) + Gamzyz1d(i) * gxzy1d(i) + &
            Gamxxy1d(i) * gyzx1d(i) + Gamyxy1d(i) * gyzy1d(i) + Gamzxy1d(i) * gyzz1d(i) )+ &
        gupxz1d(i) * ( &
            Gamxxy1d(i) * gxzz1d(i) + Gamyxy1d(i) * gyzz1d(i) + Gamzxy1d(i) * gzzz1d(i) + &
            Gamxxz1d(i) * gxzy1d(i) + Gamyxz1d(i) * gyzy1d(i) + Gamzxz1d(i) * gzzy1d(i) + &
            Gamxyz1d(i) * gxzx1d(i) + Gamyyz1d(i) * gxzy1d(i) + Gamzyz1d(i) * gxzz1d(i) + &
            Gamxyz1d(i) * gxxz1d(i) + Gamyyz1d(i) * gxyz1d(i) + Gamzyz1d(i) * gxzz1d(i) + &
            Gamxzz1d(i) * gxxy1d(i) + Gamyzz1d(i) * gxyy1d(i) + Gamzzz1d(i) * gxzy1d(i) + &
            Gamxxy1d(i) * gzzx1d(i) + Gamyxy1d(i) * gzzy1d(i) + Gamzxy1d(i) * gzzz1d(i) )+ &
        gupyy1d(i) * ( &
            Gamxyy1d(i) * gxyz1d(i) + Gamyyy1d(i) * gyyz1d(i) + Gamzyy1d(i) * gyzz1d(i) + &
            Gamxyz1d(i) * gxyy1d(i) + Gamyyz1d(i) * gyyy1d(i) + Gamzyz1d(i) * gyzy1d(i) + &
            Gamxyy1d(i) * gyzx1d(i) + Gamyyy1d(i) * gyzy1d(i) + Gamzyy1d(i) * gyzz1d(i) )+ &
        gupyz1d(i) * ( &
            Gamxyy1d(i) * gxzz1d(i) + Gamyyy1d(i) * gyzz1d(i) + Gamzyy1d(i) * gzzz1d(i) + &
            Gamxyz1d(i) * gxzy1d(i) + Gamyyz1d(i) * gyzy1d(i) + Gamzyz1d(i) * gzzy1d(i) + &
            Gamxyz1d(i) * gyzx1d(i) + Gamyyz1d(i) * gyzy1d(i) + Gamzyz1d(i) * gyzz1d(i) + &
            Gamxyz1d(i) * gxyz1d(i) + Gamyyz1d(i) * gyyz1d(i) + Gamzyz1d(i) * gyzz1d(i) + &
            Gamxzz1d(i) * gxyy1d(i) + Gamyzz1d(i) * gyyy1d(i) + Gamzzz1d(i) * gyzy1d(i) + &
            Gamxyy1d(i) * gzzx1d(i) + Gamyyy1d(i) * gzzy1d(i) + Gamzyy1d(i) * gzzz1d(i) )+ &
        gupzz1d(i) * ( &
            Gamxyz1d(i) * gxzz1d(i) + Gamyyz1d(i) * gyzz1d(i) + Gamzyz1d(i) * gzzz1d(i) + &
            Gamxzz1d(i) * gxzy1d(i) + Gamyzz1d(i) * gyzy1d(i) + Gamzzz1d(i) * gzzy1d(i) + &
            Gamxyz1d(i) * gzzx1d(i) + Gamyyz1d(i) * gzzy1d(i) + Gamzyz1d(i) * gzzz1d(i) )
    end do
  end block
!covariant second derivative of chi respect to tilted metric
  call fdderivs(ex,chi,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)

block
  integer :: i,j,k
  real(8) :: cx, cy, cz

  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)

        cx = chix(i,j,k)
        cy = chiy(i,j,k)
        cz = chiz(i,j,k)

        fxx(i,j,k) = fxx(i,j,k) - Gamxxx(i,j,k)*cx - Gamyxx(i,j,k)*cy - Gamzxx(i,j,k)*cz
        fxy(i,j,k) = fxy(i,j,k) - Gamxxy(i,j,k)*cx - Gamyxy(i,j,k)*cy - Gamzxy(i,j,k)*cz
        fxz(i,j,k) = fxz(i,j,k) - Gamxxz(i,j,k)*cx - Gamyxz(i,j,k)*cy - Gamzxz(i,j,k)*cz
        fyy(i,j,k) = fyy(i,j,k) - Gamxyy(i,j,k)*cx - Gamyyy(i,j,k)*cy - Gamzyy(i,j,k)*cz
        fyz(i,j,k) = fyz(i,j,k) - Gamxyz(i,j,k)*cx - Gamyyz(i,j,k)*cy - Gamzyz(i,j,k)*cz
        fzz(i,j,k) = fzz(i,j,k) - Gamxzz(i,j,k)*cx - Gamyzz(i,j,k)*cy - Gamzzz(i,j,k)*cz

      end do
    end do
  end do

end block
! Store D^l D_l chi - 3/(2*chi) D^l chi D_l chi in f

  f =        gupxx * ( fxx - F3o2/chin1 * chix * chix ) + &
             gupyy * ( fyy - F3o2/chin1 * chiy * chiy ) + &
             gupzz * ( fzz - F3o2/chin1 * chiz * chiz ) + &
       TWO * gupxy * ( fxy - F3o2/chin1 * chix * chiy ) + &
       TWO * gupxz * ( fxz - F3o2/chin1 * chix * chiz ) + &
       TWO * gupyz * ( fyz - F3o2/chin1 * chiy * chiz ) 
! Add chi part to Ricci tensor:

block
  integer :: i,j,k
  real(8) :: cx, cy, cz
  real(8) :: inv_chi, inv2, tmp

  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)

        cx = chix(i,j,k)
        cy = chiy(i,j,k)
        cz = chiz(i,j,k)

        inv_chi = 1.d0 / chin1(i,j,k)
        inv2     = 0.5d0 * inv_chi

        ! common scalar part: f(i,j,k) already computed
        tmp = f(i,j,k)

        Rxx(i,j,k) = Rxx(i,j,k) + &
          ( fxx(i,j,k) - cx*cx*inv2 + gxx(i,j,k)*tmp ) * inv2

        Ryy(i,j,k) = Ryy(i,j,k) + &
          ( fyy(i,j,k) - cy*cy*inv2 + gyy(i,j,k)*tmp ) * inv2

        Rzz(i,j,k) = Rzz(i,j,k) + &
          ( fzz(i,j,k) - cz*cz*inv2 + gzz(i,j,k)*tmp ) * inv2

        Rxy(i,j,k) = Rxy(i,j,k) + &
          ( fxy(i,j,k) - cx*cy*inv2 + gxy(i,j,k)*tmp ) * inv2

        Rxz(i,j,k) = Rxz(i,j,k) + &
          ( fxz(i,j,k) - cx*cz*inv2 + gxz(i,j,k)*tmp ) * inv2

        Ryz(i,j,k) = Ryz(i,j,k) + &
          ( fyz(i,j,k) - cy*cz*inv2 + gyz(i,j,k)*tmp ) * inv2

      end do
    end do
  end do

end block

! covariant second derivatives of the lapse respect to physical metric
  call fdderivs(ex,Lap,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z, &
                SYM,SYM,SYM,symmetry,Lev)

  gxxx = (gupxx * chix + gupxy * chiy + gupxz * chiz)/chin1
  gxxy = (gupxy * chix + gupyy * chiy + gupyz * chiz)/chin1
  gxxz = (gupxz * chix + gupyz * chiy + gupzz * chiz)/chin1
! now get physical second kind of connection
  Gamxxx = Gamxxx - ( (chix + chix)/chin1 - gxx * gxxx )*HALF
  Gamyxx = Gamyxx - (                     - gxx * gxxy )*HALF
  Gamzxx = Gamzxx - (                     - gxx * gxxz )*HALF
  Gamxyy = Gamxyy - (                     - gyy * gxxx )*HALF
  Gamyyy = Gamyyy - ( (chiy + chiy)/chin1 - gyy * gxxy )*HALF
  Gamzyy = Gamzyy - (                     - gyy * gxxz )*HALF
  Gamxzz = Gamxzz - (                     - gzz * gxxx )*HALF
  Gamyzz = Gamyzz - (                     - gzz * gxxy )*HALF
  Gamzzz = Gamzzz - ( (chiz + chiz)/chin1 - gzz * gxxz )*HALF
  Gamxxy = Gamxxy - (  chiy        /chin1 - gxy * gxxx )*HALF
  Gamyxy = Gamyxy - (         chix /chin1 - gxy * gxxy )*HALF
  Gamzxy = Gamzxy - (                     - gxy * gxxz )*HALF
  Gamxxz = Gamxxz - (  chiz        /chin1 - gxz * gxxx )*HALF
  Gamyxz = Gamyxz - (                     - gxz * gxxy )*HALF
  Gamzxz = Gamzxz - (         chix /chin1 - gxz * gxxz )*HALF
  Gamxyz = Gamxyz - (                     - gyz * gxxx )*HALF
  Gamyyz = Gamyyz - (  chiz        /chin1 - gyz * gxxy )*HALF
  Gamzyz = Gamzyz - (         chiy /chin1 - gyz * gxxz )*HALF

  fxx = fxx - Gamxxx*Lapx - Gamyxx*Lapy - Gamzxx*Lapz
  fyy = fyy - Gamxyy*Lapx - Gamyyy*Lapy - Gamzyy*Lapz
  fzz = fzz - Gamxzz*Lapx - Gamyzz*Lapy - Gamzzz*Lapz
  fxy = fxy - Gamxxy*Lapx - Gamyxy*Lapy - Gamzxy*Lapz
  fxz = fxz - Gamxxz*Lapx - Gamyxz*Lapy - Gamzxz*Lapz
  fyz = fyz - Gamxyz*Lapx - Gamyyz*Lapy - Gamzyz*Lapz

! store D^i D_i Lap in trK_rhs upto chi
  trK_rhs =    gupxx * fxx + gupyy * fyy + gupzz * fzz + &
        TWO* ( gupxy * fxy + gupxz * fxz + gupyz * fyz )
#if 1        
!! follow bam code
  S =  chin1 * ( gupxx * Sxx + gupyy * Syy + gupzz * Szz + &
     TWO * ( gupxy * Sxy + gupxz * Sxz + gupyz * Syz ) )
  f = F2o3 * trK * trK -(&
       gupxx * ( &
       gupxx * Axx * Axx + gupyy * Axy * Axy + gupzz * Axz * Axz + &
       TWO * (gupxy * Axx * Axy + gupxz * Axx * Axz + gupyz * Axy * Axz) ) + &
       gupyy * ( &
       gupxx * Axy * Axy + gupyy * Ayy * Ayy + gupzz * Ayz * Ayz + &
       TWO * (gupxy * Axy * Ayy + gupxz * Axy * Ayz + gupyz * Ayy * Ayz) ) + &
       gupzz * ( &
       gupxx * Axz * Axz + gupyy * Ayz * Ayz + gupzz * Azz * Azz + &
       TWO * (gupxy * Axz * Ayz + gupxz * Axz * Azz + gupyz * Ayz * Azz) ) + &
       TWO * ( &
       gupxy * ( &
       gupxx * Axx * Axy + gupyy * Axy * Ayy + gupzz * Axz * Ayz + &
       gupxy * (Axx * Ayy + Axy * Axy) + &
       gupxz * (Axx * Ayz + Axz * Axy) + &
       gupyz * (Axy * Ayz + Axz * Ayy) ) + &
       gupxz * ( &
       gupxx * Axx * Axz + gupyy * Axy * Ayz + gupzz * Axz * Azz + &
       gupxy * (Axx * Ayz + Axy * Axz) + &
       gupxz * (Axx * Azz + Axz * Axz) + &
       gupyz * (Axy * Azz + Axz * Ayz) ) + &
       gupyz * ( &
       gupxx * Axy * Axz + gupyy * Ayy * Ayz + gupzz * Ayz * Azz + &
       gupxy * (Axy * Ayz + Ayy * Axz) + &
       gupxz * (Axy * Azz + Ayz * Axz) + &
       gupyz * (Ayy * Azz + Ayz * Ayz) ) )) -1.6d1*PI*rho + EIGHT * PI * S
  f = - F1o3 *(  gupxx * fxx + gupyy * fyy + gupzz * fzz + &
        TWO* ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) + alpn1/chin1*f)
  
  fxx = alpn1 * (Rxx - EIGHT * PI * Sxx) - fxx
  fxy = alpn1 * (Rxy - EIGHT * PI * Sxy) - fxy
  fxz = alpn1 * (Rxz - EIGHT * PI * Sxz) - fxz
  fyy = alpn1 * (Ryy - EIGHT * PI * Syy) - fyy
  fyz = alpn1 * (Ryz - EIGHT * PI * Syz) - fyz
  fzz = alpn1 * (Rzz - EIGHT * PI * Szz) - fzz
#else        
! Add lapse and S_ij parts to Ricci tensor:

  fxx = alpn1 * (Rxx - EIGHT * PI * Sxx) - fxx
  fxy = alpn1 * (Rxy - EIGHT * PI * Sxy) - fxy
  fxz = alpn1 * (Rxz - EIGHT * PI * Sxz) - fxz
  fyy = alpn1 * (Ryy - EIGHT * PI * Syy) - fyy
  fyz = alpn1 * (Ryz - EIGHT * PI * Syz) - fyz
  fzz = alpn1 * (Rzz - EIGHT * PI * Szz) - fzz

! Compute trace-free part (note: chi^-1 and chi cancel!):

  f = F1o3 *(  gupxx * fxx + gupyy * fyy + gupzz * fzz + &
        TWO* ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) )
#endif

  Axx_rhs = fxx - gxx * f
  Ayy_rhs = fyy - gyy * f
  Azz_rhs = fzz - gzz * f
  Axy_rhs = fxy - gxy * f
  Axz_rhs = fxz - gxz * f
  Ayz_rhs = fyz - gyz * f

! Now: store A_il A^l_j into fij:

  fxx =       gupxx * Axx * Axx + gupyy * Axy * Axy + gupzz * Axz * Axz + &
       TWO * (gupxy * Axx * Axy + gupxz * Axx * Axz + gupyz * Axy * Axz)
  fyy =       gupxx * Axy * Axy + gupyy * Ayy * Ayy + gupzz * Ayz * Ayz + &
       TWO * (gupxy * Axy * Ayy + gupxz * Axy * Ayz + gupyz * Ayy * Ayz)
  fzz =       gupxx * Axz * Axz + gupyy * Ayz * Ayz + gupzz * Azz * Azz + &
       TWO * (gupxy * Axz * Ayz + gupxz * Axz * Azz + gupyz * Ayz * Azz)
  fxy =       gupxx * Axx * Axy + gupyy * Axy * Ayy + gupzz * Axz * Ayz + &
              gupxy *(Axx * Ayy + Axy * Axy)                            + &
              gupxz *(Axx * Ayz + Axz * Axy)                            + &
              gupyz *(Axy * Ayz + Axz * Ayy)
  fxz =       gupxx * Axx * Axz + gupyy * Axy * Ayz + gupzz * Axz * Azz + &
              gupxy *(Axx * Ayz + Axy * Axz)                            + &
              gupxz *(Axx * Azz + Axz * Axz)                            + &
              gupyz *(Axy * Azz + Axz * Ayz)
  fyz =       gupxx * Axy * Axz + gupyy * Ayy * Ayz + gupzz * Ayz * Azz + &
              gupxy *(Axy * Ayz + Ayy * Axz)                            + &
              gupxz *(Axy * Azz + Ayz * Axz)                            + &
              gupyz *(Ayy * Azz + Ayz * Ayz)

  f = chin1
! store D^i D_i Lap in trK_rhs
  trK_rhs = f*trK_rhs
          
  Axx_rhs =           f * Axx_rhs+ alpn1 * (trK * Axx - TWO * fxx)  + &
           TWO * (  Axx * betaxx +   Axy * betayx +   Axz * betazx )- &
             F2o3 * Axx * div_beta

  Ayy_rhs =           f * Ayy_rhs+ alpn1 * (trK * Ayy - TWO * fyy)  + &
           TWO * (  Axy * betaxy +   Ayy * betayy +   Ayz * betazy )- &
             F2o3 * Ayy * div_beta

  Azz_rhs =           f * Azz_rhs+ alpn1 * (trK * Azz - TWO * fzz)  + &
           TWO * (  Axz * betaxz +   Ayz * betayz +   Azz * betazz )- &
             F2o3 * Azz * div_beta

  Axy_rhs =           f * Axy_rhs+ alpn1 *( trK * Axy  - TWO * fxy )+ &
                    Axx * betaxy                  +   Axz * betazy  + &
                                     Ayy * betayx +   Ayz * betazx  + &
             F1o3 * Axy * div_beta                -   Axy * betazz

  Ayz_rhs =           f * Ayz_rhs+ alpn1 *( trK * Ayz  - TWO * fyz )+ &
                    Axy * betaxz +   Ayy * betayz                   + &
                    Axz * betaxy                  +   Azz * betazy  + &
             F1o3 * Ayz * div_beta                -   Ayz * betaxx
 
  Axz_rhs =           f * Axz_rhs+ alpn1 *( trK * Axz  - TWO * fxz )+ &
                    Axx * betaxz +   Axy * betayz                   + &
                                     Ayz * betayx +   Azz * betazx  + &
             F1o3 * Axz * div_beta                -   Axz * betayy      !rhs for Aij

! Compute trace of S_ij

  S =  f * ( gupxx * Sxx + gupyy * Syy + gupzz * Szz + &
     TWO * ( gupxy * Sxy + gupxz * Sxz + gupyz * Syz ) )

  trK_rhs = - trK_rhs + alpn1 *( F1o3 * trK * trK         + &
                gupxx * fxx + gupyy * fyy + gupzz * fzz   + &
        TWO * ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) + &
       FOUR * PI * ( rho + S ))                                !rhs for trK
  
!!!! gauge variable part

  Lap_rhs = -TWO*alpn1*trK
#if (GAUGE == 0)
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  dtSfx_rhs = Gamx_rhs - eta*dtSfx
  dtSfy_rhs = Gamy_rhs - eta*dtSfy
  dtSfz_rhs = Gamz_rhs - eta*dtSfz
#elif (GAUGE == 1)
  betax_rhs = Gamx - eta*betax
  betay_rhs = Gamy - eta*betay
  betaz_rhs = Gamz - eta*betaz

  dtSfx_rhs = ZEO
  dtSfy_rhs = ZEO
  dtSfz_rhs = ZEO
#elif (GAUGE == 2)
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  call fderivs(ex,chi,dtSfx_rhs,dtSfy_rhs,dtSfz_rhs,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  reta = gupxx * dtSfx_rhs * dtSfx_rhs + gupyy * dtSfy_rhs * dtSfy_rhs + gupzz * dtSfz_rhs * dtSfz_rhs + &
       TWO * (gupxy * dtSfx_rhs * dtSfy_rhs + gupxz * dtSfx_rhs * dtSfz_rhs + gupyz * dtSfy_rhs * dtSfz_rhs)
  reta = 1.31d0/2*dsqrt(reta/chin1)/(1-dsqrt(chin1))**2
  dtSfx_rhs = Gamx_rhs - reta*dtSfx
  dtSfy_rhs = Gamy_rhs - reta*dtSfy
  dtSfz_rhs = Gamz_rhs - reta*dtSfz
#elif (GAUGE == 3)
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  call fderivs(ex,chi,dtSfx_rhs,dtSfy_rhs,dtSfz_rhs,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  reta = gupxx * dtSfx_rhs * dtSfx_rhs + gupyy * dtSfy_rhs * dtSfy_rhs + gupzz * dtSfz_rhs * dtSfz_rhs + &
       TWO * (gupxy * dtSfx_rhs * dtSfy_rhs + gupxz * dtSfx_rhs * dtSfz_rhs + gupyz * dtSfy_rhs * dtSfz_rhs)
  reta = 1.31d0/2*dsqrt(reta/chin1)/(1-chin1)**2
  dtSfx_rhs = Gamx_rhs - reta*dtSfx
  dtSfy_rhs = Gamy_rhs - reta*dtSfy
  dtSfz_rhs = Gamz_rhs - reta*dtSfz
#elif (GAUGE == 4)
  call fderivs(ex,chi,dtSfx_rhs,dtSfy_rhs,dtSfz_rhs,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  reta = gupxx * dtSfx_rhs * dtSfx_rhs + gupyy * dtSfy_rhs * dtSfy_rhs + gupzz * dtSfz_rhs * dtSfz_rhs + &
       TWO * (gupxy * dtSfx_rhs * dtSfy_rhs + gupxz * dtSfx_rhs * dtSfz_rhs + gupyz * dtSfy_rhs * dtSfz_rhs)
  reta = 1.31d0/2*dsqrt(reta/chin1)/(1-dsqrt(chin1))**2
  betax_rhs = FF*Gamx - reta*betax
  betay_rhs = FF*Gamy - reta*betay
  betaz_rhs = FF*Gamz - reta*betaz

  dtSfx_rhs = ZEO
  dtSfy_rhs = ZEO
  dtSfz_rhs = ZEO
#elif (GAUGE == 5)
  call fderivs(ex,chi,dtSfx_rhs,dtSfy_rhs,dtSfz_rhs,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  reta = gupxx * dtSfx_rhs * dtSfx_rhs + gupyy * dtSfy_rhs * dtSfy_rhs + gupzz * dtSfz_rhs * dtSfz_rhs + &
       TWO * (gupxy * dtSfx_rhs * dtSfy_rhs + gupxz * dtSfx_rhs * dtSfz_rhs + gupyz * dtSfy_rhs * dtSfz_rhs)
  reta = 1.31d0/2*dsqrt(reta/chin1)/(1-chin1)**2
  betax_rhs = FF*Gamx - reta*betax
  betay_rhs = FF*Gamy - reta*betay
  betaz_rhs = FF*Gamz - reta*betaz

  dtSfx_rhs = ZEO
  dtSfy_rhs = ZEO
  dtSfz_rhs = ZEO
#elif (GAUGE == 6)
  if(BHN==2)then
   M = Mass(1)+Mass(2)
   A = 2.d0/M
   w1 = 1.2d1
   w2 = w1
   C1 = 1.d0/Mass(1) - A
   C2 = 1.d0/Mass(2) - A

   do k=1,ex(3)
   do j=1,ex(2)
   do i=1,ex(1)
     r1 = ((Porg(1)-X(i))**2+(Porg(2)-Y(j))**2+(Porg(3)-Z(k))**2)/ &
          ((Porg(1)-Porg(4))**2+(Porg(2)-Porg(5))**2+(Porg(3)-Porg(6))**2)
     r2 = ((Porg(4)-X(i))**2+(Porg(5)-Y(j))**2+(Porg(6)-Z(k))**2)/ &
          ((Porg(1)-Porg(4))**2+(Porg(2)-Porg(5))**2+(Porg(3)-Porg(6))**2)
     reta(i,j,k) = A + C1/(ONE+w1*r1) + C2/(ONE+w2*r2)
    enddo
    enddo
    enddo
  else
    write(*,*) "not support BH_num in Jason's form 1",BHN
  endif
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  dtSfx_rhs = Gamx_rhs - reta*dtSfx
  dtSfy_rhs = Gamy_rhs - reta*dtSfy
  dtSfz_rhs = Gamz_rhs - reta*dtSfz
#elif (GAUGE == 7)
  if(BHN==2)then
   M = Mass(1)+Mass(2)
   A = 2.d0/M
   w1 = 1.2d1
   w2 = w1
   C1 = 1.d0/Mass(1) - A
   C2 = 1.d0/Mass(2) - A

   do k=1,ex(3)
   do j=1,ex(2)
   do i=1,ex(1)
     r1 = ((Porg(1)-X(i))**2+(Porg(2)-Y(j))**2+(Porg(3)-Z(k))**2)/ &
          ((Porg(1)-Porg(4))**2+(Porg(2)-Porg(5))**2+(Porg(3)-Porg(6))**2)
     r2 = ((Porg(4)-X(i))**2+(Porg(5)-Y(j))**2+(Porg(6)-Z(k))**2)/ &
          ((Porg(1)-Porg(4))**2+(Porg(2)-Porg(5))**2+(Porg(3)-Porg(6))**2)
     reta(i,j,k) = A + C1*dexp(-w1*r1) + C2*dexp(-w2*r2)
    enddo
    enddo
    enddo
  else
    write(*,*) "not support BH_num in Jason's form 2",BHN
  endif
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  dtSfx_rhs = Gamx_rhs - reta*dtSfx
  dtSfy_rhs = Gamy_rhs - reta*dtSfy
  dtSfz_rhs = Gamz_rhs - reta*dtSfz
#endif  

  SSS(1)=SYM
  SSS(2)=SYM
  SSS(3)=SYM

  AAS(1)=ANTI
  AAS(2)=ANTI
  AAS(3)=SYM

  ASA(1)=ANTI
  ASA(2)=SYM
  ASA(3)=ANTI

  SAA(1)=SYM
  SAA(2)=ANTI
  SAA(3)=ANTI

  ASS(1)=ANTI
  ASS(2)=SYM
  ASS(3)=SYM

  SAS(1)=SYM
  SAS(2)=ANTI
  SAS(3)=SYM

  SSA(1)=SYM
  SSA(2)=SYM
  SSA(3)=ANTI

!!!!!!!!!advection term part

  call lopsided(ex,X,Y,Z,gxx,gxx_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,gxy,gxy_rhs,betax,betay,betaz,Symmetry,AAS)
  call lopsided(ex,X,Y,Z,gxz,gxz_rhs,betax,betay,betaz,Symmetry,ASA)
  call lopsided(ex,X,Y,Z,gyy,gyy_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,gyz,gyz_rhs,betax,betay,betaz,Symmetry,SAA)
  call lopsided(ex,X,Y,Z,gzz,gzz_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,Axx,Axx_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,Axy,Axy_rhs,betax,betay,betaz,Symmetry,AAS)
  call lopsided(ex,X,Y,Z,Axz,Axz_rhs,betax,betay,betaz,Symmetry,ASA)
  call lopsided(ex,X,Y,Z,Ayy,Ayy_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,Ayz,Ayz_rhs,betax,betay,betaz,Symmetry,SAA)
  call lopsided(ex,X,Y,Z,Azz,Azz_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,chi,chi_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,trK,trK_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,Gamx,Gamx_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,Gamy,Gamy_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,Gamz,Gamz_rhs,betax,betay,betaz,Symmetry,SSA)
!!
  call lopsided(ex,X,Y,Z,Lap,Lap_rhs,betax,betay,betaz,Symmetry,SSS)

#if (GAUGE == 0 || GAUGE == 1 || GAUGE == 2 || GAUGE == 3 || GAUGE == 4 || GAUGE == 5 || GAUGE == 6 || GAUGE == 7)
  call lopsided(ex,X,Y,Z,betax,betax_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,betay,betay_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,betaz,betaz_rhs,betax,betay,betaz,Symmetry,SSA)
#endif

#if (GAUGE == 0 || GAUGE == 2 || GAUGE == 3 || GAUGE == 6 || GAUGE == 7)
  call lopsided(ex,X,Y,Z,dtSfx,dtSfx_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,dtSfy,dtSfy_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,dtSfz,dtSfz_rhs,betax,betay,betaz,Symmetry,SSA)
#endif

  if(eps>0)then 
! usual Kreiss-Oliger dissipation      
  call kodis(ex,X,Y,Z,chi,chi_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,trK,trK_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dxx,gxx_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gxy,gxy_rhs,AAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gxz,gxz_rhs,ASA,Symmetry,eps)
  call kodis(ex,X,Y,Z,dyy,gyy_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gyz,gyz_rhs,SAA,Symmetry,eps)
  call kodis(ex,X,Y,Z,dzz,gzz_rhs,SSS,Symmetry,eps)
#if 0
#define i 42
#define j 40
#define k 40
if(Lev == 1)then
write(*,*) X(i),Y(j),Z(k)
write(*,*) "before",Axx_rhs(i,j,k)
endif
#undef i
#undef j
#undef k
!!stop
#endif
  call kodis(ex,X,Y,Z,Axx,Axx_rhs,SSS,Symmetry,eps)
#if 0
#define i 42
#define j 40
#define k 40
if(Lev == 1)then
write(*,*) X(i),Y(j),Z(k)
write(*,*) "after",Axx_rhs(i,j,k)
endif
#undef i
#undef j
#undef k
!!stop
#endif
  call kodis(ex,X,Y,Z,Axy,Axy_rhs,AAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Axz,Axz_rhs,ASA,Symmetry,eps)
  call kodis(ex,X,Y,Z,Ayy,Ayy_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Ayz,Ayz_rhs,SAA,Symmetry,eps)
  call kodis(ex,X,Y,Z,Azz,Azz_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamx,Gamx_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamy,Gamy_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamz,Gamz_rhs,SSA,Symmetry,eps)

#if 1 
!! bam does not apply dissipation on gauge variables
  call kodis(ex,X,Y,Z,Lap,Lap_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betax,betax_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betay,betay_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betaz,betaz_rhs,SSA,Symmetry,eps)
#if (GAUGE == 0 || GAUGE == 2 || GAUGE == 3 || GAUGE == 6 || GAUGE == 7)
  call kodis(ex,X,Y,Z,dtSfx,dtSfx_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfy,dtSfy_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfz,dtSfz_rhs,SSA,Symmetry,eps)
#endif
#endif

  endif

  if(co == 0)then
! ham_Res = trR + 2/3 * K^2 - A_ij * A^ij - 16 * PI * rho
! here trR is respect to physical metric
  ham_Res =   gupxx * Rxx + gupyy * Ryy + gupzz * Rzz + &
        TWO* ( gupxy * Rxy + gupxz * Rxz + gupyz * Ryz )

  ham_Res = chin1*ham_Res + F2o3 * trK * trK -(&
       gupxx * ( &
       gupxx * Axx * Axx + gupyy * Axy * Axy + gupzz * Axz * Axz + &
       TWO * (gupxy * Axx * Axy + gupxz * Axx * Axz + gupyz * Axy * Axz) ) + &
       gupyy * ( &
       gupxx * Axy * Axy + gupyy * Ayy * Ayy + gupzz * Ayz * Ayz + &
       TWO * (gupxy * Axy * Ayy + gupxz * Axy * Ayz + gupyz * Ayy * Ayz) ) + &
       gupzz * ( &
       gupxx * Axz * Axz + gupyy * Ayz * Ayz + gupzz * Azz * Azz + &
       TWO * (gupxy * Axz * Ayz + gupxz * Axz * Azz + gupyz * Ayz * Azz) ) + &
       TWO * ( &
       gupxy * ( &
       gupxx * Axx * Axy + gupyy * Axy * Ayy + gupzz * Axz * Ayz + &
       gupxy * (Axx * Ayy + Axy * Axy) + &
       gupxz * (Axx * Ayz + Axz * Axy) + &
       gupyz * (Axy * Ayz + Axz * Ayy) ) + &
       gupxz * ( &
       gupxx * Axx * Axz + gupyy * Axy * Ayz + gupzz * Axz * Azz + &
       gupxy * (Axx * Ayz + Axy * Axz) + &
       gupxz * (Axx * Azz + Axz * Axz) + &
       gupyz * (Axy * Azz + Axz * Ayz) ) + &
       gupyz * ( &
       gupxx * Axy * Axz + gupyy * Ayy * Ayz + gupzz * Ayz * Azz + &
       gupxy * (Axy * Ayz + Ayy * Axz) + &
       gupxz * (Axy * Azz + Ayz * Axz) + &
       gupyz * (Ayy * Azz + Ayz * Ayz) ) ))- F16 * PI * rho

! mov_Res_j = gupkj*(-1/chi d_k chi*A_ij + D_k A_ij) - 2/3 d_j trK - 8 PI s_j where D respect to physical metric
! store D_i A_jk - 1/chi d_i chi*A_jk in gjk_i
  call fderivs(ex,Axx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Axy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,0)
  call fderivs(ex,Axz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,0)
  call fderivs(ex,Ayy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Ayz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,0)
  call fderivs(ex,Azz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)

  gxxx = gxxx - (  Gamxxx * Axx + Gamyxx * Axy + Gamzxx * Axz &
                 + Gamxxx * Axx + Gamyxx * Axy + Gamzxx * Axz) - chix*Axx/chin1
  gxyx = gxyx - (  Gamxxy * Axx + Gamyxy * Axy + Gamzxy * Axz &
                 + Gamxxx * Axy + Gamyxx * Ayy + Gamzxx * Ayz) - chix*Axy/chin1
  gxzx = gxzx - (  Gamxxz * Axx + Gamyxz * Axy + Gamzxz * Axz &
                 + Gamxxx * Axz + Gamyxx * Ayz + Gamzxx * Azz) - chix*Axz/chin1
  gyyx = gyyx - (  Gamxxy * Axy + Gamyxy * Ayy + Gamzxy * Ayz &
                 + Gamxxy * Axy + Gamyxy * Ayy + Gamzxy * Ayz) - chix*Ayy/chin1
  gyzx = gyzx - (  Gamxxz * Axy + Gamyxz * Ayy + Gamzxz * Ayz &
                 + Gamxxy * Axz + Gamyxy * Ayz + Gamzxy * Azz) - chix*Ayz/chin1
  gzzx = gzzx - (  Gamxxz * Axz + Gamyxz * Ayz + Gamzxz * Azz &
                 + Gamxxz * Axz + Gamyxz * Ayz + Gamzxz * Azz) - chix*Azz/chin1
  gxxy = gxxy - (  Gamxxy * Axx + Gamyxy * Axy + Gamzxy * Axz &
                 + Gamxxy * Axx + Gamyxy * Axy + Gamzxy * Axz) - chiy*Axx/chin1
  gxyy = gxyy - (  Gamxyy * Axx + Gamyyy * Axy + Gamzyy * Axz &
                 + Gamxxy * Axy + Gamyxy * Ayy + Gamzxy * Ayz) - chiy*Axy/chin1
  gxzy = gxzy - (  Gamxyz * Axx + Gamyyz * Axy + Gamzyz * Axz &
                 + Gamxxy * Axz + Gamyxy * Ayz + Gamzxy * Azz) - chiy*Axz/chin1
  gyyy = gyyy - (  Gamxyy * Axy + Gamyyy * Ayy + Gamzyy * Ayz &
                 + Gamxyy * Axy + Gamyyy * Ayy + Gamzyy * Ayz) - chiy*Ayy/chin1
  gyzy = gyzy - (  Gamxyz * Axy + Gamyyz * Ayy + Gamzyz * Ayz &
                 + Gamxyy * Axz + Gamyyy * Ayz + Gamzyy * Azz) - chiy*Ayz/chin1
  gzzy = gzzy - (  Gamxyz * Axz + Gamyyz * Ayz + Gamzyz * Azz &
                 + Gamxyz * Axz + Gamyyz * Ayz + Gamzyz * Azz) - chiy*Azz/chin1
  gxxz = gxxz - (  Gamxxz * Axx + Gamyxz * Axy + Gamzxz * Axz &
                 + Gamxxz * Axx + Gamyxz * Axy + Gamzxz * Axz) - chiz*Axx/chin1
  gxyz = gxyz - (  Gamxyz * Axx + Gamyyz * Axy + Gamzyz * Axz &
                 + Gamxxz * Axy + Gamyxz * Ayy + Gamzxz * Ayz) - chiz*Axy/chin1
  gxzz = gxzz - (  Gamxzz * Axx + Gamyzz * Axy + Gamzzz * Axz &
                 + Gamxxz * Axz + Gamyxz * Ayz + Gamzxz * Azz) - chiz*Axz/chin1
  gyyz = gyyz - (  Gamxyz * Axy + Gamyyz * Ayy + Gamzyz * Ayz &
                 + Gamxyz * Axy + Gamyyz * Ayy + Gamzyz * Ayz) - chiz*Ayy/chin1
  gyzz = gyzz - (  Gamxzz * Axy + Gamyzz * Ayy + Gamzzz * Ayz &
                 + Gamxyz * Axz + Gamyyz * Ayz + Gamzyz * Azz) - chiz*Ayz/chin1
  gzzz = gzzz - (  Gamxzz * Axz + Gamyzz * Ayz + Gamzzz * Azz &
                 + Gamxzz * Axz + Gamyzz * Ayz + Gamzzz * Azz) - chiz*Azz/chin1
movx_Res = gupxx*gxxx + gupyy*gxyy + gupzz*gxzz &
          +gupxy*gxyx + gupxz*gxzx + gupyz*gxzy &
          +gupxy*gxxy + gupxz*gxxz + gupyz*gxyz
movy_Res = gupxx*gxyx + gupyy*gyyy + gupzz*gyzz &
          +gupxy*gyyx + gupxz*gyzx + gupyz*gyzy &
          +gupxy*gxyy + gupxz*gxyz + gupyz*gyyz
movz_Res = gupxx*gxzx + gupyy*gyzy + gupzz*gzzz &
          +gupxy*gyzx + gupxz*gzzx + gupyz*gzzy &
          +gupxy*gxzy + gupxz*gxzz + gupyz*gyzz

movx_Res = movx_Res - F2o3*Kx - F8*PI*sx
movy_Res = movy_Res - F2o3*Ky - F8*PI*sy
movz_Res = movz_Res - F2o3*Kz - F8*PI*sz
  endif

#if (ABV == 1)
  call ricci_gamma(ex, X, Y, Z,                                      &
               chi,                                                  &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,&
               Gamx   ,  Gamy    ,  Gamz    , &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,&
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,&
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,&
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,&
               Symmetry)
  call constraint_bssn(ex, X, Y, Z,&
               chi,trK, &
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Gamx,Gamy,Gamz,&
               Lap,betax,betay,betaz,rho,Sx,Sy,Sz,&
               Gamxxx, Gamxxy, Gamxxz,Gamxyy, Gamxyz, Gamxzz, &
               Gamyxx, Gamyxy, Gamyxz,Gamyyy, Gamyyz, Gamyzz, &
               Gamzxx, Gamzxy, Gamzxz,Gamzyy, Gamzyz, Gamzzz, &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz, &
               ham_Res,movx_Res,movy_Res,movz_Res,Gmx_Res,Gmy_Res,Gmz_Res, &
               Symmetry)
#endif 
#if 0
#define i 2
if(Lev == 1)then
write(*,*) X(i),Y(i),Z(i)
write(*,*) Axx(i,i,i),Axy(i,i,i),Axz(i,i,i),Ayy(i,i,i),Ayz(i,i,i),Azz(i,i,i)
write(*,*) 1+Lap(i,i,i),dtSfx(i,i,i),dtSfy(i,i,i),dtSfz(i,i,i)
write(*,*) betax(i,i,i),betay(i,i,i),betaz(i,i,i)
write(*,*) 1+chi(i,i,i),Gamx(i,i,i),Gamy(i,i,i),Gamz(i,i,i)
write(*,*) gxx(i,i,i),gxy(i,i,i),gxz(i,i,i),gyy(i,i,i),gyz(i,i,i),gzz(i,i,i)
write(*,*) trK(i,i,i)
write(*,*) "====="
write(*,*) Axx_rhs(i,i,i),Axy_rhs(i,i,i),Axz_rhs(i,i,i),Ayy_rhs(i,i,i),Ayz_rhs(i,i,i),Azz_rhs(i,i,i)
write(*,*) Lap_rhs(i,i,i),dtSfx_rhs(i,i,i),dtSfy_rhs(i,i,i),dtSfz_rhs(i,i,i)
write(*,*) betax_rhs(i,i,i),betay_rhs(i,i,i),betaz_rhs(i,i,i)
write(*,*) chi_rhs(i,i,i),Gamx_rhs(i,i,i),Gamy_rhs(i,i,i),Gamz_rhs(i,i,i)
write(*,*) gxx_rhs(i,i,i),gxy_rhs(i,i,i),gxz_rhs(i,i,i),gyy_rhs(i,i,i),gyz_rhs(i,i,i),gzz_rhs(i,i,i)
write(*,*) trK_rhs(i,i,i)
endif
#undef i
!!stop
#endif

  gont = 0

  return

  end function compute_rhs_bssn
