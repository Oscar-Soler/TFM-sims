      implicit none


      subroutine initvecpot_usr(ixImin1,ixImin2,ixImax1,ixImax2, ixCmin1,ixCmin2,ixCmax1,ixCmax2, xC, A, idir)
! initialize the vectorpotential on the edges
! used by b_from_vectorpotential()
      integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idir
      double precision, intent(in)       :: xC(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
      double precision, intent(out)      :: A(ixImin1:ixImax1,ixImin2:ixImax2)
      double precision :: rho1(ixImin1:ixImax1,ixImin2:ixImax2),&
                          rho2(ixImin1:ixImax1,ixImin2:ixImax2),&
                          rho3(ixImin1:ixImax1,ixImin2:ixImax2)

      rho1(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sqrt((xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1)-xc1)**2&
                                                +(xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)+hb)**2)
      rho2(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sqrt((xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1)-xc2)**2&
                                                +(xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)+hb2)**2)
      rho3(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sqrt((xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1)-xc3)**2&
                                                +(xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)+hb3)**2)
      if (idir==3) then
         A(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = &
                        -B0*(adipo**2-hb**2)/4.d0/adipo*log((adipo**2+rho1(ixCmin1:ixCmax1,&
                        ixCmin2:ixCmax2)**2-2.d0*adipo*(xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                        2)+hb))/(adipo**2+rho1(ixCmin1:ixCmax1,&
                        ixCmin2:ixCmax2)**2+2.d0*adipo*(xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                        2)+hb)))-B0*(adipo**2-hb2**2)/4.d0/adipo*&
                        log((adipo**2+rho2(ixCmin1:ixCmax1,ixCmin2:ixCmax2)**2&
                        -2.d0*adipo*(xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)+hb2))&
                        /(adipo**2+rho2(ixCmin1:ixCmax1,ixCmin2:ixCmax2)**2&
                        +2.d0*adipo*(xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)+hb2)))&
                        -B0*factorcentralarcade*(adipo**2-hb3**2)/4.d0/adipo&
                        *log((adipo**2+rho3(ixCmin1:ixCmax1,ixCmin2:ixCmax2)**2&
                        -2.d0*adipo*(xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)+hb3))&
                        /(adipo**2+rho3(ixCmin1:ixCmax1,ixCmin2:ixCmax2)**2&
                                +2.d0*adipo*(xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)+hb3)))
      else
         A(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = 0.d0
      end if

      end subroutine initvecpot_usr

      

