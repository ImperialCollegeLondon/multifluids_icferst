!eps=0.37
!St_gl = (k_l/d_b)*(7-10*eps+5*eps**2)*(1+0.7*Re_gl**0.2*Pr_l**0.3333)+(k_l/d_b)*(1.33-2.4*eps+1.2*eps**2)*Re_gl**0.7*Pr_l**0.3333 !+ 1.0e+7 * a_lg**10
!St_sl= (k_l/d_p)*(7-10*eps+5*eps**2)*(1+0.7*Re_sl**0.2*Pr_l**0.3333)+(k_l/d_b)*(1.33-2.4*eps+1.2*eps**2)*Re_sl**0.7*Pr_l**0.3333
!St_sg = (k_g/d_p)*(7-10*eps+5*eps**2)*(1+0.7*Re_sg**0.2*Pr_g**0.3333)+(k_l/d_b)*(1.33-2.4*eps+1.2*eps**2)*Re_sg**0.7*Pr_g**0.3333

subroutine vap(p, rho_l, rho_g, Tsat, u_l, u_g, t_l, t_g, a_l, a_g, dtr, Gamma_g,Svap_l,Svap_g,h_l,h_g)
   implicit none

   !f2py intent(in) p, rho_l, rho_g, Tsat, u_l, u_g, t_l, t_g, a_l, a_g,dtr
   !f2py intent(out) Gamma_g,Svap_l,Svap_g,h_l,h_g

   real(kind=8) p, rho_l, rho_g, Tsat, u_l, u_g, t_l, t_g, a_l, a_g
   real(kind=8) F5, Svap_l, Svap_g, Svap_l_max, Svap_g_max, Svap_l1, Svap_l2, rho_v
   real(kind=8) phi, Lh, h_g, h_l, dtr, d_b, C, a_b, a_lg, a_gl, Cp_g, Cp_l, Gamma_g, k_l, Re_lg, mu_l
   real(kind=8), parameter :: Le0 = 2375.7e3
   real(kind=8), parameter :: Pi = 3.14
   mu_l = 3.0e-4
   a_lg = a_l/(a_l + a_g); a_gl = 1.0 - a_lg
   a_b = max(a_g, 1.0e-5)
   Cp_l = 4200.0; Cp_g = 1996.0; k_l = 0.58
   if (Tsat < t_l) then
      h_l = -Le0 + Cp_l*t_l + p/rho_l
   else
      h_l = -Le0 + Cp_l*Tsat + p/rho_l
   end if

   if (Tsat < t_g) then
      h_g = Cp_g*Tsat + p/rho_g
   else
      h_g = Cp_g*t_g + p/rho_g
   end if

   Lh = h_g - h_l
   d_b = 5.0*0.06/max((rho_l*abs(u_g - u_l)**2), 1.0e-3); d_b = max(1.0e-5, d_b)

   Re_lg = rho_l*abs(u_g - u_l)*d_b/mu_l

   !dtr = 5e-2;
   rho_v = 1.0*rho_g
   if (Tsat < t_l) then
      Svap_l_max = ((min(rho_v, a_l*rho_l)/dtr)*Lh)/max(1.0e-10, abs(Tsat - t_l))
   else
      Svap_l_max = (((a_g*rho_v)/dtr)*Lh)/max(1.0e-10, abs(Tsat - t_l))
   end if
   if (Tsat < t_g) then
      Svap_g_max = ((min(rho_v, a_l*rho_l)/dtr)*Lh)/max(1.0e-10, abs(Tsat - t_g))
   else
      Svap_g_max = (((a_g*rho_v)/dtr)*Lh)/max(1.0e-10, abs(Tsat - t_g))
   end if

   if (Tsat < t_l) then
      Svap_l1 = (k_l/d_b)*(12.0/pi)*abs(Tsat - t_l)*(rho_l*Cp_l)/(rho_g*Lh)
      Svap_l2 = (k_l/d_b)*(2.0 + 0.74*(a_l*Re_lg)**0.5)
      Svap_l = max(Svap_l1, Svap_l2)*3.6*a_b/d_b
   else
      if (p <= 1.1272*1.0e6) then
         C = 65.0 - 5.69e-5*(p - 1.0e5)
      else
         C = 2.5e9*p**(-1.418)
      end if

      if (abs(u_g - u_l) <= 0.61) then
         phi = 1.0
      else
         phi = (1.639344*abs(u_g - u_l))**0.47
      end if

      if (a_g < 0.25) then
         F5 = 0.075 + 1.8*phi*C*exp(-45.0*a_b)
      else
         F5 = 0.075
      end if

      Svap_l = F5*Lh*rho_g*rho_l*a_g/(rho_l - rho_g)
      Svap_l = min(Svap_l, 17539.0*max(4.724, 472.4*a_g*a_l)*max(0.0, min(1.0, a_g/0.1)))

   end if

   Svap_g = 1.0e4*3.6*a_b/d_b

   Svap_l = min(Svap_l, Svap_l_max)
   Svap_g = min(Svap_g, Svap_g_max)

   Gamma_g = (Svap_l*(t_l - Tsat) + Svap_g*(t_g - Tsat))/Lh
   return
end subroutine vap



subroutine drag_3_phase(rho_l, rho_g, u_l, u_g, u_s, a_l, a_g, a_s, d_p, d_b, D, S_lg_l, S_lg_g, S_ls_l, S_gs_g)
!this is the drag coefficient for liquid, gas and another solid phase
!input  liquid density         (rho_l)
!       gas density            (rho_l)
!       liquid velocity        (u_l)
!       gas velocity           (u_g)
!       solid velocity         (u_s)
!       liquid volume fraction (a_l)
!       gas volume fraction    (a_g)
!       solid volume fraction  (a_s)
!       solid particle size    (d_p)
!       bubble/droplet size    (d_b)
!       Characteristic length scale (D)

!output drag coefficient
!           liquid    gas     solid
!   liquid            S_lg_l  S_ls_l
!   gas     S_lg_g            S_gs_g
!   solid   none      none    inf      
   implicit none
   !f2py intent(in) rho_l, rho_g, u_l, u_g, u_s, a_l,a_g,a_s,d_b,d_p,D
   !f2py intent(out) S_lg_l, S_lg_g, S_ls_l, S_gs_g

   real(kind=8) rho_l, rho_g, u_l, u_g, u_s, a_l, a_g, a_s, S_lg_l, S_lg_g, S_ls_l, S_gs_g
   real(kind=8) a_sg, a_gs, a_sl, a_ls, a_lg, a_gl, u_gs, u_gl, u_ls, d_b, d_p,D
   real(kind=8) mu_l, mu_g, CD_u, CD_d, A_I, aa_g, w, A_I1, A_I2, Re_gl, Re_lg,Re_sl, CD

   a_sg = a_s/(a_s + a_g); a_gs = 1.0 - a_sg
   a_sl = a_s/(a_s + a_l); a_ls = 1.0 - a_sl
   a_lg = a_l/(a_l + a_g); a_gl = 1.0 - a_lg
   mu_l = 3.0e-4 ! liquid viscosity
   mu_g = 1.0e-5 ! gas viscosity
   u_gs = abs(u_g - u_s); u_ls = abs(u_l - u_s); u_gl = abs(u_g - u_l)

   Re_gl = max(rho_l*u_gl*d_b/mu_l, 1.0e-5); Re_lg = Re_gl
   Re_sl = rho_l*abs(u_s - u_l)*d_p/mu_l
    
   S_ls_l = 150.0 * (a_ls*mu_l) / (a_sl*d_p**2*(a_l+a_s)) + 1.75 * (rho_l*u_ls) / (d_p*(a_l+a_s)) 
   S_gs_g = 150.0 * (a_gs*mu_g) / (a_sg*d_p**2*(a_g+a_s)) + 1.75 * (rho_g*u_gs) / (d_p*(a_g+a_s))

   if (a_gl <= 0.25) then
      CD_d = 1.73205/(1.0 - a_gl)**0.5
      if (Re_gl < 2.0e+5) then
         CD_u = 24.0*(1.0 + 0.1*Re_gl**0.75)/Re_gl
      else if (Re_gl < 5.0e+5) then ! linear interpolation between 2.0e+5 and 5.0e+5
         w = (Re_gl - 2.0e+5)/3.0e+5
         CD_u = (1 - w)*0.11361 + w*0.45078
      else
         CD_u = 0.45078
      end if

      CD = (1.0 - 4.0*a_gl)*CD_u + 4.0*a_gl*CD_d

      !A_I = 6.0*a_gl/d_b
      !Sigma = 0.125 * CD * A_I * rho_l * u_gl

      S_lg_l = 0.125*CD*(6.0*a_gl/d_b)/a_l*rho_l*u_gl
      S_lg_g = 0.125*CD*(6.0/d_b)/(a_g + a_l)*rho_l*u_gl

   else if (a_gl > 0.6) then

      aa_g = 0.3929 - 0.5714*a_gl
      A_I = 4.5/(0.95*D)*(a_gl - 0.05) + ((6.0 - 0.05*aa_g)/d_b)*a_lg/0.95

      !CD = (16.0/15.0)*a_lg
      !Sigma = 0.125 * CD * A_I * rho_l * u_gl

      S_lg_l = 0.125*(16.0/15.0)/(a_l + a_g)*A_I*rho_l*u_gl
      S_lg_g = 0.125*(16.0/15.0)*a_lg/a_g*A_I*rho_l*u_gl

   else

      if (a_gl > 0.55) then
         ! linear interpolation for A_I between 0.55 and 0.6 a_gl

         !A_I for 0.25<a_gl<0.6
         !A_I1 = 4.5/D*(a_gl-aa_g)/(1.0-aa_g) + (6.0*aa_g/d_b)*a_lg/(1.0-aa_g)
         aa_g = 0.3929 - 0.5714*0.55 ! aa_g for a_gl=0.55
         A_I1 = 4.5/D*(0.55 - aa_g)/(1.0 - aa_g) + (6.0*aa_g/d_b)*0.45/(1.0 - aa_g)

         !A_I for 0.6<a_gl
         !A_I2 = 4.5/(0.95*D)*(a_gl-0.05) + ((6.0-0.05*aa_g)/d_b)*a_lg/0.95
         aa_g = 0.3929 - 0.5714*0.6 ! aa_g for a_gl=0.6
         A_I2 = 4.5/(0.95*D)*(0.6 - 0.05) + ((6.0 - 0.05*aa_g)/d_b)*0.4/0.95

         w = (a_gl - 0.55)/0.05
         A_I = (1 - w)*A_I1 + w*A_I2
      else
         aa_g = 0.3929 - 0.5714*a_gl
         A_I = 4.5/D*(a_gl - aa_g)/(1.0 - aa_g) + (6.0*aa_g/d_b)*a_lg/(1.0 - aa_g)
      end if

      CD = (8.0/3.0)*(1.0 - a_gl)**2

      S_lg_l = 0.125*CD*A_I*rho_l*u_gl/a_lg
      S_lg_g = 0.125*CD*A_I*rho_l*u_gl/a_gl

   end if

   return
end subroutine drag_3_phase


subroutine drag_3_phase2(rho_l, rho_g, u_l, u_g, u_s, a_l, a_g, a_s, d_p, D, S_lg_l, S_lg_g, S_ls_l, S_gs_g)
!this is the drag correlation for air/water/porous media
!input  liquid density         (rho_l)
!       gas density            (rho_l)
!       liquid velocity        (u_l)
!       gas velocity           (u_g)
!       solid velocity         (u_s)
!       liquid volume fraction (a_l)
!       gas volume fraction    (a_g)
!       solid volume fraction  (a_s)
!       solid particle size    (d_p)
!       bubble/droplet size    (d_b)
!       Characteristic length scale (D)

!output drag coefficient
!           liquid    gas     solid
!   liquid            S_lg_l  S_ls_l
!   gas     S_lg_g            S_gs_g
!   solid   none      none    inf      
   implicit none
   !f2py intent(in) rho_l, rho_g, u_l, u_g, u_s, a_l,a_g,a_s,d_b,d_p,D
   !f2py intent(out) S_lg_l, S_lg_g, S_ls_l, S_gs_g

   real(kind=8) rho_l, rho_g, u_l, u_g, u_s, a_l, a_g, a_s, S_lg_l, S_lg_g, S_ls_l, S_gs_g
   real(kind=8) a_sg, a_gs, a_sl, a_ls, a_lg, a_gl, u_gs, u_gl, u_ls, d_b, d_p,D
   real(kind=8) mu_l, mu_g, CD_u, CD_d, A_I, aa_g, w, A_I1, A_I2, Re_gl, Re_lg,Re_sl, CD

   a_sg = a_s/(a_s + a_g); a_gs = 1.0 - a_sg
   a_sl = a_s/(a_s + a_l); a_ls = 1.0 - a_sl
   a_lg = a_l/(a_l + a_g); a_gl = 1.0 - a_lg
   mu_l = 3.0e-4 ! liquid viscosity
   mu_g = 1.0e-5 ! gas viscosity
   u_gs = abs(u_g - u_s); u_ls = abs(u_l - u_s); u_gl = abs(u_g - u_l)

   
   !Ergun equation
   S_ls_l = 150.0 * (a_ls*mu_l) / (a_sl*d_p**2*(a_l+a_s)) + 1.75 * (rho_l*u_ls) / (d_p*(a_l+a_s)) 
   S_gs_g = 150.0 * (a_gs*mu_g) / (a_sg*d_p**2*(a_g+a_s)) + 1.75 * (rho_g*u_gs) / (d_p*(a_g+a_s))

   !Chris equation
   S_lg_l = D*a_g        !D*a_g*a_l/a_l
   S_lg_g = D*a_l        !D*a_g*a_l/a_g


   return
end subroutine drag_3_phase2


subroutine tabsorb(rho_l,rho_g,u_l,u_g,u_s,Cp_l,Cp_g,mu_l,mu_g,k_l,k_g,d_p,d_b,St_gl,St_sl,St_sg)
   implicit none

   !f2py intent(in) rho_l,rho_g,u_l,u_g,u_s,Cp_l,Cp_g,mu_l,mu_g,k_l,k_g,d_p,d_b
   !f2py intent(out) St_gl,St_sl,St_sg
   !Cp_l: heat capacity of water 4200
   !Cp_l: heat capacity of steam 1996
   !k_l: heat conductivity of liquid 0.606
   !k_g: heat conductivity of steam 0.0184
   
   real(kind=8) rho_l,rho_g,u_l,u_g,u_s,Cp_l,Cp_g,mu_l,mu_g,k_l,k_g,d_p,d_b,St_gl,St_sl,St_sg
   real(kind=8) Re_sl,Re_sg,Re_gl,Pr_l,Pr_g,Nu_sl,Nu_sg,Nu_gl

   Re_sl = rho_l*abs(u_s-u_l)*d_p/mu_l
   Re_sg = rho_g*abs(u_s-u_g)*d_p/mu_g
   Re_gl = rho_l*abs(u_g-u_l)*d_b/mu_l
      
   Pr_l = Cp_l*mu_l/k_l
   Pr_g = Cp_g*mu_g/k_g
   
   Nu_sl=2.0+0.6*Re_sl**0.5*Pr_l**0.3333
   Nu_sg=2.0+0.6*Re_sg**0.5*Pr_g**0.3333
   Nu_gl=2.0+0.6*Re_gl**0.5*Pr_l**0.3333

   St_sl = (k_l/d_p)*Nu_sl
   St_sg = (k_g/d_p)*Nu_sg
   St_gl = (k_l/d_b)*Nu_gl   

   return
   end subroutine tabsorb
