program test

   implicit none

   double complex, dimension (21,21) :: s1,s2,s3,s4,s5,s6,psi4,&
          Jb,&
          J_u,J,ethb_J,J_l,eth_J,&
          sigJ,sigJb,sigu,sigr,sigrb,sigK,sigru,siguu,sigub,&
          ethb_sigJ,sigJ_u,&
          eth_sigK,ethb_sigK,sigK_u,&
          ethb_sigu,eth_sigu,&
          eth_beta,&
          U,ethb_U,U_l,eth_U,&
          Ub, &
          eth_K
   double precision, dimension (21,21) :: K,K_u,K_l,e2beta

Jb = conjg(J)
Ub = conjg(U)

!> proc_fort(psi4):
      s5 = -J*J*K*conjg(J_u)*sigK*e2beta/2.0+J*Ub*eth_beta*&
sigK*K*e2beta/4.0+Ub*J*J*K*conjg(eth_J)*sigJ*e2beta/8.0+Ub*J*K*K*&
eth_J*sigJb*e2beta/8.0-U*J*eth_K*sigJb*K*K*e2beta/4.0+Ub*J*J*K*conjg(&
eth_K)*sigK*e2beta/4.0-U*Jb*eth_J*sigK*K*K*e2beta/8.0+K*K*Ub*eth_K*&
sigJ*e2beta/8.0+Ub*J*K*ethb_sigJ*e2beta/4.0+J*K*K*K_u*sigK*e2beta+U*J&
*K*conjg(ethb_sigJ)*e2beta/4.0-Ub*J*K*conjg(ethb_J)*sigJ*e2beta/8.0+J*U*&
K*U_l*sigJb/8.0+U*J*K*K*conjg(eth_J)*sigJ*e2beta/8.0+J*U*eth_beta*&
sigJb*K*e2beta/4.0+J*U*conjg(eth_K)*sigK*K*K*e2beta/4.0+U*J*K_l*sigub&
*K/2.0+U*conjg(eth_K)*sigJ*e2beta/8.0+U*K*K*ethb_sigJ*e2beta/4.0-3.0/4.0*&
K*K*J_u*sigK*e2beta
      s4 = s5+U*U*J*conjg(eth_K)*sigrb/8.0-J*J*U*&
conjg(J_u)*sigrb/8.0-J_l*Ub*Ub*sigJ*K*K/8.0+Ub*J*K_u*sigr/4.0-sigru*&
J*conjg(ethb_U)*K/4.0+U*U*J*K_l*sigJb/4.0-K*Ub*U_l*sigJ/4.0-U*J*&
conjg(U_l)*sigK/8.0-J*eth_K*sigub*e2beta/4.0+U*J*ethb_J*sigJb*K*e2beta/&
8.0-U*J*eth_K*sigJb*K*e2beta/2.0-U*J_l*Ub*sigK*K*K/4.0-Ub*J*K*K*&
ethb_J*sigK*e2beta/8.0+J*U*conjg(eth_beta)*sigK*K*e2beta/4.0+U*Jb*ethb_J&
*sigJ*K*K*e2beta/8.0-U*Jb*eth_J*sigK*K*e2beta/4.0+U*Jb*ethb_J*sigJ*K&
*e2beta/4.0-J*Ub*conjg(eth_K)*sigJ*K*K*e2beta/4.0+J*K*conjg(J_u)*sigJ*&
e2beta/4.0+Jb*J_u*sigJ*K*K*e2beta/4.0+J*U*conjg(eth_K)*sigK*e2beta/4.0
      s5 = Ub*J*K*K*eth_K*sigK*e2beta/4.0-J*sigK_u*K*e2beta+K&
*eth_K*sigu*e2beta/4.0-U*eth_beta*sigK*e2beta/4.0-J*J*J*conjg(eth_U)*&
sigJb*e2beta/8.0+J*eth_beta*sigub*e2beta/4.0+K*eth_J*sigub*e2beta/8.0-U*&
K*K*eth_U*sigrb/8.0+Ub*J*K_l*sigu/2.0+J*sigru*K_u*K/2.0-J*Ub*conjg(&
U_l)*sigJ/8.0-U*U*ethb_J*sigrb*K*K/16.0-sigru*Ub*eth_J*K/4.0+J*J*J*&
conjg(J_u)*sigJb*e2beta/4.0-eth_sigu*e2beta/4.0+s4+J*K*K*conjg(&
J_u)*sigJ*e2beta/4.0+J*J*K*conjg(eth_K)*sigub*e2beta/4.0+K*Ub*eth_K*&
sigJ*e2beta/4.0-Jb*eth_J*sigu*e2beta/8.0-Ub*J_u*sigr*K*K/8.0
      s3 = s5-sigru*U*ethb_J*K*K/8.0-Ub*J*J*conjg(&
J_l)*sigu/4.0-U*Ub*eth_J*sigrb/16.0-Ub*J*eth_sigK*e2beta/4.0-U*J*&
ethb_sigK*e2beta/4.0-U*K*K*eth_sigK*e2beta/4.0-Ub*J*J*ethb_sigK*e2beta/&
4.0-Ub*ethb_J*sigJ*e2beta/4.0-J*J*U*U*conjg(eth_J)*sigrb/16.0+Ub*Ub*J&
*K_l*sigJ/4.0-U*J*conjg(ethb_U)*sigrb/8.0-eth_beta*sigu*e2beta/4.0-U*J*&
K*K*conjg(ethb_J)*sigK*e2beta/8.0+J*J*conjg(eth_U)*sigK*K*e2beta/4.0-U*&
J*J*conjg(eth_J)*sigK*e2beta/8.0+U*J*ethb_J*sigJb*e2beta/8.0-J*K*conjg(&
ethb_J)*sigu*e2beta/4.0+U*K*eth_J*sigJb*e2beta/8.0-J*K*conjg(ethb_U)*&
sigK*e2beta/4.0-J*K*K*conjg(ethb_U)*sigK*e2beta/4.0+K*K*K*Ub*ethb_J*&
sigJ*e2beta/8.0
      s5 = -J*J*U*conjg(eth_beta)*sigJb*e2beta/4.0-J*J*U*conjg&
(eth_K)*sigJb*e2beta/8.0-Jb*eth_J*sigu*K*e2beta/4.0-J*K*U*conjg(U_l)*&
sigK/8.0-U*J*K*ethb_sigK*e2beta/4.0+U*J*Ub*eth_K*sigrb/8.0+J*J*K*&
conjg(ethb_U)*sigJb*e2beta/8.0+Ub*Ub*J*K_l*sigJ*K/4.0-siguu*J_l/4.0-U*&
Ub*ethb_J*sigr*K*K/16.0-Ub*J_l*sigu*K/2.0-Ub*J_l*sigu*K*K/4.0-J*J*&
sigru*U*conjg(eth_J)/8.0-Ub*K*eth_U*sigr/4.0+J*K*ethb_sigu*e2beta/4.0+K&
*K*conjg(ethb_U)*sigJ*e2beta/4.0+J*J*conjg(eth_U)*sigK*e2beta/4.0-U*&
ethb_J*sigK*e2beta/8.0-J*J*conjg(eth_beta)*sigub*e2beta/4.0-K*K*K*J_u*&
sigK*e2beta/2.0+J*siguu*K_l*K/2.0
      s4 = s5-Ub*Ub*eth_J*sigr*K*K/16.0-U*U*J_l*&
sigJb*K*K/8.0+Ub*J*J*conjg(ethb_sigJ)*e2beta/4.0+U*J*J*conjg(eth_U)*&
sigrb/8.0-U*Ub*eth_J*sigrb*K*K/16.0+sigJ_u*e2beta/2.0-J*ethb_J*sigub*K&
*K*e2beta/8.0-J*K*K*conjg(eth_U)*sigJ*e2beta/8.0-Jb*eth_J*sigu*K*K*&
e2beta/8.0-U*eth_beta*sigK*K*e2beta/2.0+J*K*K*J_u*sigJb*e2beta/4.0-J*J&
*K*ethb_U*sigJb*e2beta/8.0+J*U*conjg(eth_beta)*sigK*e2beta/4.0+Jb*J_u*&
sigJ*K*e2beta/2.0+K*K*K*eth_K*sigu*e2beta/4.0-Ub*eth_beta*sigJ*e2beta/&
4.0+J*conjg(eth_beta)*sigu*e2beta/4.0+Jb*J_u*sigJ*e2beta/4.0-Ub*J_u*sigr*&
K/4.0+J*ethb_J*sigub*e2beta/8.0+K*conjg(ethb_U)*sigJ*e2beta/8.0
      s5 = s4+sigJ_u*K*K*e2beta/2.0-Ub*Ub*J*J*conjg(&
J_l)*sigJ/8.0-J*conjg(U_l)*sigu*K/8.0-Ub*Ub*eth_J*sigr*K/8.0-Ub*K*K*&
eth_U*sigr/8.0-U*J_l*sigub*K*K/4.0-3.0/8.0*U*K*ethb_J*sigK*e2beta+Ub*&
J*U*conjg(eth_K)*sigr/8.0+U*U*J*K_l*sigJb*K/4.0-J*K*K*conjg(ethb_J)*&
sigu*e2beta/8.0-U*J*eth_K*sigJb*e2beta/4.0+Ub*J*J*J*conjg(ethb_J)*&
sigJb*e2beta/8.0+J*Ub*conjg(eth_beta)*sigJ*e2beta/4.0-3.0/8.0*U*K*K*&
ethb_J*sigK*e2beta+Ub*J*U*conjg(eth_K)*sigr*K/8.0+Ub*J*K*eth_J*sigJb*&
e2beta/8.0+U*J*J*K*conjg(ethb_J)*sigJb*e2beta/8.0+J*conjg(eth_beta)*sigu&
*K*e2beta/4.0+J*Ub*eth_beta*sigK*e2beta/4.0-Ub*eth_J*sigK*K*K*K*e2beta&
/8.0
      s2 = s5-U*K*eth_sigK*e2beta/2.0-U*K*K*U_l*sigK/&
8.0+K*K*K*conjg(ethb_U)*sigJ*e2beta/8.0+U*J*K_u*sigrb/4.0+Ub*eth_J*&
sigK*K*e2beta/8.0+Ub*J*J*conjg(ethb_J)*sigK*e2beta/8.0+U*K*K*eth_K*&
sigK*e2beta/2.0-Ub*J*conjg(ethb_U)*sigr*K/8.0+J*K*K_u*sigK*e2beta-3.0/&
8.0*K*K*U*conjg(eth_K)*sigJ*e2beta-Ub*J*J*J*conjg(eth_J)*sigK*e2beta/&
8.0-Ub*eth_J*sigK*K*K*e2beta/8.0+J*J*K*conjg(eth_J)*sigu*e2beta/8.0-Ub&
*eth_beta*sigJ*K*e2beta/2.0+J*K*K*ethb_U*sigK*e2beta/4.0+Jb*K*eth_U*&
sigJ*e2beta/4.0-3.0/8.0*K*Ub*ethb_J*sigJ*e2beta+U*K*K*K*eth_K*sigK*&
e2beta/4.0+Ub*J*K_u*sigr*K/4.0+Ub*J*U*K_l*sigK/2.0-Ub*J*J*conjg(J_u)&
*sigr/8.0+s3
      s5 = U*J*conjg(ethb_sigJ)*e2beta/4.0-Ub*J*conjg(ethb_U)*&
sigr/8.0+s2+J*sigru*Ub*eth_K*K/4.0-Ub*J*J*U*conjg(eth_J)*&
sigr/16.0-U*J_u*sigrb/8.0-Ub*Ub*eth_J*sigr/16.0-Ub*eth_U*sigr/8.0-U*&
U_l*sigK/8.0-sigru*U*ethb_J/8.0+K_u*sigJ*e2beta/4.0+eth_U*sigK*e2beta/4.0-&
U*J_l*sigub/4.0-sigru*J*conjg(ethb_U)/4.0-sigru*K*K*eth_U/4.0-sigru*J_u&
/4.0-U*Ub*eth_J*sigrb*K/8.0+J*K*ethb_U*sigK*e2beta/4.0+J*U*eth_beta*&
sigJb*e2beta/4.0-J*Ub*eth_K*sigK*e2beta/4.0
      s4 = s5+J*K*Ub*U_l*sigK/8.0-J*Ub*K*conjg(U_l)*&
sigJ/8.0-J*eth_K*sigub*K*e2beta/2.0+U*Jb*ethb_J*sigJ*e2beta/8.0+J*Ub*&
ethb_J*sigK*e2beta/8.0+J*K*J_u*sigJb*e2beta/4.0+Ub*J*K*ethb_U*sigr/8.0-&
Ub*J*J*U*conjg(J_l)*sigK/4.0-U*J_l*Ub*sigK/4.0+J*sigru*U*conjg(eth_K&
)/4.0-J*J*sigru*conjg(J_u)/4.0-siguu*J_l*K/2.0+J*J*conjg(U_l)*sigub/8.0&
+J_u*sigK*e2beta/4.0-J*J*conjg(eth_sigu)*e2beta/4.0-K*eth_sigu*e2beta/2.0+&
J*ethb_sigu*e2beta/4.0-K*K*eth_sigu*e2beta/4.0+J*J*conjg(sigJ_u)*e2beta/&
2.0-J*sigK_u*e2beta+J*conjg(ethb_sigu)*e2beta/4.0
      s5 = s4-U*eth_sigK*e2beta/4.0+U*ethb_sigJ*e2beta/&
4.0+sigJ_u*K*e2beta-Ub*J_l*sigu/4.0+J*J*sigru*conjg(eth_U)/4.0-U*U*J_l*&
sigJb/8.0+U*J*Ub*eth_K*sigrb*K/8.0-U*J_l*Ub*sigK*K/2.0+U*K*K*K*&
eth_J*sigJb*e2beta/8.0+U*J*J*conjg(ethb_J)*sigJb*e2beta/4.0-U*eth_beta*&
sigK*K*K*e2beta/4.0-J*eth_K*sigub*K*K*e2beta/4.0-U*Ub*ethb_J*sigr/&
16.0-eth_beta*sigu*K*e2beta/2.0+Jb*eth_U*sigJ*e2beta/8.0-eth_beta*sigu*K*&
K*e2beta/4.0+J*J*conjg(ethb_J)*sigub*e2beta/4.0-K*K*K*ethb_U*sigJ*&
e2beta/8.0-K*K*K*K_u*sigJ*e2beta/2.0-J*J*conjg(J_u)*sigK*e2beta/4.0
      s3 = s5+K*K*K*eth_J*sigub*e2beta/8.0+Ub*Ub*J*&
eth_K*sigr/8.0-J*J*sigru*Ub*conjg(ethb_J)/8.0+J*sigru*Ub*eth_K/4.0-&
sigru*U*ethb_J*K/4.0+J*K*conjg(ethb_sigu)*e2beta/4.0+J*K*U_l*sigub/8.0-&
U*J_u*sigrb*K*K/8.0-U*U*J_l*sigJb*K/4.0-U*K*U_l*sigK/4.0+J*J*Ub*&
conjg(U_l)*sigK/8.0-J*J*conjg(J_l)*U*U*sigJb/8.0+J*J*U*conjg(U_l)*&
sigJb/8.0+J*Ub*U_l*sigK/8.0+J*Ub*ethb_U*sigr/8.0-U*K*eth_U*sigrb/4.0+&
J*U*ethb_U*sigrb/8.0+Ub*eth_J*sigK*e2beta/8.0-eth_U*sigK*K*K*K*e2beta/&
4.0+U*J*K*ethb_U*sigrb/8.0-K*K*U_l*sigu/8.0-J_l*Ub*Ub*sigJ/8.0
      s5 = s3-J*conjg(U_l)*sigu/8.0-U*U*ethb_J*sigrb/&
16.0-Ub*U_l*sigJ/8.0-Ub*J_u*sigr/8.0-K*U_l*sigu/4.0+J*sigru*K_u/2.0-U*&
eth_U*sigrb/8.0-sigru*J_u*K*K/4.0-siguu*J_l*K*K/4.0-ethb_J*sigu*e2beta/&
4.0-J*J*siguu*conjg(J_l)/4.0+J*siguu*K_l/2.0+J*sigru*ethb_U/4.0+J*U_l*&
sigub/8.0-sigru*J_u*K/2.0-sigru*Ub*eth_J/8.0-sigru*K*eth_U/2.0+J*conjg(&
eth_K)*sigu*e2beta/4.0-J*eth_U*sigJb*e2beta/8.0-J*conjg(eth_U)*sigJ*e2beta&
/8.0
      s4 = s5-3.0/8.0*K*ethb_J*sigu*e2beta-K*K*ethb_U*&
sigJ*e2beta/4.0-J*J*U*conjg(J_l)*sigub/4.0-K*K*Ub*U_l*sigJ/8.0+Ub*J*&
J*conjg(eth_U)*sigr/8.0-sigru*Ub*eth_J*K*K/8.0-U*U*ethb_J*sigrb*K/8.0&
+J*U*U_l*sigJb/8.0-U*J_u*sigrb*K/4.0-J_l*Ub*Ub*sigJ*K/4.0-U*J_l*&
sigub*K/2.0+J*sigru*K*ethb_U/4.0+U*J*K_l*sigub/2.0-Ub*Ub*J*J*conjg(&
ethb_J)*sigr/16.0-Ub*eth_beta*sigJ*K*K*e2beta/4.0-U*K*K*K*ethb_J*sigK&
*e2beta/8.0-Ub*J*J*eth_K*sigJb*e2beta/8.0+U*J*K_u*sigrb*K/4.0-Ub*J*K&
*eth_sigK*e2beta/4.0-sigru*eth_U/4.0+J*Ub*conjg(eth_K)*sigJ*e2beta/4.0
      s5 = s4+U*U*J*conjg(eth_K)*sigrb*K/8.0+J*sigru&
*U*conjg(eth_K)*K/4.0-U*J*K*conjg(ethb_U)*sigrb/8.0+Ub*eth_K*sigJ*&
e2beta/8.0+eth_U*sigK*K*e2beta/4.0-K*ethb_U*sigJ*e2beta/8.0-3.0/4.0*K*K*&
K_u*sigJ*e2beta-J*J*K_u*sigJb*e2beta/4.0+K*K*eth_K*sigu*e2beta/2.0-J*J&
*J*conjg(eth_J)*sigub*e2beta/8.0+Ub*J*ethb_sigJ*e2beta/4.0+U*K*ethb_sigJ*&
e2beta/2.0+K*K*eth_J*sigub*e2beta/4.0+K*K*K*ethb_J*sigu*e2beta/8.0-J*&
conjg(ethb_J)*sigu*e2beta/8.0-U_l*sigu/8.0+J*J*K*conjg(ethb_J)*sigub*&
e2beta/8.0-U*Jb*eth_J*sigK*e2beta/8.0-K*K*K*U*conjg(eth_K)*sigJ*e2beta/&
4.0-J*J*Ub*conjg(eth_beta)*sigK*e2beta/4.0+Ub*J*K_l*sigu*K/2.0
      s6 = s5-U*Ub*ethb_J*sigr*K/8.0+U*K*eth_K*sigK*&
e2beta/4.0-J*K*conjg(eth_U)*sigJ*e2beta/4.0+J*eth_beta*sigub*K*e2beta/4.0&
+Jb*K*K*eth_U*sigJ*e2beta/8.0-J*conjg(eth_K)*sigu*K*K*e2beta/4.0+J*&
eth_U*sigJb*K*K*e2beta/8.0-J*U*conjg(ethb_J)*sigK*e2beta/8.0-J*Ub*conjg&
(ethb_J)*sigJ*e2beta/8.0+Ub*J*U*K_l*sigK*K/2.0
      s1 = s6-J*J*U*Ub*conjg(ethb_J)*sigrb/16.0+Ub*&
Ub*J*eth_K*sigr*K/8.0-U*J*J*conjg(eth_J)*sigK*K*e2beta/8.0+J*U*conjg&
(eth_K)*sigK*K*e2beta/2.0-J*J*K*K_u*sigJb*e2beta/2.0+U*K*K*eth_J*&
sigJb*e2beta/4.0-J*J*K*Ub*eth_K*sigJb*e2beta/4.0+U*J*K*conjg(eth_J)*&
sigJ*e2beta/8.0-U*J*K*conjg(ethb_J)*sigK*e2beta/4.0-Ub*J*J*K*conjg(&
ethb_J)*sigK*e2beta/8.0+J*Ub*conjg(eth_beta)*sigJ*K*e2beta/4.0-eth_U*sigK&
*K*K*e2beta/4.0
      s2 = 1/(K+1.0)/(e2beta*e2beta)
      psi4 = s1*s2

end program test
