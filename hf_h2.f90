! simpified diagonalization for sencond order matrix
subroutine QR(F,S,C,e)
real , dimension(2,2) :: F
real , dimension(2,2) :: C
real , dimension(2,2) :: S
real , dimension(2) :: e
real p,q
e(2)=(F(1,1)+F(2,2)+((F(1,1)+F(2,2))**2-4*(F(1,1)*F(2,2)-F(1,2)*F(2,1)))**0.5)/2
e(1)=(F(1,1)+F(2,2)-((F(1,1)+F(2,2))**2-4*(F(1,1)*F(2,2)-F(1,2)*F(2,1)))**0.5)/2
 p=(-F(1,2))/(F(1,1)-e(1))
 q=(-F(1,2))/(F(1,1)-e(2))
 C(1,1)=p*(p**2*S(1,1)+p*S(1,2)+p*S(2,1)+S(2,2))**(-0.5)
 C(2,1)=(p**2*S(1,1)+p*S(1,2)+p*S(2,1)+S(2,2))**(-0.5)
 C(1,2)=q*(q**2*S(1,1)+q*S(1,2)+q*S(2,1)+S(2,2))**(-0.5)
 C(2,2)=(q**2*S(1,1)+q*S(1,2)+q*S(2,1)+S(2,2))**(-0.5)
 RETURN 
 end

! calculation of 2-electron integral
 subroutine pair(A,G,R,ge)
 REAL , DIMENSION(3,2) :: A
 REAL , DIMENSION(3,2) :: G
 REAL , DIMENSION(3,2) :: R
 REAL , DIMENSION(2,2,2,2) :: ge
 integer h,b,c,d,w,s,t,u
 real , dimension(3) :: Rp1
 real , dimension(3) :: Rp2
 real , dimension(3) :: RL
 real RLD
 real pi
 pi=dacos(-1.d0)
 Rp1=0
 Rp2=0
 RL=0
 ge=0
 RLD=0
 do w=1,2
  do s=1,2
   do t=1,2
    do u=1,2
     do h=1,3
      do b=1,3
       do c=1,3
        do d=1,3
        Rp1(1) = (A(h,w) * R(1,w) + A(b,s) * R(1,s)) / (A(h,w) + A(b,s))
        Rp1(2) = (A(h,w) * R(2,w) + A(b,s) * R(2,s)) / (A(h,w) + A(b,s))
        Rp1(3) = (A(h,w) * R(3,w) + A(b,s) * R(3,s)) / (A(h,w) + A(b,s))
        Rp2(1) = (A(c,t) * R(1,t) + A(d,u) * R(1,u)) / (A(c,t) + A(d,u))
        Rp2(2) = (A(c,t) * R(2,t) + A(d,u) * R(2,u)) / (A(c,t) + A(d,u))
        Rp2(3) = (A(c,t) * R(3,t) + A(d,u) * R(3,u)) / (A(c,t) + A(d,u))
        RL=Rp1-Rp2
        RLD=(RL(1)**2+RL(2)**2+RL(3)**2)**0.5
        if (RLD==0) THEN
        ge(w,s,t,u)=ge(w,s,t,u)+G(h,w)*G(b,s)*G(c,t)*G(d,u)*64*(exp(-A(h,w)*A(b,s)*( &
        (R(1,w)-R(1,s))**2+(R(2,w)-R(2,s))**2+(R(3,w)-R(3,s))**2)/(A(h,w)+A(b,s))))*(exp(-A(c,t)* &
        A(d,u)*((R(1,t)-R(1,u))**2+(R(2,t)-R(2,u))**2+(R(3,t)-R(3,u))**2)/(A(c,t)+A(d,u))))*((A(h,w)* &
        A(b,s)*A(c,t)*A(d,u))**0.75*(1/(4*(A(h,w)+A(b,s))*(A(c,t)+A(d,u))))**1.5* &
        (pi*((A(h,w)+A(b,s)+A(c,t)+A(d,u))/(4*(A(h,w)+A(b,s))*(A(c,t)+A(d,u)))))**(-0.5))
        ELSE
       ge(w,s,t,u)=ge(w,s,t,u)+G(h,w)*G(b,s)*G(c,t)*G(d,u)*(64/RLD)*(exp(-A(h,w)*A(b,s)*( &
        (R(1,w)-R(1,s))**2+(R(2,w)-R(2,s))**2+(R(3,w)-R(3,s))**2)/(A(h,w)+A(b,s))))*(exp(-A(c,t)* &
        A(d,u)*((R(1,t)-R(1,u))**2+(R(2,t)-R(2,u))**2+(R(3,t)-R(3,u))**2)/(A(c,t)+A(d,u))))*((A(h,w)* &
        A(b,s)*A(c,t)*A(d,u))**0.75*(1/(4*(A(h,w)+A(b,s))*(A(c,t)+A(d,u))))**1.5)&
        *erf(RLD/(2*(((A(h,w)+A(b,s)+A(c,t)+A(d,u))/(4*(A(h,w)+A(b,s))*(A(c,t)+A(d,u))))**0.5)))
        end if
        end do 
       end do
      end do
     end do
    end do
   end do
  end do 
 end do
 RETURN 
 end

! form density matrix from C
 subroutine density(C,P)
 REAL , DIMENSION(2,2) :: C
 REAL , DIMENSION(2,2) :: P
 integer t,u,j
 P=0
 do t=1,2
  do u=1,2
   do j=1,1
   P(t,u)=P(t,u)+2*C(t,j)*C(u,j)
   end do
  end do
 end  do
 RETURN
 end
 
! main program  
program hf 
 REAL , dimension(2,2) :: T, S, Vn, Hc, Sm, Fm, F , C, P
 real , dimension(2,2,2,2) :: ge
 real , dimension(3,2) :: A, R, G
 real , dimension(2) :: Z, ep, epq
 real pi,E
integer I,M,J,K,L,o,u,x,w

! parameter for sto-3g of H atom
 G(1,1) = 0.444635
 G(2,1) = 0.535328
 G(3,1) = 0.154329
 G(1,2) = 0.444635
 G(2,2) = 0.535328
 G(3,2) = 0.154329

 A(1,1) = 0.168856
 A(2,1) = 0.623913
 A(3,1) = 3.425250
 A(1,2) = 0.168856
 A(2,2) = 0.623913
 A(3,2) = 3.425250

! location of two H atom
 R(1,1)=0
 R(2,1)=0
 R(3,1)=1.4
 R(1,2)=0
 R(2,2)=0
 R(3,2)=0

 S=0
 T=0
 Vn=0
 Z=1
 pi=dacos(-1.D0)
 ep=0
 C=0
 epq=0
 E=0

! calculation of the overlap matrix S, kinetic energy matrix T and nuclear attraction matrix Vn
DO I=1,2
 DO J=1,2
  DO M=1,3
   DO K=1,3
    S(I,J)=S(I,J)+G(M,I)*G(K,J)* &
     ((4*A(M,I)*A(K,J)/((A(M,I)+A(K,J))**2))**0.75)* &
     (EXP(-(A(M,I)*A(K,J)* &
     ((R(1,I)-R(1,J))**2+(R(2,I)-R(2,J))**2+(R(3,I)-R(3,J))**2)) &
     /(A(M,I)+A(K,J))))
    T(I,J)=T(I,J)+G(M,I)*G(K,J)* &
    ((4*A(M,I)*A(K,J)/((A(M,I)+A(K,J))**2))**0.75)* &
     (A(M,I)*A(K,J)/(A(M,I)+A(K,J)))* &
     (EXP(-(A(M,I)*A(K,J)* &
     ((R(1,I)-R(1,J))**2+(R(2,I)-R(2,J))**2+(R(3,I)-R(3,J))**2)) &
     /(A(M,I)+A(K,J))))*(3-(2*A(M,I)*A(K,J)* &
     ((R(1,I)-R(1,J))**2+(R(2,I)-R(2,J))**2+(R(3,I)-R(3,J))**2))/( &
     (A(M,I)+A(K,J))))
     DO L=1,2
      IF( (((((A(M,I)*R(1,I)+A(K,J)*R(1,J))/(A(M,I)+A(K,J)))-R(1,L))**2+ &
       (((A(M,I)*R(2,I)+A(K,J)*R(2,J))/(A(M,I)+A(K,J)))-R(2,L))**2+ &
       (((A(M,I)*R(3,I)+A(K,J)*R(3,J))/(A(M,I)+A(K,J)))-R(3,L))**2)**0.5)==0) THEN
      Vn(I,J)=Vn(I,J)+G(M,I)*G(K,J)* &
       (-2*Z(L)*(EXP(-(A(M,I)*A(K,J)* &
       ((R(1,I)-R(1,J))**2+(R(2,I)-R(2,J))**2+(R(3,I)-R(3,J))**2)) &
       /(A(M,I)+A(K,J)))))*((4*A(M,I)*A(k,J))**0.75)/( &
       (A(M,I)+A(K,J))*(pi**0.5))   
      ELSE 
      Vn(I,J)=Vn(I,J)+G(M,I)*G(K,J)* &
      (-(4*A(M,I)*A(K,J)/((A(M,I)+A(K,J))**2))**0.75)*((Z(L)* &
      (EXP(-(A(M,I)*A(K,J)* &
       ((R(1,I)-R(1,J))**2+(R(2,I)-R(2,J))**2+(R(3,I)-R(3,J))**2)) &
       /(A(M,I)+A(K,J)))))/ &
       ((((A(M,I)*R(1,I)+A(K,J)*R(1,J))/(A(M,I)+A(K,J)))-R(1,L))**2+ &
       (((A(M,I)*R(2,I)+A(K,J)*R(2,J))/(A(M,I)+A(K,J)))-R(2,L))**2+ &
       (((A(M,I)*R(3,I)+A(K,J)*R(3,J))/(A(M,I)+A(K,J)))-R(3,L))**2)**0.5)* &
       ERF(((A(M,I)+A(K,J))**0.5)* &
       (((((A(M,I)*R(1,I)+A(K,J)*R(1,J))/(A(M,I)+A(K,J)))-R(1,L))**2+ &
       (((A(M,I)*R(2,I)+A(K,J)*R(2,J))/(A(M,I)+A(K,J)))-R(2,L))**2+ &
       (((A(M,I)*R(3,I)+A(K,J)*R(3,J))/(A(M,I)+A(K,J)))-R(3,L))**2)**0.5))
      END IF
     END DO
   END DO
  END DO
 END DO
END DO

! inverse matrix of the overlap matrix
Sm(1,1)=S(2,2)/((S(1,1)*S(2,2))-(S(1,2)*S(2,1)))
Sm(2,2)=S(1,1)/((S(1,1)*S(2,2))-(S(1,2)*S(2,1)))
Sm(1,2)=-S(1,2)/((S(1,1)*S(2,2))-(S(1,2)*S(2,1)))
Sm(2,1)=-S(2,1)/((S(1,1)*S(2,2))-(S(1,2)*S(2,1)))

! define Hcoll
Hc=T+Vn

! use the Hcoll as the first guess
F=Hc

Fm=MATMUL(Sm,Hc)

call QR(Fm,S,C,ep)

call pair(A,G,R,ge)

CALL density(C,P)

!SCF 
DO WHILE( abs(epq(1)-ep(1))>0.000001 .or. abs(epq(2)-ep(2))>0.000001 )
epq=ep
F=Hc
 do o=1,2
  do u=1,2
   do x=1,2
    do w=1,2
    F(o,u)=F(o,u)+P(x,w)*(ge(o,u,x,w)-0.5*ge(o,w,x,u))
    end do
   end do
  end do
 end do
Fm=MATMUL(Sm,F)
call QR(Fm,S,C,ep)
call density(C,P)
end DO
E=0

! whole energy
do o=1,2
 do u=1,2
  E=E+0.5*P(o,u)*(Hc(o,u)+F(o,u))
  end do
 end do 

! don’t forget nuclear rejection under the B-O approximation 
E=E+1/R(3,1)
  
! Output results
print*,"Combination coefficient of bonding orbitals ",C(1,1),C(2,1)
print*,"Orbital energy of bonding orbitals(hartree) ",ep(1)
print*,"Combination coefficient of anti-bonding orbitals ",C(1,2),C(2,2)
print*,"Orbital energy of anti-bonding orbitals(hartree) ",ep(2)
print*,"total energy(hartree) ",E
END program



