	SUBROUTINE NED(AJ,BJ,CJ,AM,BM,CM,CG)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION Q(100,100)
	DOUBLE PRECISION CG
	DOUBLE PRECISION AJ,BJ,CJ,AM,BM,CM
	INTEGER ZZ
	ZZ=MAX(2*AJ+1,2*BJ+1,2*CJ+1,AJ+BJ+CJ,AJ+AM,BJ+BM,CJ+CM)+2
	DO 2 I=1,ZZ
		Q(I,1)=1.D0
		Q(I,I)=1.D0
2	CONTINUE
	DO 3 I=2,ZZ-1
	DO 3 K=2,I
		Q(I+1,K)=Q(I,K-1)+Q(I,K)
3	CONTINUE
	CG=0.D0
	JA=AJ+AM+1.01D0
	MA=AJ-AM+1.01D0
	JB=BJ+BM+1.01D0
	MB=BJ-BM+1.01D0
	JC=CJ+CM+1.01D0
	MC=CJ-CM+1.01D0
	LA=BJ+CJ-AJ+1.01D0
	LB=CJ+AJ-BJ+1.01D0
	LC=AJ+BJ-CJ+1.01D0
	LT=AJ+BJ+CJ+1.01D0
	D=DABS(AM+BM-CM)-0.01D0
	IF (D) 10,10,20
10	LD=MIN0(JA,JB,JC,MA,MB,MC,LA,LB,LC)
	IF (LD) 20,20,30
30	JA2=AJ+AJ+AM+AM
	JB2=BJ+BJ+BM+BM
	JC2=CJ+CJ-CM-CM
	I2=JA2+JB2+JC2-JA2/2*2-JB2/2*2-JC2/2*2
	IF (I2) 20,40,20
40	FN=Q(JA+MA-1,LC)/Q(LT,JC+MC-1)
	FN=FN*Q(JB+MB-1,LC)/Q(LT+1,2)
	FN=FN/Q(JA+MA-1,JA)
	FN=FN/Q(JB+MB-1,JB)
	FN=FN/Q(JC+MC-1,JC)
	K0=MAX(0,LC-JA,LC-MB)+1
	K1=MIN(LC,MA,JB)
	X=0.D0
	DO 50 K=K0,K1
		X=-X-Q(LC,K)*Q(LB,MA-K+1)*Q(LA,JB-K+1)
50	CONTINUE
	IP=K1+LB+JC
	P=1-2*(IP-IP/2*2)
	CG=P*X*DSQRT(FN)
!	What we´ve calculated is a Wigner 3-j coefficient
!	Next, we´ll turn it into a Clebsch-Gordan coefficient
	CG=CG*DSQRT(2*CJ+1)*(-1)**IDNINT(AJ-BJ-CM)
20	CONTINUE
	RETURN
	END






