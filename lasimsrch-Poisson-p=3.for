	use msimsl
      IMPLICIT double precision(A-H,O-Z)
      doubleprecision XT(750),YT(750),WT(750),YO(30),R(30),xo(30),   
     cX(750),Y(750),w(750),indexa(750),indexa2(750),inew(30),
     c       wmin,denom,wbar(10000),wbar2(10000)
	real name(750),mu,lambda
	integer ir(1),itemp(1)

	OPEN(7,FILE='C:\pmedian\p-median-p=2-wsa=5.txt') 
	OPEN(8,FILE='C:\pmedian\p-median-allocations-p=2-wsa=5.txt') 
	OPEN(9,FILE='C:\pmedian\p-median-mu&var-p=2-wsa=5.txt') 

	wsa = 5.0
	irep=1
      isim=10000

	m=2
	mu=4.0
	p=20
	q=20
	n=p*q
	wmin=1.0d10
	k=0
	denom=0.0d0
	do 400 i=1,p   
	do 399 j=1,q 
	k=k+1
	name(k)=k    
	x(k)=real(i)
	y(k)=real(j)  
      denom=denom+dexp(wsa*(dsin(dacos(-1.0d0)*2.0*x(k)/real(p+1))*
     c          dsin(dacos(-1.0d0)*y(k)/real(q+1)) + 
     c          dsin(dacos(-1.0d0)*x(k)/real(p+1))*
     c		  dsin(dacos(-1.0d0)*2.0*y(k)/real(q+1)))/2.0)
  399 continue                                                
  400	continue
	write(6,*) denom
	pause

	do 998 krep=1,irep
	write(6,*) krep
	wbar(krep)=0.0d0
	wbar2(krep)=0.0d0
	do 340 i=1,n
	lambda=mu*real(n)*
     c  (dexp( wsa*(dsin(dacos(-1.0d0)*2.0*x(i)/real(p+1))*
     c              dsin(dacos(-1.0d0)*y(i)/real(q+1)) + 
     c              dsin(dacos(-1.0d0)*x(i)/real(p+1))*
     c		      dsin(dacos(-1.0d0)*2.0*y(i)/real(q+1)))/2.0 ))/denom 

	call rnpoi(1,lambda,ir)		 
	w(i)=real(ir(1)) + 1.0
	wbar(krep)=wbar(krep)+w(i)
	wbar2(krep)=wbar2(krep)+w(i)**2
 
c	w(i)=1.0

  340 CONTINUE 

	wbar(krep)=wbar(krep)/real(n) 
	wbar2(krep)=wbar2(krep)/real(n-1) - 
     c                             real(n)*(wbar(krep)**2)/real(n-1)
	write(9,*) krep,wbar(krep),wbar2(krep)

	globald=1.0e25
C   
C     ORDER THE DESTINATIONS BY THEIR X CO-ORDINATES:       
   43  NM1=N-1       
       DO 10 I=1,NM1 
       IP1=I+1       
       DO 9 J=IP1,N  
       if (X(I).LE.X(J)) GO TO 9     
       T=X(I)
       X(I)=X(J)     
       X(J)=T
       T=Y(I)
       Y(I) = Y(J)   
       Y(J)=T
       T=W(I)   
       W(I)=W(J)
       W(J)=T   
       IT=indexa(I)      
       indexa(I)=indexa(J)
       indexa(J)=IT
	 nametmp=name(i)
	 name(i)=name(j)
	 name(j)=nametmp      
    9 CONTINUE  
 10    CONTINUE 
      write(8,200) N,m,isim
 200  FORMAT (//' NUMBER OF DESTINATIONS:',I4/ 
     * ' NUMBER OF SOURCES:',I9/,' NUMBER OF RANDOM EVALUATIONS:'I9/       
     * ' INPUT DATA ORDERED ASCENDING ON X:'/' DESTINATION',6X, 
     * 'CO-ORDINATES',12X,'WEIGHT   ALLOCATION'///)   
      DO 35 I=1,N       
c      write(8,201) X(I),Y(I),W(I),indexa(I),name(i)    
  201 FORMAT(10x,5f10.0)   
   35 CONTINUE  
      if (indexa(1).NE.0) GO TO 12       
C       
C    NO INITIAL ALLOCATION IS GIVEN--GENERATE ONE:      
      KK = 0    
      DO 11 I=1,N       
      if(I.GT.M) GO TO 11       
      KK = KK + 1       
 11   indexa(I)=KK       
   12 globald = 1.0d25


      do 999 irun=1,isim
      print *,'RUN NUMBER ',irun
      icount = 0
      RMIN=1.0D25       
C       
C    USE ALTERNATE ROUTINE TO ACHIEVE A STABLE PARTITION:       
   51 icount = icount + 1
      dO 14 J=1,M       
      NN=0      
      DO 13 I=1,N       
      if (indexa(I).NE.J) GO TO 13       
      NN=NN+1   
      XT(NN)=X(I)       
      YT(NN)=Y(I)       
      WT(NN)=W(I)       
 13   CONTINUE  
 14   CALL WEBER (NN,XT,YT,WT,XO(J),YO(J),R(J)) 
      if(indexa(1).EQ.1) GO TO 50
      TEMPX = XO(indexa(1))      
      XO(indexa(1)) = XO(1)      
      XO(1) = TEMPX     
      TEMPY = YO(indexa(1))      
      YO(indexa(1)) = YO(1)      
      YO(1) = TEMPY     
 50   DO 16 I=1,N       
      RR=1.0D25 
      DO 15 J=1,M       
      D=(X(I)-XO(J))**2+(Y(I)-YO(J))**2 
      if (D.GE.RR) GO TO 15     
      RR=D      
      iiset=J    
 15   CONTINUE  
 16   indexa(I)=iiset     
      D=0.0     
      DO 17 J=1,M       
 17   D=D+R(J)  
      if(D.GE.RMIN) GO TO 18    
      RMIN=D    
      GO TO 51
	 
   18 continue
c
c record all
c
      if(d.gt.globald) go to 996
	if(xo(1).eq.xolag1.and.yo(1).eq.yolag1) goto 996
	xolag1=xo(1)
	yolag1=yo(1)

      globald = d
c      write(7,1000) irun
 1000 format(///,1x,'RUN NUMBER ',i5)
c      write(7,300) (i,xo(i),yo(i),r(i),i=1,m)
      write(7,3000) krep,irun,icount,d,(xo(i),yo(i),i=1,m)
 3000 format(3(1x,i9),2x,f25.10,10f10.6)
	do 900 i=1,n
  900 indexa2(i)=indexa(i)   
  996 if(isim.eq.0) go to 999
      call rnset(0)
      k = 0
      do 937 i=1,m
  937 inew(i) = 0
      do 940 i=1,m
  938 call rnund(1,n,itemp)
      do 939 j=1,k
  939 if(inew(j).eq.itemp(1)) go to 938
      k = k + 1
  940 inew(k) = itemp(1)
      do 60 i=1,m
      xo(i) = x(inew(i))
   60 yo(i) = y(inew(i))
      DO 160 I=1,N       
      RR=1.0D25 
      DO 150 J=1,M       
      D=(X(I)-XO(J))**2+(Y(I)-YO(J))**2 
      if (D.GE.RR) GO TO 150     
      RR=D      
      iiset=J    
 150  CONTINUE  
 160  indexa(I)=iiset     
  999 continue

  998 continue
C       
C    END OF ALTERNATE ROUTINE.  
C       
 300  FORMAT (/' GROUP     X-CENTER     Y-CENTER     DISTANCE'/     
     * (I5,1X,3G12.5))       
 301  FORMAT (' NUMBER OF EVALUATIONS:  ',I10,5X,    
     *' AGGREGATE DISTANCE:  ',G12.5,/,' CONFIGURATION:  ',
     *10(/,21X,25I3)) 
      
	do 901 i=1,n
  901 write(8,201) name(i),X(I),Y(I),W(I),indexa2(I)    
  
      STOP   
      END    
c
c
c
      SUBROUTINE WEBER (N,X,Y,W,XO,YO,SR)    
      IMPLICIT REAL*8 (A-H,O-Z)      
      DIMENSION X(N),Y(N),W(N)       
      DATA ITER,TOL /1000,1.0D-20/      
C    
C    THE FIRST ESTIMATE OF THE POINT OF MINIMUN AGGREGATE TRAVEL     
C    IS THE CENTROID:
      XO=0.0d0       
      YO=0.0d0       
      SW=0.0d0       
      DO 4 I=1,N     
      XO=XO+X(I)*W(I)
      YO=YO+Y(I)*W(I)
 4    SW=SW+W(I)     
C     if(SW.EQ.0.0) RETURN      
      XO=XO/SW       
      YO=YO/SW       
      XT=XO  
      YT=YO  
C    
      SR=0.0d0       
      if (N.EQ.1) RETURN     
C    
C     BEGIN ITERATIVE LOOP.  
C    
      RMIN=1.0D25    
c      DO 20 J=1,ITER 

	icount=0
    5 icount=icount+1
	if(icount.gt.iter) goto 21
      FX=0.0d0       
      FY=0.0d0       
      FXX=0.0d0      
      FYY=0.0d0      
      FXY=0.0d0      
      SR=0.0d0       
C    
      DO 10 I=1,N    
      DX=X(I)-XO     
      DY=Y(I)-YO     
      DX2=DX*DX      
      DY2=DY*DY      
      R2=DX2+DY2     
      if (R2.EQ.0.0d0) GO TO 10      
      R=DSQRT(R2)    
      DIV=W(I)/R     
      FX=FX-DX*DIV   
      FY=FY-DY*DIV   
      DIV=DIV/R2     
      FXX=FXX+DY2*DIV
      FYY=FYY+DX2*DIV
      FXY=FXY-DX*DY*DIV      
      SR=SR+R*W(I)   
 10   CONTINUE       
C    
C     END ITERATIVE LOOP     
C    
C    
      if (RMIN.GE.SR) GO TO 15       
C    
C    AGGREGATE DISTANCE HAS INCREASED. USE BINARY CHOP:      
      TX=XO  
      TY=YO  
      XO=0.5D0*(XO+XT)       
      YO=0.5D0*(YO+YT)       
      XT=TX  
      YT=TY  
      GO TO 20       
C    
C    PREFERRED RETURN:       
 15   if ((RMIN-SR)/SR.LT.TOL) RETURN
      XT=XO  
      YT=YO  
      RMIN=SR
      FX2=FX*FX      
      FY2=FY*FY      
      FXY=2.0d0*FXY*FX*FY    
      D=FXX*FX2+FXY+FYY*FY2  
      D1=FXX*FY2-FXY+FYY*FX2    
      if (DABS(D1).GT.DABS(D)) D=D1     
      if (D.LT.1.0D-50) RETURN  
      DIV=(FX2+FY2)/D
      XO=XO-FX*DIV   
      YO=YO-FY*DIV   

   20	continue
c	write(6,*) '# ofiterations = ',icount
	goto 5

   21 cONTINUE       
C    
C    ALGORITHM HAS EXECUTED MAXIMUM NUMBER OF ITERATIONS AND IS      
C    NOT WITHIN TOLERANCE:   
      RETURN 
      END    
