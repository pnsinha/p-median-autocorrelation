	IMPLICIT double precision (A-H,O-Z)                                  0580 
	integer infreq(2),indexa(750)                                        0590 
	real name(750),PD,nametmp
      double precision X(750),Y(750),W(750),XO(2),YO(2),sR(2),
     c       POP,AREA,wmin,temp,wtemp(750)


	OPEN(5,FILE='c:\weber\weber_example.txt')
c	OPEN(7,FILE='h:\weber\check-twain.txt') 
	OPEN(8,FILE='c:\weber-out.txt') 
c	OPEN(9,FILE='h:\rsai\rsai-2000\twain-convergence.txt')
c	open(10,file='h:\twain-obj_profile.txt')

	n=400
	write(6,*) n
	pause

	wmin=1.0d10
	do 40 i=1,n                                                          0750 
c	rEAD(5,*) name(i),skip1,pop,area,pd,skip2,skip3,x(i),y(i),skip4
c	x(i)=x(i)/10000.0d0
c	y(i)=y(i)/10000.0d0
c	w(i)=pop/area

	read(5,*) name(i),x(i),y(i),w(i)
c	name(i)=i

c	if (w(i).gt.0.0d0.and.w(i).lt.wmin) wmin=w(i)
c	if (W(I).EQ.0.0d0) W(I)=1.0d-10 
	write(6,*) i,name(i)
c     WRITE(7,201) I,X(I),Y(I),W(I),name(i)                                1070 
   40 CONTINUE  
c	do 41 i=1,n
c	w(i)=w(i)/wmin
c  41 wtemp(i)=w(i)

c     WRITE(7,200) N                                                       1010 
  200 FORMAT (' NUMBER OF DESTINATIONS:',I4/                               1020 
     * ' NUMBER OF SOURCES:',I9/,                                          1030 
     * ' INPUT DATA ORDERED ASCENDING ON X:'/' DESTINATION',6X,            1040 
     * 'CO-ORDINATES',12X,'WEIGHT   ALLOCATION'/1X,61('=')//)              1050 
      DO 35 I=1,N                                                          1060 
      WRITE(6,201) I,X(I),Y(I),W(I),name(i)                                1070 
c     WRITE(7,201) I,X(I),Y(I),W(I),name(i)                                1070 
   35 CONTINUE                                                             1090 
  201 FORMAT(1X,I7,4X,3(f12.5,2X),f10.2)                                   1080 

c                                                                          0800 
C     ORDER THE DESTINATIONS BY THEIR X CO-ORDINATES:                      0810 
      DO 10 I=1,n-1                                                       0830 
      DO 9 J=I+1,N                                                        0850 
      IF (x(I).LE.x(J)) GO TO 9                                           0860 
      Temp=X(i)                                                             0870 
      X(i)=X(j)                                                           0880 
      X(j)=Temp                                                             0890 
      Temp=Y(i)                                                             0900 
      Y(i) = Y(j)                                                         0910 
      Y(j)=Temp                                                             0920 
      Temp=W(i)                                                             0930 
      W(i)=W(j)                                                           0940 
      W(j)=Temp                                                             0950 
	nametmp=name(i)
	name(i)=name(j)
	name(j)=nametmp                                                     0980 
    9 CONTINUE                                                             0990 
   10 CONTINUE                                                             1000 

c	do 998 ii=1,n

c	do 300 i=1,n
c  300	w(i)=wtemp(i)
c 	w(ii)= 1.0d-10

  498	call twain(n,x,y,w,indexa,xo,yo,sr,rmin)
	do 499 j=1,2
  499 infreq(j) = 0	
  	do 500 i=1,n
	if(indexa(i).eq.1) infreq(1)=infreq(1)+1
  500	if(indexa(i).eq.2) infreq(2)=infreq(2)+1

   	if(xo(1).lt.xo(2)) goto 45
	temp=xo(2)
	xo(2)=xo(1)
	xo(1)=temp
	temp=yo(2)
	yo(2)=yo(1)
	yo(1)=temp
   45 continue

  998	write(8,5000) xo(1),yo(1),xo(2),yo(2),rmin,(infreq(j),j=1,2)

  999 STOP   
 5000 format(1x,2(2f10.6,2x),d10.5,2x,10i5)
      END   
c
c
c
	subroutine twain(n,x,y,w,indexa,xo,yo,sr,rmin)
	IMPLICIT double precision (A-H,O-Z) 
	integer indexa(n),infreq(2)                                       0580 
      double precision X(n),Y(n),W(n),Xt(750),Yt(750),wt(750),          0600 
     c          xo(2),YO(2),sr(2),xot(2),yot(2),srt(2)

	rmin=1.0d50
	ncount=0
c
c pick the point-pairs to generate a partition
c
	do 8 i=1,n-1
	xi=x(i)
	yi=y(i)
	do 8 j=i+1,n
	dx=x(j)-xi
	dy=y(j)-yi
c
c allocate points to groups
c
	do 3 k=1,n
	indexa(k)=1
	c=(y(k)-yi)*dx - (x(k)-xi)*dy
	if(c) 2,1,3
    1	if((y(k)-yi)**2 + (x(k)-xi)**2.gt.(y(k)-y(j))**2 + (x(k)-x(j))**2)
     c  indexa(k)=2
	goto 3
    2 indexa(k)=2
    3 continue
      indexa(i)=1
	indexa(j)=2
c
c for each group ...
c
	do 5 jj=1,2
	nn=0
c
c ... transfer co-ordinates and weights to temporary vectors
c
	do 4 k=1,n
	if(indexa(k).ne.jj) goto 4
	nn=nn+1
	xt(nn)=x(k)
	yt(nn)=y(k)
	wt(nn)=w(k)
    4 continue
c
c ... and find the minimum aggregate distance
c
    5 call weber(nn,xt,yt,wt,xot(jj),yot(jj),srt(jj))
      dist=srt(1)+srt(2)
	ncount=ncount+1

c     write(10,6000) ncount,i,j,dist,xot(1),yot(1),xot(2),yot(2)

	do 499 k=1,2
  499 infreq(k) = 0	
  	do 500 k=1,n
	if(indexa(k).eq.1) infreq(1)=infreq(1)+1
  500	if(indexa(k).eq.2) infreq(2)=infreq(2)+1

c	goto 45
   	if(xot(1).lt.xot(2)) goto 45
	xtemp=xot(2)
	xot(2)=xot(1)
	xot(1)=xtemp
	ytemp=yot(2)
	yot(2)=yot(1)
	yot(1)=ytemp

	do 44 k=1,n
	if(indexa(k).eq.1) indexa(k)=0
	if(indexa(k).eq.2) indexa(k)=1
   44	if(indexa(k).eq.0) indexa(k)=2
   45 continue

c
c save the best solution
c
	if(dist.ge.rmin) goto 8
	rmin=dist
	do 6 jj=1,2
	xo(jj)=xot(jj)
	yo(jj)=yot(jj)
    6	sr(jj)=srt(jj)

c     write(9,5000) xo(1),yo(1),xo(2),yo(2),dist,(infreq(j),j=1,2)

    8 continue
c
c reconstructe the best allocation
c
	do  9 k=1,n
	indexa(k)=1
    9 if((y(k)-yo(1))**2 + (x(k)-xo(1))**2.gt.
     c   (y(k)-yo(2))**2 + (x(k)-xo(2))**2) indexa(k)=2

	return
 5000 format(1x,2(2f10.6,2x),d10.5,2x,2i5)
 6000 format(1x,i6,2i5,5f12.6)
	end
C     
C     
C     
      SUBROUTINE WEBER (N,X,Y,W,XO,YO,SR)                                  6070 
      IMPLICIT double precision (A-H,O-Z)                                  6080 
      double precision X(N),Y(N),W(N)                                             6090 
      DATA ITER,TOL /1000,1.0D-6/                                           6100 

c	write(6,*) n

C                                                                          6110 
C    THE FIRST ESTIMATE OF THE POINT OF MINIMUN AGGREGATE TRAVEL           6120 
C    IS THE CENTROID:                                                      6130 
      XO=0.0D0                                                             6140 
      YO=0.0D0                                                             6150 
      SW=0.0D0                                                             6160 
      DO 1 I=1,N                                                           6170 
      XO=XO+X(I)*W(I)                                                      6180 
      YO=YO+Y(I)*W(I)                                                      6190 
    1 SW=SW+W(I)                                                           6200 
C     IF(SW.EQ.0.0) RETURN                                                      
      XO=XO/SW                                                             6210 
      YO=YO/SW                                                             6220 
      XT=XO                                                                6230 
      YT=YO 
	
C                                                                          6250 
      SR=0.0D0                                                             6260 
      IF (N.EQ.1) RETURN                                                   6270 
C                                                                          6280 
C     BEGIN ITERATIVE LOOP.                                                6290 
C                                                                          6300 
	icount=0
      RMIN=1.0D50                                                          6310 
c      DO 20 J=1,ITER 
    4	icount=icount+1                                                      6320
c	if(icount.gt.iter) write(6,*) 'problem'
c	if(icount.gt.iter) pause
	if(icount.gt.iter) goto 21 
      FX=0.0D0                                                             6330 
      FY=0.0D0                                                             6340 
      FXX=0.0D0                                                            6350 
      FYY=0.0D0                                                            6360 
      FXY=0.0D0                                                            6370 
      SR=0.0D0                                                             6380 
C                                                                          6390 
      DO 10 I=1,N                                                          6400 
      DX=X(I)-XO                                                           6410 
      DY=Y(I)-YO                                                           6420 
      DX2=DX*DX                                                            6430 
      DY2=DY*DY                                                            6440 
      R2=DX2+DY2                                                           6450 
      IF (R2.EQ.0.0D0) GO TO 10                                            6460 
      R=DSQRT(R2)                                                          6470 
      DIV=W(I)/R                                                           6480 
      FX=FX-DX*DIV                                                         6490 
      FY=FY-DY*DIV                                                         6500 
      DIV=DIV/R2                                                           6510 
      FXX=FXX+DY2*DIV                                                      6520 
      FYY=FYY+DX2*DIV                                                      6530 
      FXY=FXY-DX*DY*DIV                                                    6540 
      SR=SR+R*W(I)                                                         6550 
 10   CONTINUE                                                             6560 
C                                                                          6570 
C     END ITERATIVE LOOP                                                   6580 
C                                                                          6590 
C                                                                          6600 
      IF (RMIN.GE.SR) GO TO 15                                             6610 
C                                                                          6620 
C    AGGREGATE DISTANCE HAS INCREASED. USE BINARY CHOP:                    6630 
      TX=XO                                                                6640 
      TY=YO                                                                6650 
      XO=0.5D0*(XO+XT)                                                     6660 
      YO=0.5D0*(YO+YT)                                                     6670 
      XT=TX                                                                6680 
      YT=TY                                                                6690 
      GO TO 20                                                             6700 
C                                                                          6710 
C    PREFERRED RETURN:                                                     6720 

   15 IF ((RMIN-SR)/SR.LT.TOL) goto 21                                      6730 
      IF (RMIN.le.sr) goto 21                                              6730 
      XT=XO                                                                6740 
      YT=YO                                                                6750 
      RMIN=SR                                                              6760 
      FX2=FX*FX                                                            6770 
      FY2=FY*FY                                                            6780 
      FXY=2.0D0*FXY*FX*FY                                                  6790 
      Dist=FXX*FX2+FXY+FYY*FY2                                             6800 
      D1=FXX*FY2-FXY+FYY*FX2                                               6810 
      IF (DABS(D1).GT.DABS(Dist)) Dist=D1                                  6820 
      IF (Dist.LT.1.0d-50) goto 21                                          6830 
      DIV=(FX2+FY2)/Dist                                                   6840 
      XO=XO-FX*DIV                                                         6850 
      YO=YO-FY*DIV                                                         6860 
   20 goto 4                                                               6870 
C                                                                          6880 
C    ALGORITHM HAS EXECUTED MAXIMUM NUMBER OF ITERATIONS AND IS            6890 
C    NOT WITHIN TOLERANCE:                                                 6900
   21 continue
	if(icount.gt.iter)
     c       write(6,*) '# iterations = ',icount
      RETURN                                                               6910 
      END                                                                  6920 
