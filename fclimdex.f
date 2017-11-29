!       last modified 2008-05-06
! add TMAXmean and TMINmean output in qc function
! in TN10p subroutine, add an 1e-5 term on all
! thresholds, to eliminate computational error. ( 3.5 may store like
! 3.5000001 or 3.49999999 in thresholds )
! changed TN10p subroutine, set missing value for monthly output the
! same level as R, eg. >10 days missing in a month, then set this month
! missing.
!  last modified 2008-06-15
! changed percentile funcion to calculate multi-level percentiles in a single
! routine, also changed threshold.

      MODULE COMM
      use, intrinsic:: iso_fortran_env, only: stderr=>error_unit,
     &     stdin=>input_unit
      IMPLICIT NONE
      SAVE

      character(20) :: STNID
      integer(4)    :: STDSPAN, BASESYEAR, BASEEYEAR, PRCPNN,
     &         SYEAR, EYEAR, TOT, YRS, BYRS, WINSIZE, SS
!     parameter(MAXYEAR=500)
      real :: LATITUDE, PRCP(500*365), TMAX(500*365),
     &     TMIN(500*365), MISSING
      integer(4) :: YMD(500*365,3), MNASTAT(500,12,3),YNASTAT(500,3),
     &        MON(12),MONLEAP(12)

      data MON/31,28,31,30,31,30,31,31,30,31,30,31/
      data MONLEAP/31,29,31,30,31,30,31,31,30,31,30,31/
      data WINSIZE/5/

      END MODULE COMM


      Program FClimDex
      use COMM
      implicit none

      character(80) :: ifile
      character(80) :: header
      integer(4) :: stnnum

      MISSING=-99.9
      SS=int(WINSIZE/2)

      stnnum=1
      open (1, file="para.txt")
      open (2, file="infilename.txt")
      read (1, '(a80)') header
77    read (1, '(a20,f10.2,3i6,i10)',end=100) STNID, LATITUDE,
     &          STDSPAN, BASESYEAR, BASEEYEAR, PRCPNN
c     print*,'##3##',STDSPAN,BASESYEAR,BASEEYEAR,PRCPNN
      read (2, '(a80)', end=100) ifile
      if(trim(ifile) == " ") then
        write(stderr,*) "Read in data filename ERROR happen in:"//
     &             "infilename.txt, line:", stnnum
        error stop
      endif
      open (6, file=trim(ifile)//"_log")
      BYRS=BASEEYEAR-BASESYEAR+1

      call qc(ifile)
      call FD(ifile)    ! FD, SU, ID, TR
      call GSL(ifile)   ! GSL
      call TXX(ifile)   ! TXx, TXn, TNx, TNn, DTR
      call Rnnmm(ifile) ! R10mm, R20mm, Rnnmm, SDII
      call RX5day(ifile)! Rx1day, Rx5day
      call CDD(ifile)   ! CDD, CWD
      call R95p(ifile)  ! R95p, R99p, PRCPTOT
      call TX10p(ifile) ! TX10p, TN10p, TX90p, TN90p
      
      stnnum=stnnum+1
      goto 77

100   close(2)
      close(1)
      stnnum=stnnum-1
      print *, "Total ",stnnum,"stations be calculated"
      end program


      pure integer function leapyear(iyear)
      implicit none
      integer,intent(in) :: iyear

      if(mod(iyear,400) == 0) then
        leapyear=1
      else
        if(mod(iyear,100) == 0) then
          leapyear=0
        else
          if(mod(iyear,4) == 0) then
            leapyear=1
          else
            leapyear=0
          endif
        endif
      endif

      end function leapyear


      subroutine percentile(x, length, nl, per, oout)
      use COMM, only: missing, stderr
      implicit none
      integer, intent(in) :: length,nl
      real, intent(in) :: x(length), per(nl)
      real, intent(out) :: oout(nl)
      real xtos(length),bb,cc
      integer nn,i
      logical ismiss,nomiss
      
      do i=1,nl
        if(per(i) > 1.or.per(i) < 0) then
          write(stderr,*) nl,i,per(i)
          error stop "Function percentile return error: parameter perc"
        endif
      enddo

      nn=0
      do i=1, length
        if(nomiss(x(i)))then
          nn=nn+1
          xtos(nn)=x(i)
        endif
      enddo

      if(nn == 0) then
        oout=MISSING
      else
        call sort(nn,xtos)
        do i=1,nl
          bb=nn*per(i)+per(i)/3.+1/3.
          cc=real(int(bb))
          if(int(cc) >= nn) then
            oout(i)=xtos(nn)
          else
            oout(i)=xtos(int(cc))+(bb-cc)*
     &          (xtos(int(cc)+1)-xtos(int(cc)))
          endif
        enddo
      endif

      end subroutine percentile

c---Sorts an array arr(1:n) into ascending numerical order using the Quicksort
c   algorithm. n is inpu; arr is replace on output by its sorted rearrangement.
c   Parameters: M is the size of subarrays sorted by straight insertion
c   and NSTACK is the required auxiliary.
      SUBROUTINE sort(n,arr)
      use comm, only: stdin, stderr
      implicit none
      INTEGER n,M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l < M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,1,-1
            if(arr(i) <= a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=0
2         arr(i+1)=a
12      continue
        if(jstack == 0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l+1) > arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l) > arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1) > arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
3       continue
          i=i+1
        if(arr(i) < a)goto 3
4       continue
          j=j-1
        if(arr(j) > a)goto 4
        if(j < i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if (jstack > NSTACK) then
            write(stderr,*) 'NSTACK too small in sort'
            read(stdin,*) 
        endif
        if(ir-i+1 >= j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END subroutine sort
C  (C) Copr. 1986-92 Numerical Recipes Software &#5,.

      subroutine qc(ifile)
      use COMM
      implicit none
      character(80), intent(in) :: ifile
      character*80 omissf, title(3), otmpfile
      integer ios, rno, tmpymd(365*500,3), i,j
      real tmpdata(365*500,3),stddata(365,500,3),stdval(365,3),
     &   m1(365,3),stdtmp
      integer kth,month,k,trno,ymiss(3),mmiss(3),stdcnt(3),
     &        missout(500,13,3),tmpcnt,leapyear
      logical ismiss,nomiss

      data title/"PRCPMISS","TMAXMISS","TMINMISS"/
      omissf=trim(ifile)//"_NASTAT"
C     print*, BASESYEAR, BASEEYEAR
      open(10, file=ifile, STATUS="OLD", IOSTAT=ios)
c     print *, ifile, ios
      if(ios /= 0) then
        write(stderr,*)  "ERROR during opening file: ", trim(ifile)
        error stop "Program STOP!!"
      endif

      otmpfile=trim(ifile)//"_prcpQC"
      open(81, file=otmpfile)
      write(81,*) "PRCP Quality Control Log File:"
      otmpfile=trim(ifile)//"_tempQC"
      open(82, file=otmpfile)
      write(82,*) "TMAX and TMIN Quality Control Log File:"

      rno=1
88    read(10,*,end=110) (tmpymd(rno,j),j=1,3),
     &                   (tmpdata(rno,j),j=1,3)
      rno=rno+1
      goto 88
110   rno=rno-1
      close(10)

      SYEAR=tmpymd(1,1)
      EYEAR=tmpymd(rno,1)
      YRS=EYEAR-SYEAR+1

      TOT=0
      do i=SYEAR,EYEAR
        do month=1,12
          if(leapyear(i)==1) then
                kth=MONLEAP(month)
          else
                kth=MON(month)
          endif
          do k=1,kth
            TOT=TOT+1
            YMD(tot,1)=i
            YMD(tot,2)=month
            YMD(tot,3)=k
          enddo
        enddo
      enddo

      j=1
      do i=1, TOT
111     if(YMD(i,1)*10000+YMD(i,2)*100+YMD(i,3) == 
     &     tmpymd(j,1)*10000+tmpymd(j,2)*100+tmpymd(j,3)) then
          PRCP(i)=tmpdata(j,1)
          TMAX(i)=tmpdata(j,2)
          TMIN(i)=tmpdata(j,3)
          j=j+1
        elseif(YMD(i,1)*10000+YMD(i,2)*100+YMD(i,3) < 
     &     tmpymd(j,1)*10000+tmpymd(j,2)*100+tmpymd(j,3)) then
          PRCP(i)=MISSING
          TMAX(i)=MISSING
          TMIN(i)=MISSING
        elseif(YMD(i,1)*10000+YMD(i,2)*100+YMD(i,3) > 
     &     tmpymd(rno,1)*10000+tmpymd(rno,2)*100+tmpymd(rno,3)) then
          PRCP(i)=MISSING
          TMAX(i)=MISSING
          TMIN(i)=MISSING
        else
          j=j+1
          goto 111
        endif
      enddo

      trno=0
      MNASTAT=0
      YNASTAT=0
      do i=SYEAR,EYEAR
        ymiss=0
        stdcnt=0
        do month=1,12
          mmiss=0
          if(leapyear(i)==1) then
            kth=MONLEAP(month)
          else
            kth=MON(month)
          endif
          do k=1,kth
            trno=trno+1
            if(TMAX(trno) <= MISSING) TMAX(trno)=MISSING
            if(TMIN(trno) <= MISSING) TMIN(trno)=MISSING
            if(PRCP(trno) <= MISSING) PRCP(trno)=MISSING
            if(TMAX(trno) < TMIN(trno).and.nomiss(TMAX(trno))
     &                   .and.nomiss(TMIN(trno))) then
              TMAX(trno)=MISSING
              TMIN(trno)=MISSING
              write(82, *) i*10000+month*100+k, "TMAX<TMIN!!"
            endif

            if(month /= 2.or.k /= 29) then
              stdcnt=stdcnt+1
              stddata(stdcnt, i-SYEAR+1, 1)=PRCP(trno)
              stddata(stdcnt, i-SYEAR+1, 2)=TMAX(trno)
              stddata(stdcnt, i-SYEAR+1, 3)=TMIN(trno)
            endif

            if((TMAX(trno) < -70..or.TMAX(trno) > 70.).and.
     &         nomiss(TMAX(trno))) then
              TMAX(trno)=MISSING
              write(82, *) i*10000+month*100+k, "TMAX over bound!!"
            endif
            if((TMIN(trno) < -70..or.TMIN(trno) > 70.).and.
     &         nomiss(TMIN(trno))) then
              TMIN(trno)=MISSING
              write(82, *) i*10000+month*100+k, "TMIN over bound!!"
            endif
            if(PRCP(trno) < 0.and.nomiss(PRCP(trno))) then
              PRCP(trno)=MISSING
              write(81,*) i*10000+month*100+k, "PRCP less then 0!!"
            endif
            if(ismiss(PRCP(trno))) then
              mmiss(1)=mmiss(1)+1
              ymiss(1)=ymiss(1)+1
            endif
            if(ismiss(TMAX(trno))) then
              mmiss(2)=mmiss(2)+1
              ymiss(2)=ymiss(2)+1
            endif
            if(ismiss(TMIN(trno))) then
              mmiss(3)=mmiss(3)+1
              ymiss(3)=ymiss(3)+1
            endif
          enddo
          do k=1,3
            missout(i-SYEAR+1,month,k)=mmiss(k)
            if (mmiss(k) > 3) then
              MNASTAT(i-SYEAR+1,month,k)=1
            endif
          enddo
        enddo
        do k=1,3
          missout(i-SYEAR+1,13,k)=ymiss(k)
          if(ymiss(k) > 15) then
            YNASTAT(i-SYEAR+1,k)=1
          endif
        enddo
      enddo

c     do i=1,YRS
c     print *,(YNASTAT(i,k),k=1,3)
c     enddo

C Calculate STD for PRCP, TMAX and TMIN; then figure out outliers
c     stdval=0.
      do i=1,365
        do j=1,3
          stdval(i,j)=0.
        enddo
      enddo
      m1=0.
      do i=1,365
        stdcnt=0
        do j=1,YRS
          do k=2,3
            if(nomiss(stddata(i,j,k))) then
              stdcnt(k)=stdcnt(k)+1
              m1(i,k)=m1(i,k)+stddata(i,j,k)
            endif
          enddo
        enddo
        do k=2,3
          if(stdcnt(k) > 0) then
            m1(i,k)=m1(i,k)/real(stdcnt(k))
          endif
        enddo
        stdtmp=0.
        do j=1,YRS
          do k=2,3
            if(stdcnt(k) > 2.and.nomiss(stddata(i,j,k))) then
               stdval(i,k)=stdval(i,k)+
     &         (stddata(i,j,k)-m1(i,k))**2./(real(stdcnt(k))-1.)
            endif
          enddo
        enddo

        do k=2,3
          if(stdcnt(k) > 2) then
            stdval(i,k)=stdval(i,k)**0.5
          else 
            stdval(i,k)=MISSING
          endif
        enddo
      enddo

      trno=0
      do i=SYEAR,EYEAR
        tmpcnt=0
        do month=1,12
          if(leapyear(i)==1) then
            kth=MONLEAP(month)
          else 
            kth=MON(month)
          endif
          do k=1, kth
            trno=trno+1
            if(month /= 2.or.k /= 29) tmpcnt=tmpcnt+1
            if(nomiss(stdval(tmpcnt,2)))then
              if(abs(TMAX(trno)-m1(tmpcnt,2)) > 
     &          stdval(tmpcnt,2)*STDSPAN.and.nomiss(TMAX(trno)))
     &      write(82, *) "Outlier: ", i, month, k, "TMAX: ", TMAX(trno),
     &            "Lower limit:",m1(tmpcnt,2)-stdval(tmpcnt,2)*STDSPAN,
     &            "Upper limit:",m1(tmpcnt,2)+stdval(tmpcnt,2)*STDSPAN
            endif
            if(nomiss(stdval(tmpcnt,3)))then
              if(abs(TMIN(trno)-m1(tmpcnt,3)) > 
     &          stdval(tmpcnt,3)*STDSPAN.and.nomiss(TMIN(trno)))
     &      write(82, *) "Outlier: ", i, month, k, "TMIN: ", TMIN(trno),
     &            "Lower limit:",m1(tmpcnt,3)-stdval(tmpcnt,3)*STDSPAN,
     &            "Upper limit:",m1(tmpcnt,3)+stdval(tmpcnt,3)*STDSPAN
            endif
          enddo ! end do day
        enddo !end do month
      enddo ! end do year

      open(20,file=trim(omissf))
      do i=SYEAR,EYEAR
        do k=1,3
          write(20,'(i4,2x,a8,13i4)')
     &    i,title(k),(missout(i+1-SYEAR,j,k),j=1,13)
        enddo
      enddo
      close(20)
      close(81)
      close(82)

C  QC part finished, prepared data set: YMD(3), PRCP, TMAX & TMIN
C  and NASTAT dataset for missing values monthly and annual
      end subroutine qc


      subroutine FD(ifile)
      use COMM
      implicit none
      character(80), intent(in) :: ifile

      integer year, trno, kth, month,day,i,j,leapyear
      real oout(YRS,4)
      character*2 chrtmp(4)
      character*80 ofile
C oout(,1)--FD, oout(,2)--SU, oout(,3)--ID, oout(,4)--TR      
      data chrtmp/"FD","SU","ID","TR"/
      logical nomiss

      trno=0
      oout=0
      do i=1,YRS
        year=i+SYEAR-1
        do month=1,12
          if(leapyear(year)==1) then
                  kth=MONLEAP(month)
          else
                  kth=MON(month)
          endif
          do day=1,kth
            trno=trno+1
            if(YMD(trno,3) /= day) then
              error stop 'ERROR1 at FD!!!'
            endif
            if(nomiss(TMIN(trno)).and.TMIN(trno) < 0) 
     &          oout(i,1)=oout(i,1)+1
            if(nomiss(TMAX(trno)).and.TMAX(trno) > 25) 
     &          oout(i,2)=oout(i,2)+1
            if(nomiss(TMAX(trno)).and.TMAX(trno) < 0) 
     &          oout(i,3)=oout(i,3)+1
            if(nomiss(TMIN(trno)).and.TMIN(trno) > 20) 
     &          oout(i,4)=oout(i,4)+1
          enddo
        enddo
      enddo

      do i=1,YRS
        if(YNASTAT(i,2)==1) then
                oout(i,2)=MISSING  ! SU
                oout(i,3)=MISSING  ! ID
        endif
        if(YNASTAT(i,3)==1) then
                oout(i,1)=MISSING  ! FD
                oout(i,4)=MISSING  ! TR
        endif
      enddo

      do j=1,4
        ofile=trim(ifile)//"_"//chrtmp(j)
        open(22,file=ofile)
          write(22, *) "year    ", chrtmp(j)
          do i=1,YRS
            write(22, '(i8,f8.1)') i+SYEAR-1, oout(i,j)
          enddo
        close(22)
      enddo

      end subroutine fd


      subroutine GSL(ifile)
      use COMM
      implicit none
      character(80), intent(in) :: ifile

      character(80) ofile
      integer year,cnt,kth,month,day,marks,marke,leapyear,i

      real TG,oout,strt(YRS),ee(YRS)
      logical ismiss,nomiss


      strt=MISSING
      ee=MISSING
      cnt=0
      do i=1,YRS
        year=i+SYEAR-1
        marks=0
        marke=0
        do month=1,6
          if(leapyear(year) == 1) then
            kth=MONLEAP(month)
          else
            kth=MON(month)
          endif
          do day=1,kth
            cnt=cnt+1
            if(YMD(cnt,1)*10000+YMD(cnt,2)*100+YMD(cnt,3) /= 
     &         year*10000+month*100+day) then
              write(stderr,*) 'date count ERROR in GSL!'
              write(stderr,*) YMD(cnt,1)*10000 
     &                + YMD(cnt,2)*100+YMD(cnt,3),
     &                year*10000+month*100+day
              error stop
            endif
            if(nomiss(TMAX(cnt)).and.nomiss(TMIN(cnt))) then
              TG=(TMAX(cnt)+TMIN(cnt))/2.
            else
              TG=MISSING
            endif
            if(LATITUDE < 0) then
              if(nomiss(TG).and.TG < 5.)then
                marke=marke+1
              else
                marke=0
              endif
              if(marke >= 6.and.i > 1.and.ismiss(ee(i-1)))then
                ee(i-1)=cnt-5
              endif
            else
              if(nomiss(TG).and.TG > 5.)then
                marks=marks+1
              else
                marks=0
              endif
              if(marks >= 6.and.ismiss(strt(i)))then
                strt(i)=cnt-5
              endif
            endif
          enddo
        enddo 
        if(LATITUDE < 0.and.i > 1) then
          if(ismiss(ee(i-1)).and.nomiss(strt(i-1))) then
            ee(i-1)=cnt
          endif
        endif
        marks=0
        marke=0
        do month=7,12
          do day=1,MON(month)
            cnt=cnt+1
            if(nomiss(TMAX(cnt)).and.nomiss(TMIN(cnt))) then
              TG=(TMAX(cnt)+TMIN(cnt))/2.
            else
              TG=MISSING
            endif
            if(LATITUDE < 0) then
              if(nomiss(TG).and.TG > 5.)then
                marks=marks+1
              else
                marks=0
              endif
              if(marks >= 6.and.ismiss(strt(i)))then
                strt(i)=cnt-5
              endif
            else
              if(nomiss(TG).and.TG < 5.)then
                marke=marke+1
              else
                marke=0
              endif
              if(marke >= 6.and.ismiss(ee(i)))then
                ee(i)=cnt-5
              endif
            endif
          enddo
        enddo
        if(ismiss(ee(i)).and.nomiss(strt(i))) then
          ee(i)=cnt
        endif
      enddo

      ofile=trim(ifile)//"_GSL"
      open(22,file=ofile)
      write(22,*) "  year    gsl  "
      do i=1,YRS
        year=i+SYEAR-1
        if(nomiss(strt(i)).and.nomiss(ee(i)).and.
     &     YNASTAT(i,2) /= 1.and.YNASTAT(i,3) /= 1)then
          oout=ee(i)-strt(i)
        elseif(ismiss(strt(i)).or.ismiss(ee(i))) then
          oout=0.
        endif
        if(YNASTAT(i,2) == 1.or.YNASTAT(i,3) == 1) oout=MISSING
c       if(year == 1923) print *, year, YNASTAT(i,2),YNASTAT(i,3)
        write(22,'(i6,f8.1)') year, oout
      enddo
      close(22)

      end subroutine gsl


      subroutine TXX(ifile)
      use COMM
      implicit none
      character(80), intent(in) :: ifile
      character*80 ofile
      character*3 chrtmp(4)

      integer year,month,day,kth,cnt,nn,i,j,k,leapyear
      real oout(YRS,12,4),yout(YRS,4), dtr(YRS,13)
      logical ismiss,nomiss

      data chrtmp/"TXx","TXn","TNx","TNn"/

      oout=MISSING
      dtr=0.
      cnt=0
      do i=1,YRS
        year=i+SYEAR-1
        do month=1,12
          if(leapyear(year)==1) then
            kth=MONLEAP(month)
          else
            kth=MON(month)
          endif
          nn=0
          do day=1,kth
            cnt=cnt+1
            if(nomiss(TMAX(cnt)).and.nomiss(TMIN(cnt))) then
              dtr(i,month)=dtr(i,month)+(TMAX(cnt)-TMIN(cnt))
              nn=nn+1
            endif
            if(nomiss(TMAX(cnt)).and.(ismiss(oout(i,month,1)).or.
     &         TMAX(cnt) > oout(i,month,1))) then
              oout(i,month,1)=TMAX(cnt) ! TXX
            endif
            if(nomiss(TMAX(cnt)).and.(ismiss(oout(i,month,2)).or.
     &         TMAX(cnt) < oout(i,month,2))) then
              oout(i,month,2)=TMAX(cnt) ! TXN
            endif
            if(nomiss(TMIN(cnt)).and.(ismiss(oout(i,month,3)).or.
     &         TMIN(cnt) > oout(i,month,3))) then
              oout(i,month,3)=TMIN(cnt) ! TNX
            endif
            if(nomiss(TMIN(cnt)).and.(ismiss(oout(i,month,4)).or.
     &         TMIN(cnt) < oout(i,month,4))) then
              oout(i,month,4)=TMIN(cnt) ! TNN
            endif
          enddo 
          if(nn > 0.and.MNASTAT(i,month,2) == 0.and.
     &          MNASTAT(i,month,3) == 0) then
            dtr(i,month)=dtr(i,month)/nn
          else
            dtr(i,month)=MISSING
          endif
          if(MNASTAT(i,month,2) == 1)then
            oout(i,month,1)=MISSING
            oout(i,month,2)=MISSING
          endif
          if(MNASTAT(i,month,3) == 1)then
            oout(i,month,3)=MISSING
            oout(i,month,4)=MISSING
          endif
        enddo
      enddo

      yout=MISSING
      do i=1,YRS
        nn=0
        do month=1,12
          if(nomiss(oout(i,month,1)).and.(ismiss(yout(i,1)).or.
     &       oout(i,month,1) > yout(i,1))) then
            yout(i,1)=oout(i,month,1)
          endif
          if(nomiss(oout(i,month,2)).and.(ismiss(yout(i,2)).or.
     &       oout(i,month,2) < yout(i,2))) then
            yout(i,2)=oout(i,month,2)
          endif
          if(nomiss(oout(i,month,3)).and.(ismiss(yout(i,3)).or.
     &       oout(i,month,3) > yout(i,3))) then
            yout(i,3)=oout(i,month,3)
          endif
          if(nomiss(oout(i,month,4)).and.(ismiss(yout(i,4)).or.
     &       oout(i,month,4) < yout(i,4))) then
            yout(i,4)=oout(i,month,4)
          endif
          if(nomiss(dtr(i,month))) then
            dtr(i,13)=dtr(i,13)+dtr(i,month)
            nn=nn+1
          endif
        enddo
        if(nn > 0.and.YNASTAT(i,2) == 0.and.YNASTAT(i,3) == 0) then
          dtr(i,13)=dtr(i,13)/nn
        else
          dtr(i,13)=MISSING
        endif
        if(YNASTAT(i,2) == 1) then
          yout(i,1)=MISSING
          yout(i,2)=MISSING
        endif
        if(YNASTAT(i,3) == 1) then
          yout(i,3)=MISSING
          yout(i,4)=MISSING
        endif
      enddo
      do k=1,4
        ofile=trim(ifile)//"_"//chrtmp(k)
        open(22,file=ofile)
        write(22, *) "  year  jan   feb   mar   apr   may   jun  ",
     &               " jul   aug   sep   oct   nov   dec annual"
        do i=1,YRS
          write(22,'(i6,13f6.1)') i+SYEAR-1,(oout(i,j,k),j=1,12),
     &                            yout(i,k)
        enddo
        close(22)
      enddo

      ofile=trim(ifile)//"_DTR"
      open(22,file=ofile)
      write(22, *) "  year  jan   feb   mar   apr   may   jun  ",
     &             " jul   aug   sep   oct   nov   dec annual"
      do i=1,YRS
        write(22,'(i6,13f7.2)') i+SYEAR-1,(dtr(i,j),j=1,13)
      enddo
      close(22)

      end subroutine txx


      subroutine Rnnmm(ifile)
      use COMM
      implicit none
      character(80), intent(in) :: ifile
      character(80) :: ofile
      character(5) chrtmp(3)
      integer year,month,day,kth,cnt,nn,leapyear,i,k

      real oout(YRS,3),sdii(YRS)


      data chrtmp/"R10mm","R20mm","Rnnmm"/
      cnt=0
      oout=0.
      sdii=0.
      do i=1,YRS
        nn=0
        year=i+SYEAR-1
        do month=1,12
          if(leapyear(year) == 1) then
            kth=MONLEAP(month)
          else
            kth=MON(month)
          endif
          do day=1,kth
            cnt=cnt+1
            if(PRCP(cnt) >= 1.) then
              sdii(i)=sdii(i)+PRCP(cnt)
              nn=nn+1
            endif
            if(PRCP(cnt) >= 10.) oout(i,1)=oout(i,1)+1.
            if(PRCP(cnt) >= 20.) oout(i,2)=oout(i,2)+1.
            if(PRCP(cnt) >= PRCPNN) oout(i,3)=oout(i,3)+1.
          enddo
        enddo
        if(nn > 0) then
          sdii(i)=sdii(i)/nn
        endif
      enddo

      do i=1,YRS
        if(YNASTAT(i,1) == 1) then
          do k=1,3
            oout(i,k)=MISSING
          enddo
          sdii(i)=MISSING
        endif
      enddo

      do k=1,3
        ofile=trim(ifile)//"_"//chrtmp(k)
        open(22,file=ofile)
        write(22,*) "  year  ",chrtmp(k)
        do i=1,YRS
          write(22,'(i6,f8.1)') i+SYEAR-1, oout(i,k)
        enddo
        close(22)
      enddo

      ofile=trim(ifile)//"_SDII"
      open(22,file=ofile)
      write(22,*) "  year    sdii"
      do i=1,YRS
        write(22,'(i6,f8.1)') i+SYEAR-1, sdii(i)
      enddo
      close(22)

      end subroutine rnnmm


      subroutine RX5day(ifile)
      use COMM
      implicit none
      character(80), intent(in) :: ifile
      character(80) :: ofile

      integer year, month, day,cnt,i,j,k,kth,leapyear

      real r1(YRS,13), r5(YRS,13), r5prcp
      logical ismiss,nomiss

      cnt=0
      r1=MISSING
      r5=MISSING
      do i=1,YRS
        year=i+SYEAR-1
        do month=1,12
          if(leapyear(year) == 1) then
            kth=MONLEAP(month)
          else
            kth=MON(month)
          endif
          do day=1,kth
            cnt=cnt+1
c           if(year == 1904.and.month == 1) then
c             print *, year, month, day, PRCP(cnt)
c           endif
            if(cnt > 5)then
              r5prcp=0.
              do k=cnt-4,cnt
                if(nomiss(PRCP(k)))then
                  r5prcp=r5prcp+PRCP(k)
                endif
              enddo
            else
              r5prcp=MISSING
            endif
            if(nomiss(PRCP(cnt)).and.(ismiss(r1(i,month))
     &         .or.PRCP(cnt) > r1(i,month))) then
              r1(i,month)=PRCP(cnt)
            endif
            if(nomiss(PRCP(cnt)).and.r5prcp > r5(i,month)) then
              r5(i,month)=r5prcp
            endif
          enddo
          if(MNASTAT(i,month,1) == 1) then
            r1(i,month)=MISSING
            r5(i,month)=MISSING
          endif
          if(nomiss(r1(i,month)).and.(ismiss(r1(i,13))
     &       .or.r1(i,month) > r1(i,13))) then
            r1(i,13)=r1(i,month)
          endif
          if(nomiss(r5(i,month)).and.(ismiss(r5(i,13))
     &       .or.r5(i,month) > r5(i,13))) then
            r5(i,13)=r5(i,month)
          endif
        enddo
        if(YNASTAT(i,1) == 1) then
          r1(i,13)=MISSING
          r5(i,13)=MISSING
        endif
      enddo

      ofile=trim(ifile)//"_RX1day"
      open(22,file=ofile)
      write(22, *) "  year  jan    feb    mar    apr    may    jun  ",
     &             "  jul    aug    sep    oct    nov    dec  annual"
      do i=1,YRS
        write(22, '(i6,13f7.1)') i+SYEAR-1,(r1(i,j),j=1,13)
      enddo
      close(22)
      ofile=trim(ifile)//"_RX5day"
      open(22,file=ofile)
      write(22, *) "  year  jan    feb    mar    apr    may    jun  ",
     &             "  jul    aug    sep    oct    nov    dec  annual"
      do i=1,YRS
        write(22, '(i6,13f7.1)') i+SYEAR-1,(r5(i,j),j=1,13)
      enddo
      close(22)

      end subroutine rx5day


      subroutine CDD(ifile)
      use COMM
      implicit none
      character(80) :: ifile, ofile

      integer year, month, day, kth, cnt,i,leapyear

      real ocdd(YRS), ocwd(YRS), nncdd, nncwd
      logical ismiss

      cnt=0
      ocdd=0.
      ocwd=0.
      do i=1,YRS
        if(i==1)nncdd=0.
        if(i==1)nncwd=0.
        year=i+SYEAR-1
        do month=1,12
          if(leapyear(year) == 1) then
            kth=MONLEAP(month)
          else
            kth=MON(month)
          endif
          do day=1,kth
            cnt=cnt+1
            if(ismiss(PRCP(cnt))) then
              nncdd=0.
              nncwd=0.
            elseif(PRCP(cnt) < 1) then
              nncdd=nncdd+1.
              if(nncwd > ocwd(i)) ocwd(i)=nncwd
              nncwd=0.
            else
              nncwd=nncwd+1.
              if(nncdd > ocdd(i)) ocdd(i)=nncdd
              nncdd=0.
            endif
c           if(year == 1959.and.month == 12) then 
c                   print *, month,day,nncdd, ocdd(i)
c           endif
          enddo
        enddo

        if(ocwd(i) < nncwd) then
          if(year == EYEAR) then
                  ocwd(i)=nncwd
          elseif(PRCP(cnt+1) < 1..or.ismiss(PRCP(cnt+1)))then
                  ocwd(i)=nncwd
          endif
        endif

        if(ocdd(i) < nncdd) then
          if(year == EYEAR) then
                  ocdd(i)=nncdd
          elseif(PRCP(cnt+1) >= 1..or.ismiss(PRCP(cnt+1)))then
                  ocdd(i)=nncdd
          endif
          if(ocdd(i) == 0) ocdd(i)=MISSING
        endif

        if(YNASTAT(i,1) == 1) then
          ocdd(i)=MISSING
          ocwd(i)=MISSING
        endif
      enddo

      ofile=trim(ifile)//"_CDD"
      open(22,file=ofile)
      write(22,*) "  year   cdd  "
      do i=1,YRS
        write(22,'(i6,f8.1)') i+SYEAR-1, ocdd(i)
      enddo
      close(22)

      ofile=trim(ifile)//"_CWD"
      open(22,file=ofile)
      write(22,*) "  year   cwd  "
      do i=1,YRS
        write(22,'(i6,f8.1)') i+SYEAR-1, ocwd(i)
      enddo
      close(22)

      end subroutine cdd


      subroutine R95p(ifile)
      use COMM
      implicit none
      character*80 ifile
      character*80 ofile
      integer year, month, day, kth,cnt,leng,i,leapyear

      real r95out(YRS), prcptmp(TOT),r99out(YRS), prcpout(YRS), p95, 
     &     p99,rlev(2),rtmp(2)
      logical nomiss

      cnt=0
      leng=0
      prcptmp=MISSING
      do i=1,YRS
        year=i+SYEAR-1
        do month=1,12
          if(leapyear(year) == 1)then
            kth=MONLEAP(month)
          else
            kth=MON(month)
          endif
          do day=1,kth
            cnt=cnt+1
            if(year >= BASESYEAR.and.year <= BASEEYEAR.and.
     &         nomiss(PRCP(cnt)).and.PRCP(cnt) >= 1.)then
              leng=leng+1
              prcptmp(leng)=PRCP(cnt)
            endif
          enddo
        enddo
      enddo
      rlev(1)=0.95
      rlev(2)=0.99
      call percentile(prcptmp,leng,2,rlev,rtmp)
      p95=rtmp(1)
      p99=rtmp(2)
c     p95=percentile(prcptmp,leng,0.95)
c     p99=percentile(prcptmp,leng,0.99)

      cnt=0
      r95out=0.
      r99out=0.
      prcpout=0.
      do i=1,YRS
        year=i+SYEAR-1
        do month=1,12
          if(leapyear(year) == 1)then
            kth=MONLEAP(month)
          else
            kth=MON(month)
          endif
          do day=1,kth
            cnt=cnt+1
            if(PRCP(cnt) >= 1..and.nomiss(PRCP(cnt)))then
              prcpout(i)=prcpout(i)+PRCP(cnt)
              if(PRCP(cnt) > p95) r95out(i)=r95out(i)+PRCP(cnt)
              if(PRCP(cnt) > p99) r99out(i)=r99out(i)+PRCP(cnt)
            endif
          enddo
        enddo
        if(YNASTAT(i,1) == 1) then
          prcpout(i)=MISSING
          r95out(i)=MISSING
          r99out(i)=MISSING
        endif
      enddo

      ofile=trim(ifile)//"_PRCPTOT"
      open(22,file=ofile)
      write(22,*) "  year prcptot "
      do i=1,YRS
        write(22,'(i6,f8.1)') i+SYEAR-1, prcpout(i)
      enddo
      close(22)
      ofile=trim(ifile)//"_R95p"
      open(22,file=ofile)
      write(22,*) "  year  r95p  "
      do i=1,YRS
        write(22,'(i6,f8.1)') i+SYEAR-1, r95out(i)
      enddo
      close(22)
      ofile=trim(ifile)//"_R99p"
      open(22,file=ofile)
      write(22,*) "  year  r99p  "
      do i=1,YRS
        write(22,'(i6,f8.1)') i+SYEAR-1, r99out(i)
      enddo
      close(22)

      end subroutine r95p


      subroutine TX10p(ifile)
      use COMM
      implicit none
      character(80) :: ifile, ofile
      integer :: leapyear

      integer year, month, day, kth, cnt, nn,  missxcnt, missncnt,
     &        iter, cntx, cntn,i,j,byear,flgtn,flgtx,flg,idum

      real tminbase(TOT),txdata(BYRS,365+2*SS),
     &     tndata(BYRS,365+2*SS),thresan10(365),txdtmp(BYRS,365),
     &     tndtmp(BYRS,365),tnboot(BYRS,365+2*SS),
     &     txboot(BYRS,365+2*SS),thresan90(365),thresax10(365),
     &     thresax90(365),tx10out(YRS,13),tx90out(YRS,13),
     &     thresax50(365),thresan50(365),tx50out(YRS,13),
     &     tn50out(YRS,13),thresbx50(365,BYRS,BYRS-1),
     &     thresbn50(365,BYRS,BYRS-1),threstmp(365,3),rlevs(3),
     &     tn10out(YRS,13),tn90out(YRS,13),thresbn90(365,BYRS,BYRS-1),
     &     thresbn10(365,BYRS,BYRS-1),thresbx90(365,BYRS,BYRS-1),
     &     thresbx10(365,BYRS,BYRS-1),wsdi(YRS),csdi(YRS)
      logical nomiss

      data rlevs/0.1,0.5,0.9/

      cnt=0
      nn=0
      txdtmp=MISSING
      tndtmp=MISSING
      do i=1,YRS
        year=i+SYEAR-1
        nn=0
        do month=1,12
          if(leapyear(year) == 1) then
            kth=MONLEAP(month)
          else
            kth=MON(month)
          endif
          do day=1,kth
            cnt=cnt+1
            if(year >= BASESYEAR.and.year <= BASEEYEAR.and.(month /= 2
     &         .or.day /= 29))then
              nn=nn+1
              txdtmp(i+SYEAR-BASESYEAR,nn)=TMAX(cnt)
              tndtmp(i+SYEAR-BASESYEAR,nn)=TMIN(cnt)
            endif
          enddo
        enddo
        if(year >= BASESYEAR.and.year <= BASEEYEAR.and.nn /= 365)then
          write(stderr,*) "date count error in TX10p!", nn
          error stop
        endif
      enddo
      
      do i=1,BYRS
        do j=1,SS
          if(i == 1) then
            tndata(i,j)=tndtmp(i,1)
            txdata(i,j)=txdtmp(i,1)
          else 
            tndata(i,j)=tndtmp(i-1,365+j-SS)
            txdata(i,j)=txdtmp(i-1,365+j-SS)
          endif
        enddo
        do j=1,365
          tndata(i,j+SS)=tndtmp(i,j)
          txdata(i,j+SS)=txdtmp(i,j)
        enddo
        do j=1,SS
          if(i == BYRS)then
            tndata(i,j+365+SS)=tndtmp(i,365)
            txdata(i,j+365+SS)=txdtmp(i,365)
          else
            tndata(i,j+365+SS)=tndtmp(i+1,j)
            txdata(i,j+365+SS)=txdtmp(i+1,j)
          endif
        enddo
      enddo

      flgtn=0
      flgtx=0
c     call threshold(tndata,.1,thresan10, flgtn)
      call threshold(tndata,rlevs,3,threstmp,flgtn)
      do i=1,365
        thresan10(i)=threstmp(i,1)-1e-5
        thresan50(i)=threstmp(i,2)+1e-5
        thresan90(i)=threstmp(i,3)+1e-5
      enddo
c     thresan10=thresan10-1e-5
      if(flgtn == 1) then
        write(stderr,*) "TMIN Missing value overflow in exceedance rate"
        tn10out=MISSING
        tn50out=MISSING
        tn90out=MISSING
c     else
c       call threshold(tndata,.5,thresan50, flgtn)
c       thresan50=thresan50+1e-5
c       call threshold(tndata,.9,thresan90, flgtn)
c       thresan90=thresan90+1e-5
      endif

      call threshold(txdata,.1,thresax10, flgtx)
      call threshold(txdata,rlevs,3,threstmp,flgtx)
      do i=1,365
        thresax10(i)=threstmp(i,1)-1e-5
        thresax50(i)=threstmp(i,2)+1e-5
        thresax90(i)=threstmp(i,3)+1e-5
      enddo
c     thresax10=thresax10-1e-5
      if(flgtx == 1) then
        write(stderr,*) "TMAX Missing value overflow in exceedance rate"
        tx10out=MISSING
        tx50out=MISSING
        tx90out=MISSING
c     else
c       call threshold(txdata,.5,thresax50, flgtx)
c       thresax50=thresax50+1e-5
c       call threshold(txdata,.9,thresax90, flgtx)
c       thresax90=thresax90+1e-5
      endif

      do i=1,BYRS
        txboot=txdata
        tnboot=tndata
        nn=0
        do iter=1,BYRS
          if(iter /= i) then
            nn=nn+1
            do day=1,365+2*SS
              if(flgtx == 0) txboot(i,day)=txboot(iter,day)
              if(flgtn == 0) tnboot(i,day)=tnboot(iter,day)
            enddo
            if(flgtx == 0)then
              call threshold(txboot,rlevs,3,threstmp,flg)
              do day=1,365
                thresbx90(day,i,nn)=threstmp(day,3)+1e-5
                thresbx50(day,i,nn)=threstmp(day,2)+1e-5
                thresbx10(day,i,nn)=threstmp(day,1)-1e-5
              enddo
            endif

            if(flgtn == 0) then
              call threshold(tnboot,rlevs,3,threstmp,flg)
              do day=1,365
                thresbn90(day,i,nn)=threstmp(day,3)+1e-5
                thresbn50(day,i,nn)=threstmp(day,2)+1e-5
                thresbn10(day,i,nn)=threstmp(day,1)-1e-5
              enddo
            endif
          endif
        enddo
      enddo

      if(flgtx == 0)then
        tx10out=0.
        tx50out=0.
        tx90out=0.
      endif
      if(flgtn == 0)then
        tn10out=0.
        tn50out=0.
        tn90out=0.
      endif
      cnt=0
      do i=1,YRS
        year=i+SYEAR-1
        byear=year-BASESYEAR+1
        nn=0
        do month=1,12
          missncnt=0
          missxcnt=0
          if(leapyear(year) == 1)then
            kth=MONLEAP(month)
          else
            kth=MON(month)
          endif
          do day=1,kth
            if(month /= 2.or.day /= 29) nn=nn+1
            cnt=cnt+1
            if(nomiss(TMAX(cnt)))then
              if(year < BASESYEAR.or.year > BASEEYEAR) then
                if(TMAX(cnt) > thresax90(nn)) tx90out(i,month)=
     &                   tx90out(i,month)+1
                if(TMAX(cnt) > thresax50(nn)) tx50out(i,month)=
     &                   tx50out(i,month)+1
                if(TMAX(cnt) < thresax10(nn)) tx10out(i,month)=
     &                   tx10out(i,month)+1
              else
                do iter=1,BYRS-1
c                   if(byear > 30.or.iter > 29) then
c                     print *, i, year,month,day
c                     stop
c                   endif
                  if(TMAX(cnt) > thresbx90(nn,byear,iter))then
                    tx90out(i,month)=tx90out(i,month)+1
                  endif
                  if(TMAX(cnt) > thresbx50(nn,byear,iter))then
                    tx50out(i,month)=tx50out(i,month)+1
                  endif
                  if(TMAX(cnt) < thresbx10(nn,byear,iter))
     &              tx10out(i,month)=tx10out(i,month)+1
                enddo
              endif
            else
              missxcnt=missxcnt+1
            endif
            if(nomiss(TMIN(cnt)))then
              if(year < BASESYEAR.or.year > BASEEYEAR) then
                if(TMIN(cnt) > thresan90(nn)) tn90out(i,month)=
     &                 tn90out(i,month)+1
                if(TMIN(cnt) > thresan50(nn)) tn50out(i,month)=
     &                 tn50out(i,month)+1
                if(TMIN(cnt) < thresan10(nn)) tn10out(i,month)=
     &                   tn10out(i,month)+1
              else
                do iter=1,BYRS-1
                  if(TMIN(cnt) > thresbn90(nn,byear,iter))
     &              tn90out(i,month)=tn90out(i,month)+1
                  if(TMIN(cnt) > thresbn50(nn,byear,iter))
     &              tn50out(i,month)=tn50out(i,month)+1
                  if(TMIN(cnt) < thresbn10(nn,byear,iter))
     &              tn10out(i,month)=tn10out(i,month)+1
                enddo
              endif
            else
              missncnt=missncnt+1
            endif
          enddo ! do day=1,kth

c         if(year >= BASESYEAR.and.year <= BASEEYEAR)then
c           print *, year,month,tx10out(i,month),tx90out(i,month),
c    &        tn10out(i,month),tn90out(i,month),missxcnt,missncnt
c         endif

          if(year >= BASESYEAR.and.year <= BASEEYEAR)then
            tn90out(i,month)=tn90out(i,month)/(BYRS-1.)
            tn50out(i,month)=tn50out(i,month)/(BYRS-1.)
            tn10out(i,month)=tn10out(i,month)/(BYRS-1.)
            tx90out(i,month)=tx90out(i,month)/(BYRS-1.)
            tx50out(i,month)=tx50out(i,month)/(BYRS-1.)
            tx10out(i,month)=tx10out(i,month)/(BYRS-1.)
          endif

c         if(year == 1952) then
c           print *,year,month,tn10out(i,month),tn90out(i,month)
c         endif

          if(missxcnt <= 10.and.flgtx == 0)then
            tx90out(i,13)=tx90out(i,13)+tx90out(i,month)
            tx90out(i,month)=tx90out(i,month)*100./(kth-missxcnt)
            tx50out(i,13)=tx50out(i,13)+tx50out(i,month)
            tx50out(i,month)=tx50out(i,month)*100./(kth-missxcnt)
            tx10out(i,13)=tx10out(i,13)+tx10out(i,month)
            tx10out(i,month)=tx10out(i,month)*100./(kth-missxcnt)
          else
            tx90out(i,month)=MISSING
            tx50out(i,month)=MISSING
            tx10out(i,month)=MISSING
          endif
          if(missncnt <= 10.and.flgtn == 0)then
            tn90out(i,13)=tn90out(i,13)+tn90out(i,month)
            tn90out(i,month)=tn90out(i,month)*100./(kth-missncnt)
            tn50out(i,13)=tn50out(i,13)+tn50out(i,month)
            tn50out(i,month)=tn50out(i,month)*100./(kth-missncnt)
            tn10out(i,13)=tn10out(i,13)+tn10out(i,month)
            tn10out(i,month)=tn10out(i,month)*100./(kth-missncnt)
          else
            tn90out(i,month)=MISSING
            tn50out(i,month)=MISSING
            tn10out(i,month)=MISSING
          endif
        enddo ! do month=1,12
        if(YNASTAT(i,3) == 1.or.flgtn == 1) then
          tn10out(i,13)=MISSING
          tn50out(i,13)=MISSING
          tn90out(i,13)=MISSING
        else
          tn10out(i,13)=tn10out(i,13)*100/365.
          tn50out(i,13)=tn50out(i,13)*100/365.
          tn90out(i,13)=tn90out(i,13)*100/365.
        endif
        if(YNASTAT(i,2) == 1.or.flgtx == 1) then
          tx10out(i,13)=MISSING
          tx50out(i,13)=MISSING
          tx90out(i,13)=MISSING
        else
          tx10out(i,13)=tx10out(i,13)*100/365.
          tx50out(i,13)=tx50out(i,13)*100/365.
          tx90out(i,13)=tx90out(i,13)*100/365.
        endif
      enddo

130   continue
      if(flgtx == 0)then
      ofile=trim(ifile)//"_TX90p"
      open(22,file=ofile)
      write(22, *) "  year  jan   feb   mar   apr   may   jun  ",
     &             " jul   aug   sep   oct   nov   dec annual"
      do i=1,YRS
        write(22,'(i6,13f7.2)') i+SYEAR-1,(tx90out(i,j),j=1,13)
      enddo
      close(22)

      ofile=trim(ifile)//"_TX50p"
      open(22,file=ofile)
      write(22, *) "  year  jan   feb   mar   apr   may   jun  ",
     &             " jul   aug   sep   oct   nov   dec annual"
      do i=1,YRS
        write(22,'(i6,13f7.2)') i+SYEAR-1,(tx50out(i,j),j=1,13)
      enddo
      close(22)

      ofile=trim(ifile)//"_TX10p"
      open(22,file=ofile)
      write(22, *) "  year  jan   feb   mar   apr   may   jun  ",
     &             " jul   aug   sep   oct   nov   dec annual"
      do i=1,YRS
        write(22,'(i6,13f7.2)') i+SYEAR-1,(tx10out(i,j),j=1,13)
      enddo
      close(22)
      endif

      if(flgtn == 0)then
      ofile=trim(ifile)//"_TN90p"
      open(22,file=ofile)
      write(22, *) "  year  jan   feb   mar   apr   may   jun  ",
     &             " jul   aug   sep   oct   nov   dec annual"
      do i=1,YRS
        write(22,'(i6,13f7.2)') i+SYEAR-1,(tn90out(i,j),j=1,13)
      enddo
      close(22)

      ofile=trim(ifile)//"_TN50p"
      open(22,file=ofile)
      write(22, *) "  year  jan   feb   mar   apr   may   jun  ",
     &             " jul   aug   sep   oct   nov   dec annual"
      do i=1,YRS
        write(22,'(i6,13f7.2)') i+SYEAR-1,(tn50out(i,j),j=1,13)
      enddo
      close(22)

      ofile=trim(ifile)//"_TN10p"
      open(22,file=ofile)
      write(22, *) "  year  jan   feb   mar   apr   may   jun  ",
     &             " jul   aug   sep   oct   nov   dec annual"
      do i=1,YRS
        write(22,'(i6,13f7.2)') i+SYEAR-1,(tn10out(i,j),j=1,13)
      enddo
      close(22)
      endif

      cnt=0
      wsdi=0.
      csdi=0.
c     if(flg == 1) then
c       wsdi=MISSING
c       csdi=MISSING
c       goto 140
c     endif

      do i=1,YRS
        cntx=0
        cntn=0
        nn=0
        year=i+SYEAR-1
        do month=1,12
          if(leapyear(year) == 1)then
            kth=MONLEAP(month)
          else
            kth=MON(month)
          endif
          do day=1,kth
            if(month /= 2.or.day /= 29) nn=nn+1
            cnt=cnt+1
            if(TMAX(cnt) > thresax90(nn).and.nomiss(TMAX(cnt))) then
              cntx=cntx+1
              if(month == 12.and.day == 31.and.cntx >= 6)
     &                  wsdi(i)=wsdi(i)+cntx
            elseif(cntx >= 6)then
              wsdi(i)=wsdi(i)+cntx
              cntx=0
            else
              cntx=0
            endif
            if(TMIN(cnt) < thresan10(nn).and.nomiss(TMIN(cnt))) then
              cntn=cntn+1
              if(month == 12.and.day == 31.and.cntn >= 6)
     &                  csdi(i)=csdi(i)+cntn
            elseif(cntn >= 6)then
              csdi(i)=csdi(i)+cntn
              cntn=0
            else
              cntn=0
            endif
          enddo  ! day
        enddo    ! month
        if(YNASTAT(i,3) == 1) csdi(i)=MISSING
        if(YNASTAT(i,2) == 1) wsdi(i)=MISSING
      enddo      ! year

      if(flgtx == 0) then
      ofile=trim(ifile)//"_WSDI"
      open(22,file=ofile)
      write(22,*) " year     wsdi"
      do i=1,YRS
        write(22,'(i6,f6.1)') i+SYEAR-1,wsdi(i)
      enddo
      close(22)
      endif

      if(flgtn == 0)then
      ofile=trim(ifile)//"_CSDI"
      open(22,file=ofile)
      write(22,*) " year     csdi"
      do i=1,YRS
        write(22,'(i6,f6.1)') i+SYEAR-1,csdi(i)
      enddo
      close(22)
      endif

      end subroutine tx10p


      subroutine threshold(idata, lev, nl, odata, flg)
      use COMM, only: byrs, ss, winsize
      implicit none
      integer flg,nl,i,j,k
      real, intent(in) :: idata(BYRS,365+2*SS), lev(nl)
      real, intent(out) :: odata(365,nl)

      real tosort(BYRS*WINSIZE),rtmp(nl)
      integer nn
      logical ismiss,nomiss

      do i=1,365
        nn=0
        do j=1,BYRS
          do k=i,i+2*SS
c           if(j == 1.and.k == 1) print*,'##2##',idata(j,k),MISSING
            if(nomiss(idata(j,k))) then
              nn=nn+1
              tosort(nn)=idata(j,k)
            endif
          enddo
        enddo
        if(nn < int(BYRS*WINSIZE*.85)) then
c         print*,"##1##",nn
          flg=1
          return
        endif
        call percentile(tosort,nn,nl,lev,rtmp)
        do j=1,nl
          odata(i,j)=rtmp(j)
        enddo
      enddo

      end subroutine threshold


      pure logical function ismiss(a)
      use COMM, only: missing
      implicit none
      real, intent(in) :: a
      real rmiss

      rmiss=MISSING+1.
      if(a > rmiss) then
        ismiss=.FALSE.
      else
        ismiss=.TRUE.
      endif

      end function ismiss


      pure logical function nomiss(a)
      use COMM, only: missing
      implicit none
      real,intent(in) :: a
      real rmiss
      rmiss=MISSING+1.
      if(a < rmiss) then
        nomiss=.FALSE.
      else
        nomiss=.TRUE.
      endif

      end function nomiss


      impure FUNCTION ran2(idum)
      implicit none
      INTEGER, intent(inout) :: idum
      REAL ran2
      integer, parameter :: IM1=2147483563,IM2=2147483399,IMM1=IM1-1,
     & IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     & NTAB=32,NDIV=1+IMM1/NTAB
      real, parameter :: AM=1./real(IM1),EPS=1.2e-7,RNMX=1.-EPS

     
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum <= 0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum < 0) idum=idum+IM1
          if (j <= NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum < 0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2 < 0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy < 1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)

      END function ran2
