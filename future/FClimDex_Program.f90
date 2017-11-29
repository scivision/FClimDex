Program FClimDex
use COMM
implicit none

character(80) :: ifile, header,fsite
integer :: stnnum,upara, ulog,argc,i

MISSING=-99.9
SS=int(WINSIZE/2)

argc = command_argument_count()
if (argc<2) error stop "must specify sitefile datafile"

call get_command_argument(1,fsite)
call get_command_argument(2,ifile)

stnnum=1
open(newunit=upara, file=fsite, status='old',action='read')

!      open(newunit=uin, file=trim(path)//"infilename.txt",
!     &     status='old',action='read')

read(upara, '(a80)') header
read(upara, '(a20,f10.2,3i6,i10)') STNID, LATITUDE, STDSPAN, BASESYEAR, BASEEYEAR, PRCPNN
!     print*,'##3##',STDSPAN,BASESYEAR,BASEEYEAR,PRCPNN
!      read(uin, '(a80)', end=100) ifile
if(trim(ifile) == " ") then
!        write(stderr,*) "in: infilename.txt, line:", stnnum
error stop "Read in data filename ERROR "
endif

open(newunit=ulog, file=trim(ifile)//"_log", status='unknown',   action='write')
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
!      goto 77

!100   close(uin)
close(upara)
stnnum=stnnum-1
print *, "Total ",stnnum,"stations calculated"
write(ulog,*) "Total ",stnnum,"stations calculated"
close(ulog)
end program
