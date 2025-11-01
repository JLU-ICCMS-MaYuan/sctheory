!****g* src/xpk0
! 
!  NAME
!   xpk0
!  AUTHOR
!   W. Hergert (modification)
!  TYPE
!   Test program
!  PACKAGE
!   Pikaia
!  USAGE
!   Test of genetic algorithm
!  DESCRIPTION
!   * The originalcode from Pikaia distribution is used.
!   * An input subprogram is used with the tokens of the TB fit program.
!     This allows to modify easily the control parameters. 
!  MODIFICATION HISTORY
!   14th July 2022  : modification of the original test
!
! SOURCE
! 
      program xpkaia
      implicit none
      character(LEN=60)              :: problem
      integer                        :: seed, i, status, nparm,nfou
      real                           :: ctrl(12), ctrl0(12), f, twod
      logical                        :: l_final,l_data,l_scale
      real, allocatable,dimension(:) :: x
      external twod    
      
!
!     Set control variables (use defaults)
!
      do  i=1,12
         ctrl0(i) = -1
      enddo 
      ctrl0(2)=50
      ctrl= ctrl0
! 
!     call input routine    to overwrite ctrl   
!
      write(*,*) 'input'
      CALL ga_input(ctrl,problem,seed,nparm,nfou,l_final,l_data,l_scale)
      ALLOCATE(x(nparm))
!
!     First, initialize the random-number generator
!
      call rninit(seed)
!      
!     Now call pikaia
!     write(*,*) 'call pikaia'
      call pikaia(twod,nparm,ctrl,x,f,status)
!
!     Print the results
!    
      write(*,*) problem
      write(*,*) ' status: ',status
      write(*,*) '      x: ',x
      write(*,*) '      f: ',f
      write(*,*) '   ctrl:'
      write(*,20) ctrl
!
!     Export protocol to a file
!  
       WRITE(8,*) '!----------------------------------------------!'
       WRITE(8,*) '!  smooth 2d-landscape (see 5.1 p 35 ff.       !'
       WRITE(8,*) '!----------------------------------------------!'
       WRITE(8,*)
       WRITE(8,21) problem
21     FORMAT(3x,A60)
       WRITE(8,*)       
       WRITE(8,*) '   vector of control data  (standard test) :'
       WRITE(8,*)
       WRITE(8,20) ctrl0
       WRITE(8,*)
       WRITE(8,*) '   vector of control data  ( actually used):'
       WRITE(8,*)
       WRITE(8,20) ctrl
 20    FORMAT(3x,6f12.5/3x,6f12.5)
       WRITE(8,*)
       WRITE(8,*) '   status of calculation : ', status
       WRITE(8,*) '   result of the fit     : ', x
!
!      Output fit
!
       IF(l_final) THEN
          DO i=1,nparm
             WRITE(10,30) i, x(i)
          ENDDO
       ENDIF   
 30    FORMAT(I5,2x,F14.6)               
      end
