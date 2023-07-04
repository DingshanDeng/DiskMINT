MODULE auxilary
   
   use global_variables
   
   PUBLIC
   
   CONTAINS
   
   !########################################################
   SUBROUTINE GetNline(filename, nline)
   !########################################################
     
      CHARACTER(len=300), INTENT(IN)  :: filename
      INTEGER, INTENT(OUT)            :: nline
   
      CHARACTER(len=300)              :: dummy
      INTEGER                         :: error
      INTEGER                         :: tmp
   
      OPEN(10, file = filename, status = 'old')
   
      nline = 0
      DO
         READ(10, '(a)', iostat=error) dummy
         IF (error /= 0) EXIT

         tmp = index(dummy, '!')

         IF (tmp.NE.0) THEN
            dummy = dummy(1:tmp - 1)
         ENDIF

         IF (dummy.NE.'') THEN
            nline = nline + 1
         ENDIF
      ENDDO
   
      CLOSE(10)
   
      RETURN
    
   END SUBROUTINE GetNline
   !########################################################

   !########################################################
   SUBROUTINE ReadFile(filename, nline, line)
     
      CHARACTER(len=300), INTENT(IN)  :: filename
      INTEGER, INTENT(IN)             :: nline
      CHARACTER(len=300),DIMENSION(nline), INTENT(OUT) :: line
   
      CHARACTER(len=300)              :: dummy
      INTEGER                         :: error
      INTEGER                         :: tmp
      INTEGER                         :: i
   
      OPEN(10, file = filename, status = 'old')
   
      i = 0
      DO
         READ(10, '(a)', iostat=error) dummy
         IF (error /= 0) EXIT

         tmp = INDEX(dummy, '!')

         IF (tmp.NE.0) THEN
            dummy = dummy(1:tmp - 1)
         ENDIF

         IF (dummy.NE.'') THEN
            i = i + 1
            line(i) = dummy
         ENDIF
      ENDDO
   
      IF (i.EQ.nline) THEN
         !PRINT*, 'Succesffuly read : ', trim(filename)
      ELSE
         PRINT*,'Problem encountered while reading : ', filename
         STOP
      END IF
   
      CLOSE(10)
   
      RETURN
    
   END SUBROUTINE ReadFile
   !########################################################
   
   !########################################################
   subroutine cubic_root(a,b,c,x)
   
      ! Find the real roots of the cubic equation:
      ! x^3 + ax^2 + bx + c = 0
      ! Inspired form Numerical recipe
   
      implicit none
   
      ! In/Out
      real(kind=dp), intent(in) :: a, b, c
      real(kind=dp), intent(out) :: x
   
      ! Local
      real(kind=dp) :: Q, R, D1, D2, theta
      real(kind=dp) :: x1,x2,x3
   
      Q = a**2.0_dp - 3.0_dp*b
      Q = Q/9.0_dp
   
      R = 2.0_dp*a**3.0_dp - 9.0_dp*a*b + 27.0_dp*c
      R = R/54.0_dp
   
      if (R**2_dp>Q**3.0_dp) then
         D1 = abs(R) + sqrt(R**2.0_dp - Q**3.0_dp)
         D1 = -sign(1.0_dp,R) * D1**(1.0_dp/3.0_dp)
         if (D1.ne.0.0_dp) then
            D2 = Q/D1
         else
            D2 = 0.0_dp
         end if
         x1 = (D1+D2) - a/3.0_dp
      else
         print*, 'Three real roots...'
         theta = acos(R/sqrt(Q**3.0_dp))
         x1 = -2.0_dp*sqrt(Q)*cos(theta/3.0_dp) - a/3.0_dp
         x2 = -2.0_dp*sqrt(Q)*cos((theta + 2.0_dp*pi)/3.0_dp) - a/3.0_dp
         x3 = -2.0_dp*sqrt(Q)*cos((theta - 2.0_dp*pi)/3.0_dp) - a/3.0_dp
         print*, 'x1 = ', x1
         print*, 'x2 = ', x2
         print*, 'x3 = ', x3
         stop
      end if
   
      x = x1
      
      return
   
   end subroutine cubic_root
   !########################################################
   
   !########################################################
   subroutine get_sorted_index(arr,index)
   !########################################################
   
   implicit none
   REAL(kind=dp), DIMENSION(:), intent(in) :: arr !<[in] the input array for the sorted process
   INTEGER, DIMENSION(:), intent(out) :: index !<[out] the index of the elements, from lowest to highest corresponding values
   INTEGER, parameter :: nn=15, nstack=50
   REAL(kind=dp) :: a
   INTEGER :: n,k,i,j,indext,jstack,l,r
   INTEGER, DIMENSION(nstack) :: istack
   INTEGER :: tmp_index
   
   if (size(index).ne.size(arr)) then
     write(*,*) 'in get_sorted_index. size are different for the two arguments'
     call exit(24)
   endif
   
   n = size(index)
   
   ! initialize list of index
   do i=1,n
     index(i) = i
   enddo
   
   jstack=0
   l=1
   r=n
   do
     if (r-l < nn) then
       do j=l+1,r
         indext=index(j)
         a=arr(indext)
         do i=j-1,l,-1
           if (arr(index(i)) <= a) exit
           index(i+1)=index(i)
         end do
         index(i+1)=indext
       end do
       if (jstack == 0) return
       r=istack(jstack)
       l=istack(jstack-1)
       jstack=jstack-2
     else
       k=(l+r)/2
       
       ! swaping indexes
       tmp_index = index(k)
       index(k) = index(l+1)
       index(l+1) = tmp_index
       
       call icomp_xchg(arr,index(l),index(r))
       call icomp_xchg(arr,index(l+1),index(r))
       call icomp_xchg(arr,index(l),index(l+1))
       i=l+1
       j=r
       indext=index(l+1)
       a=arr(indext)
       do
         do
           i=i+1
           if (arr(index(i)) >= a) exit
         end do
         do
           j=j-1
           if (arr(index(j)) <= a) exit
         end do
         if (j < i) exit
         tmp_index = index(i)
         index(i) = index(j)
         index(j) = tmp_index
       end do
       index(l+1)=index(j)
       index(j)=indext
       jstack=jstack+2
       if (jstack > nstack) then
         write(*, *) 'indexx: nstack too small'
         call exit(24)
       endif
       if (r-i+1 >= j-l) then
         istack(jstack)=r
         istack(jstack-1)=i
         r=j-1
       else
         istack(jstack)=j-1
         istack(jstack-1)=l
         l=i
       end if
     end if
   end do
   end subroutine get_sorted_index
   !########################################################
   
   !########################################################
   subroutine icomp_xchg(arr,i,j)
   !########################################################
   
   REAL(kind=dp), DIMENSION(:), intent(in) :: arr !<[in] the reference array
   INTEGER, intent(inout) :: i,j !<[in,out] index we will swap if necessary
   INTEGER :: swp
   if (arr(j) < arr(i)) then
     swp=i
     i=j
     j=swp
   end if
   end subroutine icomp_xchg
   !########################################################
   
END MODULE auxilary