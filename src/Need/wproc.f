      subroutine wproc(string)

      character * (*) string
      character * (*) stri
      integer jj

      jj=istrl(string)
      stri=string(1:jj)      

      write(*,*)string,'C'

      write(6,*)'C                                                  C'
      write(6,*)'C      Computing ',stror,
     . ' cross section for            C'    
      write(6,*)'C                                                  C' 
      write(*,*)'C ',stri,'C'
      write(6,*)'C                                                  C'
      write(6,*)'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC' 

      return
      end 
