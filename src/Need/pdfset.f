      subroutine pdfset
      implicit none
      include 'nlooprun.f'
      include 'pdfiset.f'
      double precision amz
      common/couple/amz

      character *50 prefix
      character *36 pdfstring
      integer nset
      common/prefix/nset,prefix
      common/pdfstring/pdfstring


      if     (iset.eq.61) then
      amz=0.1197d0 
      nlooprun=2
      pdfstring='MRST2002 NLO'
      elseif     (iset.eq.62) then
      amz=0.1154d0
      nlooprun=3
      pdfstring='MRST2002 NNLO'
      elseif     (iset.eq.49) then    
      amz=0.130d0
      nlooprun=1    
      pdfstring='MRST2002 LO'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      elseif     (iset.eq.71) then
      amz=0.1205d0
      nlooprun=2
      pdfstring='MRST2004 NLO'
      elseif     (iset.eq.72) then
      amz=0.1167d0
      nlooprun=3
      pdfstring='MRST2004 NNLO'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      elseif     (iset.eq.41) then
      amz=0.119d0
      nlooprun=2
      pdfstring='MRST2001 NLO'
      elseif     (iset.eq.42) then
      amz=0.117d0
      nlooprun=2
      pdfstring='MRST2001 NLO lower alphas'
      elseif     (iset.eq.43) then
      amz=0.121d0
      nlooprun=2
      pdfstring='MRST2001 NLO higher alphas'
      elseif     (iset.eq.44) then
      amz=0.121d0
      nlooprun=2
      pdfstring='MRST2001 NLO better fit to jet data'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      elseif     (iset.eq.45) then
      amz=0.1155d0
      nlooprun=3
      pdfstring='MRST2001 NNLO'
      elseif     (iset.eq.46) then
      amz=0.1155d0
      nlooprun=3
      pdfstring='MRST2001 NNLO fast evolution'
      elseif     (iset.eq.47) then
      amz=0.1155d0
      nlooprun=3
      pdfstring='MRST2001 NNLO slow evolution'
      elseif     (iset.eq.48) then
      amz=0.1180d0
      nlooprun=3
      pdfstring='MRST2001 NNLO better fit to jet data'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      elseif (iset.eq.2) then !CTEQ4 NLO
      amz=0.116d0
      nlooprun=2
      pdfstring='CTEQ4 NLO'
      elseif (iset.eq.1) then !CTEQ4 LO
      amz=0.132d0
      nlooprun=1
      pdfstring='CTEQ4 LO'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      elseif (iset.eq.11) then !MRS98 NLO central gluon
      amz=0.1175
      nlooprun=2
      pdfstring='MRST98 NLO'
      elseif (iset.eq.12) then !MRS98 NLO higher gluon
      amz=0.1175
      nlooprun=2
      pdfstring='MRST98 NLO higher gluon'
      elseif (iset.eq.13) then !MRS98 NLO lower gluon
      amz=0.1175
      nlooprun=2
      pdfstring='MRST98 NLO lower gluon'
      elseif (iset.eq.14) then !MRS98 NLO lower as
      amz=0.1125
      nlooprun=2
      pdfstring='MRST98 NLO lower alphas'
      elseif (iset.eq.15) then !MRS98 NLO higher as
      amz=0.1225
      nlooprun=2
      pdfstring='MRST98 NLO higher alphas'
      elseif (iset.eq.16) then !MRST98 LO central gluon
      amz=0.125
      nlooprun=1
      pdfstring='MRST98 LO'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      elseif (iset.eq.21) then
      Call SetCtq5(1)
      amz=0.118d0
      nlooprun=2
      pdfstring='CTEQ5M NLO'
      elseif (iset.eq.22) then
      Call SetCtq5(2)
      amz=0.118d0
      nlooprun=2
      pdfstring='CTEQ5D NLO DIS'
      elseif (iset .eq.23) then
      Call SetCtq5(3)
      amz=0.127d0
      nlooprun=1
      pdfstring='CTEQ5L LO'
      elseif (iset .eq.24) then
      Call SetCtq5(4)
      amz=0.118d0
      nlooprun=2
      pdfstring='CTEQ5HJ NLO large x gluon enhanced'
      elseif (iset .eq.25) then
      Call SetCtq5(5)
      amz=0.118d0
      nlooprun=2
      pdfstring='CTEQ5HQ NLO Heavy quark'
      elseif (iset .eq.28) then
      Call SetCtq5(8)
      amz=0.118d0
      nlooprun=2
      pdfstring='CTEQ5M1 NLO improved'
      elseif (iset .eq.29) then
      Call SetCtq5(9)
      amz=0.118d0
      nlooprun=2
      pdfstring='CTEQ5HQ1 NLO improved'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      elseif (iset.eq.30) then !MRS99 NLO central gluon
      amz=0.1175
      nlooprun=2
      pdfstring='MRST99 NLO'
      elseif (iset.eq.31) then !MRS99 NLO higher gluon
      amz=0.1175
      nlooprun=2
      pdfstring='MRST99 higher gluon'
      elseif (iset.eq.32) then !MRS99 NLO lower gluon
      amz=0.1175
      nlooprun=2
      pdfstring='MRST99 lower gluon'
      elseif (iset.eq.33) then !MRS99 NLO lower as
      amz=0.1125 
      nlooprun=2
      pdfstring='MRST99 lower alphas'
      elseif (iset.eq.34) then !MRS99 NLO higher as
      amz=0.1225
      nlooprun=2
      pdfstring='MRST99 higher alphas'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      elseif (iset .eq.53) then
      amz=0.118d0
      Call SetCtq6(1)
      nlooprun=2
      pdfstring='CTEQ6M NLO'
      elseif (iset.eq.51) then
      amz=0.118d0
      Call SetCtq6(3)
      nlooprun=1
      pdfstring='CTEQ6L LO'
      elseif (iset.eq.52) then
      amz=0.130d0
      Call SetCtq6(4)
      nlooprun=1
      pdfstring='CTEQ6L1 LO'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC                                                     CCCCCCC
CCCCC                NEW: MSTW2008                        CCCCCCC
      
      elseif (iset.eq.90) then
      amz=0.13939d0
      nlooprun=1
      pdfstring='MSTW2008 LO'
      prefix = "Pdfdata/mstw2008/mstw2008lo"
      elseif (iset.eq.91) then
      amz=0.12018
      nlooprun=2
      pdfstring='MSTW2008 NLO'
      prefix = "Pdfdata/mstw2008/mstw2008nlo"
      elseif (iset.eq.92) then
      amz=0.11707
      nlooprun=3
      pdfstring='MSTW2008 NNLO'
      prefix = "Pdfdata/mstw2008/mstw2008nnlo"
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      else
        write(6,*) 'Unimplemented distribution= ',iset

      stop
      endif      
      return
      end
 

