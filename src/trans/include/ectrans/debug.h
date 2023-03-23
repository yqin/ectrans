!#define COLLECT_DATAPATTERN
!#define BLOCK_PACKING
!#define BLOCK_UNPACKING
!#define DUMP_STDERR

#ifdef DUMP_STDERR
#define DEBUG_RAW(fun,line,proc,var,val) write(0,'(A,":",I0,"(",I0,"): ",A," = ",I0)') fun,line,proc,var,val
#else
#define DEBUG_RAW(fun,line,proc,var,val)
#endif

SUBROUTINE DEBUG(OPROC,FUNC,LINE,VAR,VAL)
USE TPM_DISTR,ONLY: MYPROC

IMPLICIT NONE

CHARACTER(*),INTENT(IN) :: FUNC,VAR
INTEGER,INTENT(IN) :: OPROC,LINE,VAL

! 0 to debug all ranks, n to debug a specific rank, -1 to disable dubugging
IF (OPROC == MYPROC .OR. OPROC == 0) THEN
  DEBUG_RAW(FUNC,LINE,MYPROC,VAR,VAL)
ENDIF

END SUBROUTINE DEBUG
