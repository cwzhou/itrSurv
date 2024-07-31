!itrSurv.f90 in testpackage/itrSurv/src
MODULE INNERS
  IMPLICIT NONE

  PUBLIC

  INTEGER, PARAMETER :: dp = selected_real_kind(15,307)

  INTEGER, SAVE :: ERT ! 0/1 1 = use extremely randomized tree
  INTEGER, SAVE :: minEvent ! minimum number of events in a node
  INTEGER, SAVE :: mTry ! maximum number of covariates to try
  INTEGER, SAVE :: n  ! number of sampled cases
  INTEGER, SAVE :: nAll ! number of cases
  INTEGER, SAVE :: nLevs ! maximum number of levels
  INTEGER, SAVE :: nodeSize ! minimum number of cases in a node
  INTEGER, SAVE :: np ! number of covariates
  INTEGER, SAVE :: nrNodes ! maximum number of nodes in a tree
  INTEGER, SAVE :: nt ! number of time points
  INTEGER, SAVE :: nTree ! number of trees
  INTEGER, SAVE :: replace
  INTEGER, SAVE :: rule ! 1 if logrank 0 if truncated mean
  INTEGER, SAVE :: sampleSize
  ! for survival probability, the index of nearest time point
  INTEGER, SAVE :: sIndex
  INTEGER, SAVE :: uniformSplit ! 0/1 1 = random cutoff comes from values

  ! censoring indicator 1 = not censored
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: deltaAll
  ! censoring indicator 1 = not censored
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: deltaAll_m
  ! censoring indicator 1 = not censored
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: delta
  ! censoring indicator 1 = not censored for cause m for CR 
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: delta_m
 ! number of categories in each np covariate
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: nCat

  ! the probability for a random split
  REAL(dp), SAVE :: rs
  ! the fraction above sIndex time that survival time lies
  REAL(dp), SAVE :: sFraction
  ! the stratified random split coefficient
  REAL(dp), SAVE :: stratifiedSplit

  ! time differences
  REAL(dp), DIMENSION(:), ALLOCATABLE, SAVE :: dt
  ! probability mass vector of survival function for sampled cases
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: pr
  ! probability mass vector of survival function for all cases
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: prAll
  ! covariates to be considered for split for sampled cases
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: x
  ! covariates to be considered for split for all cases
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: xAll

  ! TRUE = time differences vector has been allocated
  LOGICAL, SAVE :: dtAllocated = .FALSE.
  ! TRUE = all other allocatables have been allocated
  LOGICAL, SAVE :: isAllocated = .FALSE.
  ! TRUE = using survival probability
  LOGICAL, SAVE :: isSurvival = .FALSE.
  LOGICAL, SAVE :: isCIF = .FALSE.
  LOGICAL, SAVE :: isPhase1 = .FALSE.
  LOGICAL, SAVE :: isPhase2CR = .FALSE.
  LOGICAL, SAVE :: isPhase2 = .FALSE.

  TYPE NODE
    INTEGER :: nNode
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Func
    REAL(dp), DIMENSION(:), ALLOCATABLE :: mean
    REAL(dp), DIMENSION(:), ALLOCATABLE :: Prob
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: matrix
  END TYPE

  TYPE(NODE), DIMENSION(:), ALLOCATABLE, SAVE :: trees

  TYPE ForestValues
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Func
    REAL(dp), DIMENSION(:), ALLOCATABLE :: mean
    REAL(dp), DIMENSION(:), ALLOCATABLE :: Prob
  END TYPE

  TYPE(ForestValues), SAVE :: forest
  ! =================================================================================

  CONTAINS

! sample an array of indices allowing for duplicates
!   nCases: integer, the total number of indices to sample
!   n: integer, the number of indices to draw
! returns an array (1:n) of the indices sampled
FUNCTION sampleWithReplace(nCases, n) RESULT(array)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nCases
  INTEGER, INTENT(IN) :: n

  INTEGER, DIMENSION(1:n) :: array

  INTEGER :: i

  REAL(dp) :: rnd

  EXTERNAL :: rnd

  DO i = 1, n
    array(i) = 1 + floor(rnd(0.d0, 1.d0)*nCases)
  END DO

END FUNCTION sampleWithReplace
! =================================================================================

! sample an array of indices without allowing duplicates
!   nCases: integer, the total number of indices to sample
!   n: integer, the number of indices to draw
! returns an array (1:n) of the indices sampled
FUNCTION sampleWithOutReplace(nCases, n) RESULT(array)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nCases
  INTEGER, INTENT(IN) :: n

  INTEGER, DIMENSION(1:n) :: array

  INTEGER :: j, cnt

  LOGICAL, DIMENSION(1:nCases) :: used

  REAL(dp) :: rnd

  EXTERNAL :: rnd

  used(:) = .FALSE.

  cnt = 1
  DO WHILE (cnt <= n)
    ! draw a number
    j = 1 + floor(rnd(0.d0, 1.d0)*nCases)
    ! if the number has already been drawn cycle
    if (used(j)) CYCLE
    ! if the number has not previously been drawn, store and set used
    array(cnt) = j
    used(j) = .TRUE.
    cnt = cnt + 1
  END DO

END FUNCTION sampleWithOutReplace
! =================================================================================

! Identify the optimal split
!   nCases : integer, the number of elements in input casesIn
!   casesIn : integer(:), the indices of the cases in this node
!   nv : integer, the number of covariates to include in search
!   varsIn : integer(:), covariates to include in search
!   splitVar : integer, the index of the selected variable for splitting
!   cutoffBest : real(:), the cutoff (<= go to 'left' node)
!   splitFound : integer, 0 = no split; 1 = found a split
!   casesOut : integer(:), elements of casesIn that go left; ind if yes, 0
!     otherwise
!   nCuts : integer, the number of cutoff values returned
!   lft : integer, the number of cases in the left node
! CALL tfindSplit(size(ind), ind, size(pind), pind, splitVar, cutoffBest, splitFound, indOut, nc, lft)
!
SUBROUTINE tfindSplit(nCases, casesIn, nv, varsIn, &
                    & splitVar, cutoffBest, splitFound, casesOut, nCuts, lft)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nCases
  INTEGER, DIMENSION(1:nCases), INTENT(IN) :: casesIn
  INTEGER, INTENT(IN) :: nv
  INTEGER, DIMENSION(1:nv), INTENT(IN) :: varsIn
  INTEGER, INTENT(OUT) :: splitVar
  REAL(dp), DIMENSION(1:nLevs), INTENT(OUT) :: cutoffBest
  INTEGER, INTENT(OUT) :: splitFound
  INTEGER, DIMENSION(1:nCases), INTENT(OUT) :: casesOut
  INTEGER, INTENT(OUT) :: nCuts
  INTEGER, INTENT(OUT) :: lft

  INTEGER :: cnt, cnt_m, i, ikv, j, jj, k, kv, l, nUncensored, nUncensored_m, ptr, rightNode, rightNode_m
  INTEGER :: rUnifSet, rUnifSet_m, set, splitLeft, splitLeftFinal, tieCovariate
  INTEGER :: splitLeft_m, splitLeftFinal_m
  INTEGER :: tieValue, variablesTried
  INTEGER, DIMENSION(1:nCases) :: cases, dSorted, dSorted_m, tcases
  INTEGER, DIMENSION(1:nv) :: variables
  INTEGER, DIMENSION(:), ALLOCATABLE :: ind, ind_m, indSingles, indSingles_m, leftCases, leftCases_m, rightCases, rightCases_m
  INTEGER, DIMENSION(:), ALLOCATABLE :: uncensoredIndices, uncensoredIndices_m

  REAL(dp) :: cutoff, maxValueSplit, maxValueXm, rUnif, rUnif_m, valuej, tester3a, tester3b
  REAL(dp), DIMENSION(1:nt) :: atRiskLeft, atRiskRight, D, D_m, denJ
  REAL(dp), DIMENSION(1:nt) :: eventsLeft, eventsLeft_m, eventsRight, eventsRight_m, numJ, pd1, pd2, Rcum, pd1_m, pd2_m
  REAL(dp), DIMENSION(1:nCases) :: xSorted
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: prl, prr, prl_m, prr_m

  LOGICAL :: randomSplit
  LOGICAL, DIMENSION(:), ALLOCATABLE :: singles, singles_m

  LOGICAL :: are_equal
  INTEGER :: index_delta

  LOGICAL :: notEqual


  REAL(dp) :: rnd, random1

  EXTERNAL :: rnd

  are_equal = .TRUE.

  ! write(*,'(/,A)') '============================ tfindSplit ============================'
  ! PRINT *, "******************** tfindSplit ********************"
  ! determine if this is to be a random split
  randomSplit = rnd(0.d0, 1.d0) <= rs
  ! PRINT *, "rs =", rs

  ! splitFound is flag indicating if a split was found
  ! 0 = no split found; 1 = split found
  splitFound = 0

  !! initialize tracking variables

  ! index of variable with largest critical value
  ! negative number indicates no split found
  splitVar = -1

  ! largest critical value
  maxValueSplit = 0.0

  ! cutoff that gives largest critical value or is randomly selected
  cutoffBest = 0.0

  ! set casesOut to 0
  casesOut = 0

  ! set number of cutoffs to 0
  nCuts = 0

  ! indices of variables to be explored
  variables = varsIn

  ! location of last available parameter in variables
  ptr = nv

  ! tracks the number of variables tried
  variablesTried = 0

  tcases = (/(i,i=1,nCases)/)
  ! PRINT *, "tcases"
  ! PRINT *, tcases

  IF (nCases < 2*nodeSize) RETURN

  DO i = 1, nv

    ! if mTry successful splits already explored, exit
    IF (variablesTried .EQ. mTry) EXIT

    ! randomly select a covariate on which to split
    ikv = 1 + floor(rnd(0.d0, 1.d0)*ptr)
    kv = variables(ikv)

    ! move this index to the end of the list; shift pointer down
    j = variables(ptr)
    variables(ptr) = kv
    variables(ikv) = j

    ptr = ptr - 1

    ! pull appropriate covariate. If un-ordered factor use mean survival/cif time
    ! for each factor level
    IF (nCat(kv) .GT. 1) THEN
      ! PRINT *, "nCat(kv) > 1 so using mean survival/cif time for each factor level"
      CALL getCovariate(nCases, casesIn, kv, xSorted)
    ELSE
      ! PRINT *, "nCat(kv) is either 0 = continuous or 1 = ordered factors"
      xSorted = x(casesIn,kv)
    END IF

    ! PRINT *, "xSorted"
    ! PRINT *, SIZE(xSorted)
    ! PRINT *, xSorted

    cases = casesIn

    ! sort the covariate and track the indices
    CALL qsort4(xSorted, cases, 1, nCases)
!    CALL hpsort_eps_epw(nCases, xSorted, cases, 1d-8)

    ! sort event indicator data accordingly
    dSorted = delta(cases)
    dSorted_m = delta_m(cases)

    ! ******************** splitBoundaries ********************
    ! identify minimum cases for left and right splits based on minimum uncensored cases, minimum node size, and assurance that all equal valued cases are included in the minimum nodes
    ! PRINT *, "******************** splitBoundaries ********************"

    rUnif = 0.d0
    rUnif_m = 0.d0
    rUnifSet = -1
    rUnifSet_m = -1
    ! PRINT *, "test1"

    ! Do below for both isPhase1 and isPhase2CR
    ! cases that are not-censored
    IF (isPhase1) THEN
    	! PRINT *, "PHASE 1"
    	uncensoredIndices = pack(tcases, dSorted .EQ. 1)
    	! PRINT *, "uncensoredIndices: ", uncensoredIndices
    	nUncensored = size(uncensoredIndices)
    	! PRINT *, "nUncensored: ", nUncensored
    END IF

    IF (isPhase2CR) THEN
    	! PRINT *, "PHASE 2 CR"
    	uncensoredIndices = pack(tcases, dSorted_m .EQ. 1)
    	! PRINT *, "uncensoredIndices: ", uncensoredIndices
    	nUncensored = size(uncensoredIndices)
    	! PRINT *, "nUncensored: ", nUncensored
    END IF

    ! PRINT *, "test2a"
    ! PRINT *, "minEvent: ", minEvent 

    ! if too few cases to meet minimum number of uncensored cases, CYCLE
    IF (nUncensored .LT. (minEvent * 2)) CYCLE
   
    !! able to split and satisfy minimum number of events in each node

    ! cases to left include all indices up to and including minEvent case
    ! must have at least nodeSize cases
    splitLeft = max(uncensoredIndices(minEvent), nodeSize)
    ! PRINT *, "test3a"

    ! move splitLeft up to include cases with equivalent values of x
    splitLeft = count(xSorted .LE. (xSorted(splitLeft) + 1e-8))
    ! PRINT *, "test4a"

    ! cases to right
    ! include all indices down to and including nUncensored - minEvent + 1 case
    ! must have at least nodeSize cases
    rightNode = min(uncensoredIndices(nUncensored - minEvent + 1), &
                  & nCases - nodeSize + 1)
    ! PRINT *, "test4.5a"

    ! move rightNode down to include cases with equivalent values of x
    ! splitLeftFinal is the last possible case for the left node
    splitLeftFinal = count(xSorted .LT. xSorted(rightNode))
    
    ! if the splitLeft index is above the splitLeftFinal index cycle,
    ! split is not possible
    IF (splitLeft .GT. splitLeftFinal) CYCLE
    ! PRINT *, "test6a"

  ! IF (isPhase2CR) THEN
      ! cases that are not-censored
      uncensoredIndices_m = pack(tcases, dSorted_m .EQ. 1)
      nUncensored_m = size(uncensoredIndices_m)
      ! PRINT *, "test2b"

      ! if too few cases to meet minimum number of uncensored cases, CYCLE
      IF (nUncensored_m .LT. (minEvent * 2)) CYCLE
     
      !! able to split and satisfy minimum number of events in each node

      ! cases to left include all indices up to and including minEvent case
      ! must have at least nodeSize cases
      splitLeft_m = max(uncensoredIndices_m(minEvent), nodeSize)
      ! PRINT *, "test3b"

      ! move splitLeft up to include cases with equivalent values of x
      splitLeft_m = count(xSorted .LE. (xSorted(splitLeft_m) + 1e-8))
      ! PRINT *, "test4b"

      ! cases to right
      ! include all indices down to and including nUncensored - minEvent + 1 case
      ! must have at least nodeSize cases
      rightNode_m = min(uncensoredIndices_m(nUncensored_m - minEvent + 1), &
                    & nCases - nodeSize + 1)
      ! PRINT *, "test5b"

      ! move rightNode down to include cases with equivalent values of x
      ! splitLeftFinal is the last possible case for the left node
      splitLeftFinal_m = count(xSorted .LT. xSorted(rightNode_m))

      ! if the splitLeft index is above the splitLeftFinal index cycle,
      ! split is not possible
      IF (splitLeft_m .GT. splitLeftFinal_m) CYCLE

    ! END IF

    ! PRINT *, "test7"
    rUnifSet = 0
    rUnifSet_m = 0
    IF ((.NOT. randomSplit) .AND. (ERT .EQ. 1)) THEN

      !************* getUniformSplit *****************
      ! PRINT *, "************* getUniformSplit *****************"

      !! split based on extremely randomized trees and uniform split inputs

      ! if ERT splitLeft = splitLeftFinal, which is the splitting point
      ! and cutoff is set to random value or mid-point depending on uniERT

      !! extremely randomized tree methods used

      IF (uniformSplit .EQ. 0) THEN
      	! PRINT *, "TESTER_UNIF1"

        ! if the cutoff is not determined from a uniform distribution
        ! randomly sample available indices to identify the last case
        ! of the left split to define the minimum; cutoff is the
        ! mid-point {x[r] + x[r+1]} / 2

        ! only indices for which x[r] < x[r+1] can be sampled
        ind = (/ (i, i = splitLeft, splitLeftFinal) /)
        singles = xSorted(ind) .LT. xSorted(ind + 1)
        indSingles = pack(ind, singles)
        cnt = size(indSingles)

        splitLeftFinal = indSingles(1 + floor(rnd(0.d0, 1.d0)*cnt))
        ! the last required case in the left node is now splitLeftFinal
        splitLeft = splitLeftFinal
        rUnif = (xSorted(splitLeftFinal) + xSorted(splitLeftFinal+1))/2.0
        rUnifSet = 1

        ! IF (isPhase2CR) THEN
          ind_m = (/ (i, i = splitLeft_m, splitLeftFinal_m) /)
          singles_m = xSorted(ind_m) .LT. xSorted(ind_m + 1)
          indSingles_m = pack(ind_m, singles_m)
          cnt_m = size(indSingles_m)

          splitLeftFinal_m = indSingles_m(1 + floor(rnd(0.d0, 1.d0)*cnt_m))
          splitLeft_m = splitLeftFinal_m
          rUnif_m = (xSorted(splitLeftFinal_m) + xSorted(splitLeftFinal_m+1))/2.0
          rUnifSet_m = 1
        ! END IF

      ELSE IF (uniformSplit .EQ. 1) THEN
      	! PRINT *, "TESTER_UNIF2"
        ! randomly select a value in the range of values that satisfy the allowed cases in the left/right nodes
        random1 = rnd(0.d0, 1.d0)
        rUnif = random1 * (xSorted(splitLeftFinal+1) - &
              & xSorted(splitLeft)) + xSorted(splitLeft)
        rUnifSet = 1

        ! identify the first case that splits to the right of this value
        ! the preceding case is the last case to the left node
        splitLeftFinal = nCases - count(xSorted > rUnif)

        ! the last required case in the left node is now splitLeftFinal
        splitLeft = splitLeftFinal

        ! IF (isPhase2CR) THEN
          rUnif_m = random1 * (xSorted(splitLeftFinal_m+1) - &
                & xSorted(splitLeft_m)) + xSorted(splitLeft_m)
          rUnifSet_m = 1
          splitLeftFinal_m = nCases - count(xSorted > rUnif_m)
          splitLeft_m = splitLeftFinal_m
        ! END IF

         ! tester3a = rUnif
         ! tester3b = rUnif_m
         ! PRINT *, "TESTER CHECKING: splitLeft_m: ", splitLeft_m, "splitLeft: ", splitLeft
      	 ! PRINT *, "tester3a:", tester3a, "tester3b", tester3b
         ! IF (tester3a /= tester3b) THEN
			! PRINT *, "ERROR: TESTER: SPLITLEFT NEQ SPLITLEFT_M"
			! PRINT *, tester3b
			! PRINT *, tester3a
		! END IF

      END IF

    END IF

    ! PRINT *, "TESTER_UNIF-END"

    ! -1 is returned if cannot satisfy minimum requirements for nodes
    ! cycle to next covariate
    ! PRINT *, "rUnifSet:", rUnifSet
    IF (rUnifSet .EQ. -1) CYCLE

    ! increment the number of covariates that have been explored
    variablesTried = variablesTried + 1

    ! write(*,'(A)') '********************************* maxValue *********************************'
    !***************** maxValue ***************

    ! set initial values for outputs
    set = 0
    maxValueXm = 0.d0
    cutOff = 0.d0

    leftCases = cases(1:(splitLeft-1))
    rightCases = cases(splitLeft:nCases)
    prl = pr(leftCases,:)
    prr = pr(rightCases,:)
    eventsLeft = sum(prl * &
                   & spread(dSorted(1:(splitLeft-1)), 2, nt), DIM = 1)
    eventsRight = sum(prr * &
                    & spread(dSorted(splitLeft:nCases), 2, nt), DIM = 1)

    ! PRINT *, "tt1"
    ! IF (isPhase2CR) THEN 
      leftCases_m = cases(1:(splitLeft_m-1))
      rightCases_m = cases(splitLeft_m:nCases)
      prl_m = pr(leftCases_m,:)
      prr_m = pr(rightCases_m,:)

      eventsLeft_m = sum(prl_m * &
                   & spread(dSorted_m(1:(splitLeft_m-1)), 2, nt), DIM = 1)
      eventsRight_m = sum(prr_m * &
                    & spread(dSorted_m(splitLeft_m:nCases), 2, nt), DIM = 1)
    ! END IF

    IF (isPhase1) THEN
		IF (ANY(leftCases_m /= leftCases)) THEN
		    PRINT *, "LEFTCASES NEQ LEFTCASES_M"
		END IF
		IF (ANY(rightCases_m /= rightCases)) THEN
		    PRINT *, "rightCASES NEQ rightCASES_M"
		END IF
		IF (ANY(dSorted_m /= dSorted)) THEN
		    PRINT *, "dSorted NEQ dSorted_M"
		END IF
	END IF

    ! PRINT *, "tt2"
    ! looking at "at risk" set now (we want overall failure for both Step1 and Step2CR)
    ! PRINT *, "prl: ", prl
    ! PRINT *, "prr: ", prr
    pd1 = sum(prl, DIM = 1) ! for group 1
    pd2 = sum(prr, DIM = 1) ! for group 2
    ! PRINT *, "group1: pd1 = sum(prl, DIM = 1): ", pd1
    ! PRINT *, "group2: pd2 = sum(prr, DIM = 1) ", pd2
    ! pd1_m = sum(prl_m, DIM = 1) ! for group 1
    ! pd2_m = sum(prr_m, DIM = 1) ! for group 2

    ! PRINT *, "tt3"
    atRiskLeft(1) = splitLeft - 1
    atRiskRight(1) = nCases - splitLeft + 1
    ! PRINT *, "tt4"

    DO j = 2, nt
      atRiskLeft(j) = atRiskLeft(j-1) - pd1(j-1)
      atRiskRight(j) = atRiskRight(j-1) - pd2(j-1)
    END DO
    ! PRINT *, "tt5"

    ! if logrank, do calculations that do not depend on node occupancy
    IF (rule == 1) THEN
      CALL logrankSetUp(atRiskLeft, atRiskRight, eventsLeft, eventsRight, &
                      & numJ, denJ)
    END IF

  DO index_delta = 1, nCases
    IF (delta(index_delta) /= delta_m(index_delta)) THEN
      are_equal = .FALSE.
    END IF
  END DO

  ! IF (are_equal) THEN
  !  PRINT *, "tfindSplit: Arrays delta and delta_m are equal."
  ! ELSE
  !  PRINT *, "tfindSplit: Arrays delta and delta_m are not equal."
  ! END IF

    cnt = 1
    cnt_m = 1
    DO j = splitLeft, splitLeftFinal
      ! PRINT *, "j: ", j, "/ncnt: ", cnt, "/ncnt_m: ", cnt_m

      ! at risk indicators for jth case
      ! number of events for jth case
      pd1 = prr(cnt,:) ! at risk set for group 1 for overall survival
      ! PRINT *, "pd1: ", pd1
      
      cnt = cnt + 1
      cnt_m = cnt_m + 1
      
      D = pd1*delta(cases(j))
      ! PRINT *, "D: ", D
      D_m = pd1*delta_m(cases(j))
      ! PRINT *, "D_m: ", D_m
      
      Rcum(1) = 0.0
      Rcum(2) = pd1(1)
      ! PRINT *, "Rcum for timepoints 1 and 2: ", Rcum
      ! from time point 3 to nt
      DO k = 3, nt
        IF (pd1(k-1) .GT. 1d-8) THEN
          Rcum(k) = Rcum(k-1) + pd1(k-1)
        ELSE
          Rcum(k) = Rcum(k-1)
        END IF
      END DO

      ! PRINT *, "Rcum:", Rcum

      Rcum = 1.d0 - Rcum
      ! number at risk
      ! PRINT *, "Rcum = 1.d0 - Rcum:", Rcum

      ! add the jth case to the left node
      atRiskLeft = atRiskLeft + Rcum

      ! remove the jth case from the right node
      atRiskRight = atRiskRight - Rcum

      ! write(*,'(/,/,A)'), '######### Line 558'
      ! number of events
      ! PRINT *, "=======TEST_R1: eventsRight:", eventsRight
      ! PRINT *, "=======TEST_R1: eventsRight for cause m:", eventsRight_m
      ! PRINT *, "=======TEST_L1: eventsLeft:", eventsLeft
      ! PRINT *, "=======TEST_L1: eventsLeft for cause m:", eventsLeft_m
      ! PRINT *, "D: ", D
      ! PRINT *, "D_m: ", D_m

      ! PRINT *, "! add the jth case to the left node"
      ! add the jth case to the left node
      eventsLeft = eventsLeft + D
      eventsLeft_m = eventsLeft_m + D_m

      ! PRINT *, '! remove the jth case from the right node'
      ! remove the jth case from the right node
      eventsRight = eventsRight - D
      eventsRight_m = eventsRight_m - D_m

      ! if the case is not the last case with this covariate value, cycle
      IF (xSorted(j) .GE. (xSorted(j+1) - 1d-8)) CYCLE

      ! PRINT *, "atRiskLeft:", atRiskLeft
      ! PRINT *, "atRiskRight:", atRiskRight
      ! PRINT *, "=======TEST_R2: eventsRight:", eventsRight
      ! PRINT *, "=======TEST_R2: eventsRight for cause m:", eventsRight_m
      ! PRINT *, "=======TEST_L2: eventsLeft:", eventsLeft
      ! PRINT *, "=======TEST_L2: eventsLeft for cause m:", eventsLeft_m
      ! PRINT *, "numJ:", numJ
      ! PRINT *, "denJ:", denJ
      ! write(*,'(/,/)')

      ! PRINT *, "RULE IS:", rule

      ! calculate test statistic
      IF (rule == 1) THEN
        ! PRINT *, "&&&&&&& logrank test &&&&&&&&&"
        CALL logrank(atRiskLeft, atRiskRight, eventsLeft, numJ, &
                   & denJ, valuej)
        ! PRINT *, "logrank statistic valuej = ", valuej
      ELSE IF (rule == 2) THEN
        ! PRINT *, "&&&&&&& truncated mean test &&&&&&&&&"
        CALL meanSplit(atRiskLeft, atRiskRight, eventsLeft, eventsRight, valuej)
      ELSE IF (rule == 3) THEN
        ! PRINT *, "&&&&&&& gray's test &&&&&&&&&"
        CALL Gray_m(nt, nt, atRiskLeft, eventsLeft_m, &
        & atRiskRight, eventsRight_m, valuej)
      END IF
      ! PRINT *, "test statistic: ", valuej

      IF ((set .EQ. 0) .OR. (valuej .GT. maxValueXm)) THEN
        ! PRINT *, "&&&&&&& TEST 5 &&&&&&&&&"

        ! if first value or value > current max, save
        IF (rUnifSet .EQ. 1) THEN
          ! PRINT *, "rUnifSet: ", rUnifSet
          ! PRINT *, "rUnifSet_m: ", rUnifSet_m
          ! PRINT *, "rUnif: ", rUnif
          ! PRINT *, "rUnif_m: ", rUnif_m
          cutoff = rUnif
        ELSE
          cutoff = (xSorted(j) + xSorted(j+1))/2.d0
        END IF
        ! PRINT *, "cutoff: ", cutoff
        maxValueXm = valuej
        tieValue = 1
        set = 1

      ELSE IF (valuej > (maxValueXm - 1d-8)) THEN
        ! PRINT *, "&&&&&&& TEST 6 &&&&&&&&&"
        ! if value is a tie, randomly determine if cutoff should be taken
        tieValue = tieValue + 1

        IF (rnd(0.d0, 1.d0) < (1.d0 / REAL(tieValue))) THEN
          cutoff = (xSorted(j) + xSorted(j+1))/2.d0
        END IF

      END IF

    END DO

    ! if not successful, cycle to next covariate
    ! this condition should never be true
    IF (set .EQ. 0) CYCLE

    ! if successful, determine if it yields the maximum value of the
    ! covariates considered
    IF ((splitVar .EQ. -1) .OR. (maxValueXm .GT. maxValueSplit)) THEN
      ! PRINT *, "SuCCESSful."

      ! if first non-zero or largest value, keep cutoff and value and
      ! reset tie counter to 1
      splitVar = kv
      maxValueSplit = maxValueXm
      tieCovariate = 1

      ! count the number of cases in the left node
      lft = count(xSorted .LE. cutoff)
      casesOut = cases

      cutoffBest = 0

      IF (nCat(kv) .LE. 1) THEN
        cutoffBest(1) = cutoff
        nCuts = 1
      ELSE
        ! for factors, identify factor values contained in left cases
        l = 1
        DO j = 1, nCat(kv)
          DO jj = 1, lft
            IF (nint(x(casesOut(jj),kv)) .EQ. j) THEN
              cutoffBest(l) = j
              l = l + 1
              EXIT
            END IF
          END DO
        END DO
        nCuts = l - 1
      END IF

      ! if a random split was triggered, break out of loop over covariates
      IF (randomSplit) EXIT

    ELSE IF (maxValueXm .GT. (maxValueSplit-1d-8)) THEN
      ! PRINT *, "&&&&&&& TEST 7 &&&&&&&&&"

      ! if equal to current maximum value, increment tie counter and randomly
      ! select the cutoff with equal probability for each tie
      tieCovariate = tieCovariate + 1

      IF (rnd(0.d0, 1.d0) < (1.0 / tieCovariate)) THEN
        splitVar = kv
        maxValueSplit = maxValueXm
        ! count the number of cases in the left node
        lft = count(xSorted .LE. cutoff)
        casesOut = cases

        cutoffBest = 0

        IF (nCat(kv) .LE. 1) THEN
          cutoffBest(1) = cutoff
          nCuts = 1
        ELSE
          ! for factors, identify factor values contained in left cases
          l = 1
          DO j = 1, nCat(kv)
            DO jj = 1, lft
              IF (nint(x(casesOut(jj),kv)) .EQ. j) THEN
                cutoffBest(l) = j
                l = l + 1
                EXIT
              END IF
            END DO
          END DO
          nCuts = l - 1
        END IF
      END IF
    END IF

  END DO

  ! if no split was possible return
  if (splitVar .EQ. -1) RETURN

  ! if successful at finding a split set flag and return
  splitFound = 1
  ! PRINT *, "splitFound: ", splitFound
  ! PRINT *, "splitVar: ", splitVar
  ! PRINT *, "cutoffBest: ", cutoffBest
  ! PRINT *, "splitFound: ", splitFound
  ! PRINT *, "casesOut: ", casesOut
  ! PRINT *, "nCuts: ", nCuts
  ! PRINT *, "lft: ", lft

  ! PRINT *, "============================END OF tfindSplit============================"

  RETURN

END SUBROUTINE tfindSplit
! =================================================================================

! Calculate the Kaplan Meier estimator
! ns integer, the number of time points
! nj real(:), at risk
! oj real(:), events
! z real(:), estimator
SUBROUTINE kaplan(ns, nj, oj, z)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ns
  REAL(dp), DIMENSION(1:ns), INTENT(IN) :: nj
  REAL(dp), DIMENSION(1:ns), INTENT(IN) :: oj
  REAL(dp), DIMENSION(1:ns), INTENT(OUT) :: z

  INTEGER :: i

  REAL(dp), DIMENSION(1:ns) :: num

  ! PRINT *, "******************** kaplan ********************"
  num = nj - oj

  z(1) = num(1) / nj(1)
  ! PRINT * , "z(1)", z(1)

  IF (ns .LT. 2) RETURN

  DO i = 2, ns
    IF (nj(i) > 1d-8) THEN
      z(i) = (num(i)/nj(i)) * z(i-1)
    ELSE
      z(i) = z(i-1)
    END IF
  END DO

  DO i = 1, ns
    IF (z(i) < 0.0 .OR. z(i) > 1.0) THEN
        PRINT *, "CHECK Kaplan Meier estimate!!"
        PRINT *, "z(i) wrong at time point ", i, ": ", z(i)
        PRINT *, "num(i): ", num(i)
        PRINT *, "nj(i): ", nj(i)
        PRINT *, "oj(i): ", oj(i)
        PRINT *, "delta", delta
        PRINT *, "delta_m", delta_m
        STOP  ! Terminate the entire program if any element        
    END IF
END DO

  ! PRINT *, "Kaplan Meier estimator for ns timepoints: ", z

END SUBROUTINE kaplan
! =================================================================================

SUBROUTINE nelsonAalenRecurrent(ns, nj, oj, h)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ns
  ! Number of distinct event times
  REAL(dp), DIMENSION(1:ns), INTENT(IN) :: nj
  ! nj(i): Number of subjects at risk just before time t_i, including those who have had previous events and are still under observation
  REAL(dp), DIMENSION(1:ns), INTENT(IN) :: oj
  ! oj(i): Number of events at time t_i, including recurrent events from the same subjects
  REAL(dp), DIMENSION(1:ns), INTENT(OUT) :: h
  ! h(i): Nelson-Aalen cumulative hazard estimator at time t_i

  INTEGER :: i

  ! Initialize the cumulative hazard estimator
  h(1) = oj(1) / nj(1)
  ! The initial cumulative hazard is the hazard at the first time point

  IF (ns .LT. 2) RETURN
  ! If there is only one time point, the initial value is the result

  ! Calculate the Nelson-Aalen estimator for recurrent events
  DO i = 2, ns
    IF (nj(i) > 1d-8) THEN
      h(i) = h(i-1) + oj(i) / nj(i)
      ! Add the hazard contribution at time t_i to the cumulative hazard from the previous time points
    ELSE
      h(i) = h(i-1)
      ! If nj(i) is very small, carry forward the previous cumulative hazard
    END IF
  END DO

END SUBROUTINE nelsonAalenRecurrent
! =================================================================================


! Truncated Mean test
! N1j: real(:), at risk in group 1
! N2j: real(:), at risk in group 2
! O1j: real(:), events in group 1
! O2j: real(:), events in group 2
! Z: real, truncated mean
SUBROUTINE meanSplit(N1j, N2j, O1j, O2j, Z)
  IMPLICIT NONE

  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: N1j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: N2j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: O1j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: O2j
  REAL(dp), INTENT(OUT) :: Z

  REAL(dp), DIMENSION(1:nt) :: E1, E2

  ! PRINT *, "******************** meanSplit ********************"
  CALL kaplan(nt, N1j, O1j, E1)
  ! PRINT *, "Kaplan Meier Estimate for Group 1: ", E1

  CALL kaplan(nt, N2j, O2j, E2)
  ! PRINT *, "Kaplan Meier Estimate for Group 2: ", E1

  Z = sum((E1 - E2) * dt)
  Z = Z*Z
  ! PRINT *, "Z statistic for truncated mean test: ", Z

END SUBROUTINE meanSplit
! =================================================================================

! Log rank test set up
! N1j: real(:), at risk in group 1
! N2j: real(:), at risk in group 2
! O1j: real(:), events in group 1
! O2j: real(:), events in group 2
! numJ: real(:), numerator
! denJ: real(:), denominator
SUBROUTINE logRankSetUp(N1j, N2j, O1j, O2j, numJ, denJ)

  IMPLICIT NONE

  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: N1j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: N2j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: O1j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: O2j
  REAL(dp), DIMENSION(1:nt), INTENT(OUT) :: numJ
  REAL(dp), DIMENSION(1:nt), INTENT(OUT) :: denJ

  INTEGER :: i

  REAL(dp) :: Nj, Oj

  ! PRINT *, "******************** LogRankSetup ********************"
  numJ = 0.d0
  denJ = 0.d0

  DO i = 1, nt
  	! if number at risk in group1 is equal to 0, then skip to next time point index (i)
    IF (N1j(i) .LT. 1d-8) CYCLE 
    ! if number at risk in group2 is equal to 0, then skip to next time point index (i) 
    IF (N2j(i) .LT. 1d-8) CYCLE
    ! time points for which both groups have individuals at risk

    ! number of individuals at risk for groups 1 and 2
    Nj = N1j(i) + N2j(i)

    ! number of events in group 1 and 2 at time point i
    Oj = O1j(i) + O2j(i)

    ! average (b/c combining two groups) event rate (total number of events in both groups divided by total number of people at risk in both groups at time i)
    numJ(i) = Oj / Nj

    ! variance?
    denJ(i) = numJ(i) * (Nj - Oj) / (Nj * Nj)

  END DO

END SUBROUTINE logrankSetUp
! =================================================================================

! Log rank test
! N1j: real(:), at risk in group 1
! N2j: real(:), at risk in group 2
! O1j: real(:), events in group 1
! O2j: real(:), events in group 2
! numJ: real(:), numerator
! denJ: real(:), denominator
! Z: real, test value
SUBROUTINE logRank(N1j, N2j, O1j, numJ, denJ, Z)

  IMPLICIT NONE

  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: N1j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: N2j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: O1j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: numJ
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: denJ
  REAL(dp), INTENT(OUT) :: Z

  INTEGER :: i

  REAL(dp) :: den, num

  ! PRINT *, "******************** logrank ********************"
  num = 0.d0
  den = 0.d0
  DO i = 1, nt
    IF (N1j(i) .LT. 1d-8) CYCLE
    IF (N2j(i) .LT. 1d-8) CYCLE
    ! time points for which both groups have individuals at risk

    num = num + O1j(i) - numJ(i) * N1j(i)
    den = den + denJ(i) * N1j(i) * N2j(i)

    ! Oj/Nj * [(Nj - Oj) / (Nj * Nj)] * N1j(i) * N2j(i)

  END DO

  IF (den .GT. 1d-8) THEN
    Z = num * num / den
  ELSE
    Z = 0.d0
  END IF
  ! PRINT *, "logrank test statistic Z = ", Z

END SUBROUTINE logrank
! =================================================================================

! Calculate the Aalen-Johansen CIF estimator for cause m for group k
! nsk integer, the number of time points for group k
! delta_mkj INTEGER(:), indicator vector for failure time for cause m in group k(=1 for cause m in group k; =0 otherwise)
! Nkj real(:), at risk in group k
! Okj real (:), events in group k
! CIF real(:), estimator
SUBROUTINE CIF_mk(nsk, Nkj, Okj, CIF, JumpCIF)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nsk
  ! INTEGER, DIMENSION(1:nsk), INTENT(IN) :: delta_mkj
  REAL(dp), DIMENSION(1:nsk), INTENT(IN) :: Nkj
  REAL(dp), DIMENSION(1:nsk), INTENT(IN) :: Okj
  REAL(dp), DIMENSION(1:nsk), INTENT(OUT) :: CIF
  REAL(dp), DIMENSION(1:nsk), INTENT(OUT) :: JumpCIF

  INTEGER :: i

  REAL(dp), DIMENSION(1:nsk) :: Zk

  ! Zk = KM estimator for group k at failure time from any cause
  call kaplan(nsk, Nkj, Okj, Zk) 

  ! if failure due to cause m in group k, then contribute Nkj(1), otherwise 0 for time point 1. ! this is Z=1 b/c time point 1 is always 0 here.
  CIF(1) = 1/Nkj(1) ! * Okj(1)
  JumpCIF(1) = CIF(1)

  ! if nsk = 1 then CIF = CIF(1)
  IF (nsk .LT. 2) RETURN 

  DO i = 2, nsk
      IF (Nkj(i) > 1d-8) THEN
        ! we loop through all the failure times due to any cause and only the cause m will contribute to CIF summation
        CIF(i) = CIF(i-1) + Zk(i-1) / Nkj(i) * Okj(i)
        JumpCIF(i) = Zk(i-1) / Nkj(i) * Okj(i)
      ELSE
        CIF(i) = CIF(i-1)
        JumpCIF(i) = 0
      END IF
  END DO

  DO i = 1, nsk
    IF (JumpCIF(i) < 0.0 .OR. JumpCIF(i) > 1.0) THEN
        PRINT *, "******************** CIF_mk (CIF estimator for cause m group k) ********************"
        PRINT *, "Zk = KM estimator for group k at failure time from any cause", Zk(i)
        PRINT *, "JUMPSIZE IS WRONG at TIME POINT: ", i
        PRINT *, "CIF(i): ", CIF(i)
        PRINT *, "JumpCIF(i): ", JumpCIF(i)
        PRINT *, "Nkj(i): ", Nkj(i)
        STOP  ! Terminate the entire program if any element satisfies the condition
    END IF
  END DO


END SUBROUTINE CIF_mk
!  CALL CIF_mk(nt1, delta_11j, N1j, O1j, CIF11, JumpCIF11) ! cause 1 group 1 
!  CALL CIF_mk(nt2, delta_12j, N2j, O2j, CIF12, JumpCIF12) ! cause 1 group 2 
!-----
!  CALL CIF_mk(nt1, delta_21j, N1j, O1j, CIF21, JumpCIF21) ! cause 2 group 1 
!  CALL CIF_mk(nt2, delta_22j, N2j, O2j, CIF22, JumpCIF22) ! cause 2 group 2 
! =================================================================================

! Gray's Test for cause m between groups 1 and 2
! ns1 integer, the number of time points for group 1
! ns2 integer, the number of time points for group 2
! delta_m1j integer(:), indicator vector for failure time for cause m in group 1(=1 for cause m in group 1; =0 otherwise)
! delta_m2j integer(:), indicator vector for failure time for cause m in group 2(=1 for cause m in group 2; =0 otherwise)
! N1j real(:), at risk in group 1
! O1j real (:), events in group 1
! N2j real(:), at risk in group 2
! O2j real (:), events in group 2
! Gray_CIF real(:), estimator
SUBROUTINE Gray_m(ns1, ns2, N1j, O1j, N2j, O2j, Gray_CIF)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ns1
  INTEGER, INTENT(IN) :: ns2
  ! INTEGER, DIMENSION(1:ns1), INTENT(IN) :: delta_m1j
  REAL(dp), DIMENSION(1:ns1), INTENT(IN) :: N1j
  REAL(dp), DIMENSION(1:ns1), INTENT(IN) :: O1j
  ! INTEGER, DIMENSION(1:ns2), INTENT(IN) :: delta_m2j
  REAL(dp), DIMENSION(1:ns2), INTENT(IN) :: N2j
  REAL(dp), DIMENSION(1:ns2), INTENT(IN) :: O2j
  REAL(dp), INTENT(OUT) :: Gray_CIF ! this is summation from 0 to tau we only need to create one number

  INTEGER :: i, j, time

  REAL(dp) :: sum1, sum2
  REAL(dp), DIMENSION(1:ns2) :: CIFm1, CIFm2, JumpCIFm1, JumpCIFm2

  ! PRINT *, "******************** Gray_m (test statistic for m cause) ********************"

  ! Check and print negative values in N1j
    DO i = 1, ns1
        IF (N1j(i) < 0.0) THEN
            PRINT *, "Negative values in N1j:"
            PRINT *, "N1j(", i, ") = ", N1j(i)
        END IF
    END DO

    ! Check and print negative values in O1j
    
    DO i = 1, ns1
        IF (O1j(i) < 0.0) THEN
            PRINT *, "Negative values in O1j:"
            PRINT *, "O1j(", i, ") = ", O1j(i)
        END IF
    END DO

    ! Check and print negative values in N2j
    
    DO j = 1, ns2
        IF (N2j(j) < 0.0) THEN
            PRINT *, "Negative values in N2j:"
            PRINT *, "N2j(", j, ") = ", N2j(j)
        END IF
    END DO

    ! Check and print negative values in O2j
    
    DO j = 1, ns2
        IF (O2j(j) < 0.0) THEN
            PRINT *, "Negative values in O2j:"
            PRINT *, "O2j(", j, ") = ", O2j(j)
        END IF
    END DO


  Gray_CIF = 0.d0
  CALL CIF_mk(ns1, N1j, O1j, CIFm1, JumpCIFm1) ! cause m group 1 
  CALL CIF_mk(ns2, N2j, O2j, CIFm2, JumpCIFm2) ! cause m group 2 
  
  DO time = 1, ns1
    IF (JumpCIFm1(time) <0.0 .OR. JumpCIFm1(time) > 1.0) THEN
      PRINT *, "******************** Gray_m (test statistic for m cause) ********************"
      PRINT *, "time point: ", time
      PRINT *, "JUMPCIFM1(TIME) IS WRONG."
      PRINT *, "N1j(time):", N1j(time)
      PRINT *, "O1j(time):", O1j(time)
      PRINT *, "CIFm1(time):", CIFm1(time)
    END IF
  END DO

  ! PRINT *, "CIFm1: m-cause CIF for group 1: ", CIFm1
  ! PRINT *, "CIFm2: m-cause CIF for group 2: ", CIFm2
  sum1 = 1*JumpCIFm1(1) !CIF(1) = 0 so 1-CIFm1(1)=1
  sum2 = 1*JumpCIFm2(1)

  DO i = 2, ns1
      IF (1-CIFm1(i) .LT. 1d-8) CYCLE

      sum1 = sum1 + 1/(1-CIFm1(i-1))*JumpCIFm1(i)*O1j(i) !only want those from group1 cause m
      
  END DO

  DO i = 2, ns2
      IF (1-CIFm2(i) .LT. 1d-8) CYCLE

      sum2 = sum2 + 1/(1-CIFm2(i-1))*JumpCIFm2(i)*O2j(i) !only want those from group2 cause m

  END DO

  Gray_CIF = sum1 - sum2
  ! PRINT *, "Gray CIF statistic: ", Gray_CIF

END SUBROUTINE Gray_m
! # comparing group 1 vs group 2 for node splitting
! CALL Gray_m(ns1, ns2, delta_11j, N1j, O1j, delta_12j, N2j, O2j, Gray_CIF1) ! for cause 1
! CALL Gray_m(ns1, ns2, delta_21j, N1j, O1j, delta_22j, N2j, O2j, Gray_CIF2) ! for cause 2
! =================================================================================

! estimate the survival/cif function and mean survival/cif time
!   nCases: integer, the number of elements in casesIn
!   casesIn: integer, the indices of the subset for which the value is
!     calculated
!   Func: real(:), the estimated survival/CIF function
!   mean: real, the estimated mean/CIF survival
SUBROUTINE calcValueSingle(nCases, casesIn, Func, mean)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nCases
  INTEGER, DIMENSION(1:nCases), INTENT(IN) :: casesIn
  REAL(dp), DIMENSION(1:nt), INTENT(OUT) :: Func
  REAL(dp), INTENT(OUT) :: mean

  INTEGER :: i, index1

  REAL(dp), DIMENSION(1:nt) :: Nj, Oj, Oj_m, Rb
  REAL(dp), DIMENSION(1:nt) :: jumpCIF

  ! PRINT *, "******************** calcValueSingle ********************"
  ! PRINT *, "estimate the survival/cif function and mean survival/cif time"
  Func = 0.d0
  mean = 0.d0

  ! PRINT *, "nCases: ", nCases
  ! number of at risk cases at each time point
  ! {nt}
  Rb = sum(pr(casesIn,:), DIM = 1)
  ! PRINT *, "Number of at risk cases at each time point: ", Rb

  Nj(1) = nCases !number at risk at first time point is everyone
  DO i = 2, nt
    ! at each time point after first tp, we take # at risk at previous tp and subtract the number of events that happened
    Nj(i) = Nj(i-1) - Rb(i-1)
  END DO
  ! PRINT *, "Number at Risk Cases at Each Time Point: Nj", Nj

  ! number of events at each time point
  ! {nt}
  DO i = 1, nt
    ! here: delta depends on if you are doing survival or CR
    Oj_m(i) = sum(pr(casesIn, i)*delta_m(casesIn))
    Oj(i) = sum(pr(casesIn, i)*delta(casesIn))
    ! PRINT *, "((((((((((((( CHECKING OJ(i) ))))))))))))): ", Oj(i)
    ! PRINT *, "((((((((((((( CHECKING OJ_m(i) ))))))))))))): ", Oj_m(i)

    IF (isPhase1) THEN
	    IF (Oj_m(i) /= Oj(i)) THEN
	      PRINT *, "ERROR: LINE 1146"
	      PRINT *, "time: ", i
	      PRINT *, "pr(casesIn, i): ", pr(casesIn, i)
	      PRINT *, "pr(casesIn, i)*delta_m(casesIn): ", pr(casesIn, i)*delta_m(casesIn)
	      PRINT *, "pr(casesIn, i)*delta(casesIn): ", pr(casesIn, i)*delta(casesIn)     
	      ! PRINT *, "casesIn", casesIn
	      PRINT *, "delta_m(casesIn): ", delta_m(casesIn)
	      PRINT *, "delta(casesIn): ", delta(casesIn)
	      PRINT *, "Oj_m(i) = ", Oj_m(i)
	      PRINT *, "Oj(i) = ", Oj(i)
	    END IF
	END IF

  END DO

  IF (isPhase1) THEN
    ! PRINT *, "STEP1 estimate KM SURVIVAL FUNCTION"
    ! Kaplan-Meier estimate survival function
    ! {nt}
    ! PRINT *, "Number of Events at Each Time Point: Oj", Oj
  	! PRINT *, "Number of Events at Each Time Point for Cause M: Oj_m", Oj_m
  	! PRINT *, "Number of time points: nt", nt
  	! PRINT *, "Number at Risk: Nj", Nj
    
    CALL kaplan(nt, Nj, Oj, Func)
    ! PRINT *, "Kaplan Meier Estimate Survival Function: ", Func
    ! mean survival time
    mean = sum(Func * dt)
    ! PRINT *, "mean survial time: ", mean
    IF (mean < 0.0) PRINT *, "ERROR: mean < 0"
  END IF

  IF (isPhase2CR) THEN
    ! PRINT *, "STEP2CR estimate CIF"
    ! AJ estimate of CIF at each time point {nt}
    CALL CIF_mk(nt, Nj, Oj_m, Func, jumpCIF)
    ! PRINT *, "end of STEP2CR estimate CIF"

    ! mean CIF time
    mean = sum(Func * dt)
    ! PRINT *, "mean cumulative incidence time: ", mean

    ! Check if the mean is greater than 5
    IF (mean < 0.0) THEN
      ! Print the arguments
      PRINT *, "AJ Estimate CIF: "
      PRINT *, Func
      PRINT *, "Jumpsize: "
      PRINT *, jumpCIF
      PRINT *, "dt: "
      PRINT *, dt
      PRINT *, "nt = ", nt
      PRINT *, "nCases (# of elements in casesIn): ", nCases
    END IF

  END IF

  RETURN
END SUBROUTINE calcValueSingle
! =================================================================================

! For factor covariates, calculated the mean survival/cif time and
! use as covariate values for splitting
! nCases integer, number of cases under consideration
! casesIn, integer(:), cases under consideration
! kv, integer, covariate under consideration
! array, real(:), covariate vector as mean survival times
SUBROUTINE getCovariate(nCases, casesIn, kv, array)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nCases
  INTEGER, DIMENSION(1:nCases), INTENT(IN) :: casesIn
  INTEGER, INTENT(IN) :: kv
  REAL(dp), DIMENSION(1:nCases), INTENT(OUT) :: array

  INTEGER :: i
  INTEGER, DIMENSION(1:nCases) :: ix
  INTEGER, DIMENSION(:), ALLOCATABLE :: ind

  REAL(dp) :: mean
  REAL(dp), DIMENSION(1:nt) :: Func

  LOGICAL, DIMENSION(1:nCases) :: inSubset

  ! PRINT *, "******************** getCovariate ********************"
  array = 0.d0

  ! convert covariate to integers to ensure correct equality tests
  ix = nint(x(casesIn,kv))

  DO i = 1, nCat(kv)

    ! identify cases that have this level
    inSubset = ix .EQ. i !indicator that covariate = i level

    ! ensure that there are individuals in the ith level
    IF (count(inSubset) .EQ. 0) CYCLE

    ! pack indices and estimate survival
    ind = pack(casesIn, inSubset)

    ! calculate the mean survival/cif time for each individual in this subset
    ! PRINT *, "getCovariate: CALL calcValueSingle"
    CALL calcValueSingle(size(ind), ind, Func, mean)

    WHERE (inSubset) array = mean

  END DO

END SUBROUTINE getCovariate
! =================================================================================

! Calculate survival/cif function, mean survival/cif time, and if appropriate
! survival/cif probability
! nCases, integer, number of cases under consideration
! casesIn, integer(:), indices of cases under consideration
! Func, real(:), estimated survival/cif function
! mean, real, estimated mean survival/cif
! Prob, real, estimated survival/cif probability or zero
SUBROUTINE tcalculateValue(nCases, casesIn, Func, mean, Prob)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nCases
  INTEGER, DIMENSION(1:nCases), INTENT(IN) :: casesIn
  REAL(dp), DIMENSION(1:nt), INTENT(OUT) :: Func
  REAL(dp), INTENT(OUT) :: mean
  REAL(dp), INTENT(OUT) :: Prob

  ! PRINT *, "******************** tcalculateValue ********************"
  Func = 0.d0
  mean = 0.d0
  Prob = 0.d0

  ! PRINT *, "tcalculateValue: CALL calcValueSingle"
  CALL calcValueSingle(nCases, casesIn, Func, mean)

  IF (.NOT. isSurvival) RETURN ! if not doing survival/cif probabilities then exit

  ! estimate survival/cif probability at SurvivalTime/CIFTime
  Prob = Func(sIndex) * (1.d0 - sFraction) + &
           & Func(sIndex+1) * sFraction
  IF (Prob .LT. 1d-8) Prob = 0.d0

  RETURN
END SUBROUTINE tCalculateValue
! =================================================================================

! grow each tree
! forestSurvFunc, real(:), survival function averaged over forest
! forestMean, real, mean survival averaged over forest
! forestSurvProb, real, survival probability averaged over forest
SUBROUTINE tsurvTree(forestSurvFunc, forestMean, forestSurvProb)
  IMPLICIT NONE

  REAL(dp), DIMENSION(1:nAll*nt), INTENT(OUT) :: forestSurvFunc
  REAL(dp), DIMENSION(1:nAll), INTENT(OUT) :: forestMean
  REAL(dp), DIMENSION(1:nAll), INTENT(OUT) :: forestSurvProb

  INTEGER :: i, iTree, j, k, lft, m, nc, ncur, splitFound, splitVar
  INTEGER, DIMENSION(1:sampleSize) :: indices, jdex, xrand
  INTEGER, DIMENSION(1:np) :: newstat, pindices
  INTEGER, DIMENSION(1:np, 1:nrNodes) :: cstat
  INTEGER, DIMENSION(1:nrNodes, 1:2) :: stm
  INTEGER, DIMENSION(1:nAll) :: allStatus
  INTEGER, DIMENSION(:), ALLOCATABLE :: ind, indOut, leftCases, rightCases, pind

  REAL(dp) :: srs
  REAL(dp), DIMENSION(1:nLevs) :: cutoffBest
  REAL(dp), DIMENSION(1:nrNodes) :: mean, Prob
  REAL(dp), DIMENSION(1:nAll) :: xm
  REAL(dp), DIMENSION(1:nt, 1:nrNodes) :: Func
  REAL(dp), DIMENSION(1:nrNodes, 1:(5+nLevs)) :: nMatrix
  REAL(dp), DIMENSION(1:nt, 1:nAll) :: tforestSurvFunc

  LOGICAL, DIMENSION(1:np) :: cand
  LOGICAL, DIMENSION(1:nAll) :: tst

  LOGICAL :: are_equal
  INTEGER :: index_delta

  are_equal = .TRUE.

  ! PRINT *, "******************** tsurvTree ********************"
  ! PRINT *, "nAll: ", nAll
  tforestSurvFunc = 0.d0
  forestMean = 0.d0
  forestSurvProb = 0.d0

  DO iTree = 1, nTree
    ! PRINT *, "iTree: ", iTree
    ! PRINT *, "################# Line 1319"
    ! WRITE(*, '(/,A,A,/)') 'iTree: ', iTree

    Func = 0.d0
    mean = 0.d0
    Prob = 0.d0
    nMatrix = 0.0
    allStatus = 1

    ! sample data and set local variables x, pr, and delta to the selected
    ! subset
    IF (replace .EQ. 1) THEN
      PRINT *, "sample with replacement"
      xrand = sampleWithReplace(nAll, sampleSize)
      n = sampleSize
      x = xAll(xrand,:)
      pr = prAll(xrand,:)
      delta = deltaAll(xrand)
      delta_m = deltaAll_m(xrand)
    ELSE IF (nAll .NE. sampleSize) THEN
      PRINT *, "sample without replacement"
      xrand = sampleWithoutReplace(nAll, sampleSize)
      n = sampleSize
      x = xAll(xrand,:)
      pr = prAll(xrand,:)
      delta = deltaAll(xrand)
      delta_m = deltaAll_m(xrand)
    ELSE
      ! PRINT *, "replace \neq 1 and nAll = sampleSize"
      n = sampleSize
      x = xAll
      pr = prAll
      delta = deltaAll
      delta_m = deltaAll_m
    END IF

    ! cutoff for identifying covariates to be explored
    srs = stratifiedSplit / REAL(np)

    ! indices for all cases
    indices = (/(i,i=1,n)/)
    jdex = indices

    ! indices for all covariates
    pindices = (/(i,i=1,np)/)

    !! initialize first node

    ! calculate survival function and mean survival of the first node
    ! IF (isPhase1) PRINT *, "calculate survival function and mean survival of the first node"
    ! IF (isPhase2CR) PRINT *, "calculate CIF function and mean CIF of the first node"
    ! PRINT *, "tsurvTree: CALL calcValueSingle"
    CALL calcValueSingle(n, indices, Func(:,1), mean(1))
    ! PRINT *, "--------calcValueSingle OUTPUT mean(1): ", mean(1)

    IF (isSurvival) THEN
      ! estimate survival/cif probability at SurvivalTime/CIFTime
      Prob(1) = Func(sIndex,1) * (1.d0 - sFraction) + &
                  & Func(sIndex+1,1) * sFraction
      IF (Prob(1) .LT. 1d-8) Prob(1) = 0.d0
    END IF

    ! determine if the node can split based on basic minimum requirements
    if (n .LE. nodeSize .OR. sum(delta_m(indices)) .LE. 1) THEN
      nMatrix(1,1) = -1
    ELSE
      nMatrix(1,1) = -2
    END IF

    cstat(:,1) = 0

    ! start and finish locations of indices in node
    stm(1,1) = 1
    stm(1,2) = n

    ! location of most recent storage location in matrices/vectors
    ! ncur is incremented when a node successfully splits indicating the
    ! location in the nodes list where base information for the each daughter
    ! is stored
    ncur = 1

    ! --===----= here

    ! PRINT *, "nrNodes: ", nrNodes
    DO k = 1, nrNodes

      ! if k is beyond current node count or
      ! current node count at limit, break from loop
      IF (k .GT. ncur .OR. ncur .GT. (nrNodes - 2)) EXIT
      ! PRINT *, "TEST7"

      ! if node is not to be split, cycle to next node
      IF (nint(nMatrix(k,1)) .EQ. -1) CYCLE
      ! PRINT *, "TEST8"

      ! indices for cases contained in node
      ind = jdex(stm(k,1):stm(k,2))

      ! if there are deficient variables, use only these variables
      cand = cstat(:,k) .LT. floor(srs * sum(cStat(:,k)))
      pind = pack(pindices, cand)
      IF (size(pind) .EQ. 0) pind = pindices
      ! PRINT *, "TEST9"

      ! split cases
      indOut = ind

      ! PRINT *, "################# Line 1441"
      CALL tfindSplit(size(ind), ind, size(pind), pind, splitVar, cutoffBest, &
                    & splitFound, indOut, nc, lft) 

      IF (splitFound .EQ. 0 ) THEN
        ! if no split available, set node k as terminal node
        nMatrix(k,1) = -1
        CYCLE
      END IF

      ! set node k to be interior (i.e. has split)
      nMatrix(k,1) = -3

      ! add split information to node
      nMatrix(k,4) = pindices(splitVar)
      nMatrix(k,5) = nc
      nMatrix(k,6:(6+nc-1)) = cutoffBest(1:nc)
      nMatrix(k,2) = ncur + 1
      nMatrix(k,3) = ncur + 2

      ! PRINT *, "increment the times the variable was used in a split."
      ! increment the times the variable was used in a split
      newstat = cStat(:,k)
      newstat(nint(nMatrix(k,4))) = newstat(nint(nMatrix(k,4))) + 1

      ! store new case order in jdex
      jdex(stm(k,1):stm(k,2)) = indOut

      !! left node
      ! PRINT *, "!!!!! left node !!!!!"

      ncur = ncur + 1

      ! index boundaries for cases in left node
      stm(ncur,1) = stm(k,1)
      stm(ncur,2) = stm(k,1) + lft - 1

      leftCases = jdex(stm(ncur,1):stm(ncur,2))
      ! PRINT *, "leftCases: ", leftCases

      ! get basic node information for left daughter
      ! IF (isPhase1) PRINT *, "get basic node information for left daughter."
      ! PRINT *, "tsurvTree: Line 1494: CALL calcValueSingle"
      CALL calcValueSingle(size(leftCases), leftCases, Func(:,ncur), &
                         & mean(ncur))

      ! estimate survival probability at SurvivalTime
      ! IF (isPhase1) THEN
      Prob(ncur) = Func(sIndex,ncur) * (1.d0 - sFraction) + &
             & Func(sIndex+1,ncur) * sFraction
      IF (Prob(ncur) .LT. 1d-8) Prob(1) = 0.d0
      ! END IF

      IF (size(leftCases) .LE. nodeSize .OR. sum(delta_m(leftCases)) .LE. 1) THEN
        ! if the number of cases in the node is at or below the minimum required
        ! or the number of uncensored event is only 1
        ! status is terminal
        nMatrix(ncur,1) = -1
      ELSE
        nMatrix(ncur,1) = -2
      END IF

      cstat(:,ncur) = newstat

      ! PRINT *, "################# Line 1524"
      !! right node
      ! PRINT *, "!!!!! right node !!!!!"

      ncur = ncur + 1

      ! index boundaries for cases in right node
      stm(ncur,1) = stm(k,1) + lft
      stm(ncur,2) = stm(k,2)

      ! retrieve left and right cases
      rightCases = jdex(stm(ncur,1):stm(ncur,2))

      ! calculate survival function and mean survival time for right node
      ! IF (isPhase1) PRINT *, "calculate survival function and mean survival time for right node."
      ! PRINT *, "tsurvTree: Line 1531: CALL calcValueSingle"
      CALL calcValueSingle(size(rightCases), rightCases, Func(:,ncur), &
                         & mean(ncur))

      ! PRINT *, "################# Line 1545"
      IF (isSurvival) THEN
        Prob(ncur) = Func(sIndex,ncur) * (1.d0 - sFraction) + &
                       & Func(sIndex+1,ncur) * sFraction
        IF (Prob(ncur) .LT. 1d-8) Prob(1) = 0.d0
      END IF

      IF (size(rightCases) .LE. nodeSize .OR. &
        & sum(delta_m(rightCases)) .LE. 1) THEN
        ! if the number of cases in the node is at or below the minimum required
        ! or the number of uncensored event (from cause M) is only 1
        ! status is terminal
        nMatrix(ncur,1) = -1
      ELSE
        nMatrix(ncur,1) = -2
      END IF

      cstat(:,ncur) = newstat

      ! retrieve the variable on which data is split
      m = nint(nMatrix(k,4))

      ! retrieve the covariate
      xm = xAll(:,m)

      IF (nCat(m) .LE. 1) THEN
        ! if a numeric variable, use the cutoff value to
        ! identify if individual i goes left
        ! PRINT *, "if a numeric variable, use the cutoff value to identify if individual i goes left"
        tst = xm .LE. nMatrix(k,6)
      ELSE
        ! if an unordered factor, use category to identify if individual
        ! i goes left
        ! PRINT *, "if an unordered factor, use category to identify if individual i goes left"
        tst = .FALSE.
        DO j = 1, nint(nMatrix(k,5))
          tst = tst .OR. nint(xm) .EQ. nint(nMatrix(k,j+5))
        END DO
      END IF

      DO j = 1, nAll
        IF (allStatus(j) .NE. k) CYCLE
        IF (tst(j)) THEN
          allStatus(j) = nint(nMatrix(k,2))
        ELSE
          allStatus(j) = nint(nMatrix(k,3))
        END IF
      END DO

    END DO

    ! ensure that all nodes that are not "interior" are "terminal"
    WHERE (nMatrix(:,1) .EQ. -2) nMatrix(:,1) = -1

    trees(iTree)%Func = Func(:,1:ncur)
    trees(iTree)%mean = mean(1:ncur)
    trees(iTree)%Prob = Prob(1:ncur)
    trees(iTree)%matrix = nMatrix(1:ncur,:)
    trees(iTree)%nNode = ncur

    tforestSurvFunc = tforestSurvFunc + Func(:,allStatus)
    forestMean = forestMean + mean(allStatus)
    forestSurvProb = forestSurvProb + Prob(allStatus)

  END DO

  ! -- to here
  !PRINT *, "################# Line 1617"

  ! This change is to eliminate a strange lto warning from R
  ! forestSurvFunc = reshape(tforestSurvFunc, (/nt*nAll/)) / nTree
  j = 0
  DO i = 1, SIZE(tforestSurvFunc,2)
    forestSurvFunc(j+1:j+SIZE(tforestSurvFunc,1)) = tforestSurvFunc(:,i)
    j = j + SIZE(tforestSurvFunc,1)
  END DO

  forestSurvFunc = forestSurvFunc / nTree
  forestMean = forestMean / nTree
  forestSurvProb = forestSurvProb / nTree

  ! PRINT *, "forestSurvFunc: ", forestSurvFunc
  ! PRINT *, "forestMean: ", forestMean
  ! PRINT *, "forestSurvProb: ", forestSurvProb

  ! PRINT *, "END OF SUBROUTINE TSURVTREE"

END SUBROUTINE tSurvTree

! =================================================================================

! use tree structure to predict for newData
SUBROUTINE predict(iTree, nr, nc, newData)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iTree
  INTEGER, INTENT(IN) :: nr
  INTEGER, INTENT(IN) :: nc
  REAL(dp), DIMENSION(:,:) :: newData

  INTEGER :: i, j, m
  INTEGER, DIMENSION(1:nr) :: stat

  REAL(dp), DIMENSION(1:nr) :: xm
  REAL(dp), DIMENSION(1:nr,1:nc) :: nData

  LOGICAL, DIMENSION(1:nr) :: sti, tst

  nData = reshape(newData, (/nr,nc/))

  ! column 1 is 0/1 indicating interior/terminal
  ! column 2 is the index of left
  ! column 3 is the index of right
  ! column 4 is the split variable
  ! column 5 is the number of cutoffs
  ! column 6:nCol are the cutoff values

  ! begin with every case in node 1
  stat = 1

  PRINT *, "******************** predict ********************"
  DO i = 1, trees(iTree)%nNode

    ! WRITE(*, '(/,A,A,/)') 'trees(iTree)%nNode number: ', i


    ! if terminal cycle
    IF (nint(trees(iTree)%matrix(i,1)) .EQ. -1) CYCLE

    ! identify individuals in this node
    sti = stat .EQ. i

    ! retrieve the variable on which data is split
    m = nint(trees(iTree)%matrix(i,4))

    ! retrieve the covariate
    xm = nData(:,m)

    IF (nCat(m) .LE. 1) THEN
      ! if a numeric variable, use the cutoff value to
      ! identify if individual i goes left
      tst = xm .LE. trees(iTree)%matrix(i,6)
    ELSE
      ! if an unordered factor, use category to identify if individual
      ! i goes left
      tst = .FALSE.
      DO j = 1, nint(trees(iTree)%matrix(i,5))
        tst = tst .OR. nint(xm) .EQ. nint(trees(iTree)%matrix(i,j+5))
      END DO
    END IF

    WHERE (tst .AND. sti) stat = nint(trees(iTree)%matrix(i,2))
    WHERE ( (.NOT. tst) .AND. sti) stat = nint(trees(iTree)%matrix(i,3))

  END DO

  forest%Func = forest%Func + trees(iTree)%Func(:,stat)
  forest%mean = forest%mean + trees(iTree)%mean(stat)
  ! PRINT *, "hi"
  ! PRINT *, "forest%Prob: ", forest%Prob
  ! PRINT *, "stat:", stat
  ! PRINT *, "trees(iTree)%Prob(stat): ", trees(iTree)%Prob(stat)
  forest%Prob = forest%Prob + trees(iTree)%Prob(stat)
  ! PRINT *, "end: forest%Prob: ", forest%Prob

  RETURN

END SUBROUTINE
! =================================================================================

END MODULE INNERS
! =================================================================================! =================================================================================

! n, integer, number of cases in data
! np, integer, number of covariates in data
! xt, integer(:), covariates
! nCat, integer(:), number of levels in each covariate (0 = continuous,
!   1 = ordered factors, 2+ = unordered factor)
! nt, integer, number of time points
! nNodes, integer, number of nodes
! tFunc, real(:), survival functions for each node
! mean, real(:), mean survival of each node
! Prob, real(:), survival probability of each node
! nCols, integer, number of columns in node matrix
! tnodes, real(:), node information
! predFunc, real(:), predicted survival functions
! predMean, real(:), predicted mean survival
! predProb, real(:), predicted survival probabilities
SUBROUTINE predictSurvTree(n, np, xt, nCat, nt, nNodes, tFunc, mean, &
  & Prob, nCols, tnodes, predFunc, predMean, predProb)
  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = selected_real_kind(15,307)

  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(IN) :: np
  REAL(dp), DIMENSION(1:n*np), INTENT(IN) :: xt
  INTEGER, DIMENSION(1:np), INTENT(IN) :: nCat
  INTEGER, INTENT(IN) :: nt
  INTEGER, INTENT(IN) :: nNodes
  REAL(dp), DIMENSION(1:nt*nNodes), INTENT(IN) :: tFunc
  REAL(dp), DIMENSION(1:nNodes), INTENT(IN) :: mean
  REAL(dp), DIMENSION(1:nNodes), INTENT(IN) :: Prob
  INTEGER, INTENT(IN) :: nCols
  REAL(dp), DIMENSION(1:nNodes*nCols), INTENT(IN) :: tnodes
  REAL(dp), DIMENSION(1:n*nt), INTENT(OUT) :: predFunc
  REAL(dp), DIMENSION(1:n), INTENT(OUT) :: predMean
  REAL(dp), DIMENSION(1:n), INTENT(OUT) :: predProb

  INTEGER :: i, j, m
  INTEGER, DIMENSION(1:n) :: stat

  REAL(dp), DIMENSION(1:n) :: xm
  REAL(dp), DIMENSION(1:nt, 1:n) :: tsurv
  REAL(dp), DIMENSION(1:n, 1:np) :: x
  REAL(dp), DIMENSION(1:nNodes, 1:nCols) :: nodes
  REAL(dp), DIMENSION(1:nt, 1:nNodes) :: Func

  LOGICAL, DIMENSION(1:n) :: sti, tst

  ! PRINT *, "******************** predictSurvTree ********************"
  ! PRINT *, "size of (mean survival of each node) is: ", SIZE(mean)
  ! PRINT *, mean 
  x = reshape(xt, (/n, np/))
  nodes = reshape(tnodes, (/nNodes, nCols/))
  Func = reshape(tFunc, (/nt, nNodes/))

  tsurv = 0.d0
  
  ! column 1 is indicator of interior/terminal
  ! column 2 is the index of left
  ! column 3 is the index of right
  ! column 4 is the split variable
  ! column 5 is the number of cutoffs
  ! column 6:nCol are the cutoff values

  ! begin with every case in node 1
  ! stat will contain the terminal node to which each case belongs
  stat = 1
  ! PRINT *, "TESTING0 stat: ", stat
  ! PRINT *, "nNodes:", nNodes

  DO i = 1, nNodes
    ! PRINT *, "NODE #: ", i 
    ! PRINT *, "TESTING1 stat: ", stat

    ! if terminal cycle
    IF (nint(nodes(i,1)) .EQ. -1) CYCLE

    ! identify individuals in this node
    sti = stat .EQ. i
    ! PRINT *, "indiv in this node:", sti
    ! PRINT *, "TESTING2 stat: ", stat

    ! retrieve the variable on which data is split
    m = nint(nodes(i,4))
    ! PRINT *, "retrieve the variable on which data is split: m = ", m

    ! retrieve the covariate
    xm = x(:,m)

    ! PRINT *, "nCat(m): ", nCat(m)

    IF (nCat(m) .LE. 1) THEN
      ! if a numeric variable, use the cutoff value to
      ! identify if individual i goes left
      ! PRINT *, "if a numeric variable, use the cutoff value to identify if inidv i goes left"
      tst = xm .LE. nodes(i,6)
    ELSE
      ! if an unordered factor, use category to identify if individual
      ! i goes left
      ! PRINT *, "if an unordered factor, use category to identify if indiv i goes left"
      tst = .FALSE.
      DO j = 1, nint(nodes(i,5))
        tst = tst .OR. nint(xm) .EQ. nint(nodes(i,j+5))
      END DO
    END IF
    !PRINT *, "TESTING3 stat: ", stat

    ! PRINT *,"nodes(i,2)", nodes(i,2)

    WHERE (tst .AND. sti) stat = nint(nodes(i,2))
    WHERE ( (.NOT. tst) .AND. sti) stat = nint(nodes(i,3))
    ! PRINT *, "TESTING4 stat: ", stat

  END DO

  ! retrieve appropriate values based on terminal node
  ! stat contains the terminal node each person is in
  tsurv = Func(:,stat)
  ! print *, "tsurv (first ", 5, " subjects)"
    ! Print the first N subjects of tsurv
    ! do i = 1, 5
       ! PRINT *, tsurv(:,i)
    ! end do

  ! This change is to eliminate a strange lto warning from R
!  predFunc = reshape(tsurv,(/n*nt/))
  j = 0
  DO i = 1, SIZE(tsurv,2)
    predFunc(j+1:j+SIZE(tsurv,1)) = tsurv(:,i)
    j = j + SIZE(tsurv,1)
  END DO

  ! print *, "predFunc (first ", 5, " subjects)"
    ! Print the first N subjects of predFunc
   !  do i = 1, 5
    !     PRINT *, predFunc(i)
    ! end do
  
  ! PRINT *, "TESTING100 stat: ", stat
  ! PRINT *, "size(stat):", size(stat)
  ! PRINT *, "mean(stat) which is terminal node mean surv time for each person:"
  ! PRINT *, mean(stat)
  ! PRINT *, "size(mean(stat)):", size(mean(stat))
  ! mean = mean survival of each node (size is node size)
  ! stat = terminal node index (size is sample size)
  ! mean(stat) takes the the mean value for each index for each person because stat is 1:n
  ! This means if person one is in 8th terminal node, then predMean has the 8th value of mean for its first index.
  predMean = mean(stat) ! selecting index of terminal node in mean array
  predProb = Prob(stat) ! selecting index of terminal node in Prob array
  ! PRINT *, "End of predictSurvTree"

END SUBROUTINE predictSurvTree
! =================================================================================

! set up basic information for the module
! t_nt, integer, the number of time points
! t_dt, real(:), the time differences between time points
! t_rs, real, the probability for a random split
! t_ERT, integer, the indicator of extremely randomized trees
! t_uniformSplit, integer, the indicator of method for determining cut-off
!   when using ERT
! t_nodeSize, integer, the minimum number of cases in each node
! t_minEvent, integer, the minimum number of events in each node
! t_rule, integer, 0 = mean, 1 = logrank
! t_sIndex, integer, the indices of time points that is closest to the
!   requested survival time
! t_sFraction, real, the fractional distance between time points the the
!   requested survival time
! t_stratifiedSplit, real, the coefficient for determining stratification
! t_replace, integer, indicator of sampling with replacement
SUBROUTINE setUpBasics(t_nt, t_dt, t_rs, t_ERT, t_uniformSplit, t_nodeSize, &
                     & t_minEvent, t_rule, t_sIndex, t_sFraction, &
                     & t_stratifiedSplit, t_replace)

  USE INNERS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: t_nt
  REAL(dp), DIMENSION(1:t_nt), INTENT(IN) :: t_dt
  REAL(dp), INTENT(IN) :: t_rs
  INTEGER, INTENT(IN) :: t_ERT
  INTEGER, INTENT(IN) :: t_uniformSplit
  INTEGER, INTENT(IN) :: t_nodeSize
  INTEGER, INTENT(IN) :: t_minEvent
  ! we changed t_rule from intent so that we can manipulate it so rule = 4 --> rule = 2 for CR (mean test)
  ! INTEGER, INTENT(IN) :: t_rule
  INTEGER :: t_rule
  INTEGER, INTENT(IN) :: t_sIndex
  REAL(dp), INTENT(IN) :: t_sFraction
  REAL(dp), INTENT(IN) :: t_stratifiedSplit
  INTEGER, INTENT(IN) :: t_replace

  ! PRINT *, "******************** setUpBasics ********************"
  nt = t_nt
  sIndex = t_sIndex
  sFraction = t_sFraction

  isSurvival = sIndex > 0

  isPhase1 = t_rule < 3
  isPhase2CR = t_rule > 2
  isPhase2 = isPhase2CR

  ! PRINT *, "isPhase1:", isPhase1
  ! PRINT *, "isPhase2CR:", isPhase2CR

  ! we set t_rule = 2 if in Phase 2 after assigning isPhase2CR so that we do truncated mean test
  ! below is why we changed t_rule from intent(in) to just integer
  IF (isPhase2CR) t_rule = 2

  IF (dtAllocated) DEALLOCATE(dt)

  ALLOCATE(dt(1:nt))

  dtAllocated = .TRUE.
  ! write(*, *) 'Check2 Fortran'

  dt = t_dt

  rs = t_rs
  ERT = t_ERT
  uniformSplit = t_uniformSplit
  nodeSize = t_nodeSize
  minEvent = t_minEvent
  rule = t_rule
  stratifiedSplit = t_stratifiedSplit
  replace = t_replace
  ! write(*,*) 'End of setUpBasics Fortran'

END SUBROUTINE setUpBasics
! =================================================================================

! set up basic information for the module that is step dependent
! t_n, integer, the number of cases under consideration
! t_np, integer, the number of covariates
! t_x, real(:), the covariates
! t_pr, real(:), the probability mass vector of survival function
! t_delta, integer(:), the indicator of censoring
! t_delta_m, integer(:), the indicator of censoring for cause m (for CR endpoint)
! t_mTry, integer, the maximum number of covariates to try for splitting
! t_nCat, integer(:), the number of categories in each covariate
! t_sampleSize, integer, the number of cases to sample for each tree
! t_ntree, integer, the number of trees in the forest
! t_nrNodes, integer, the maximum number of nodes
SUBROUTINE setUpInners(t_n, t_np, t_x, t_pr, t_delta, t_delta_m, t_mTry, t_nCat, &
                      & t_sampleSize, t_nTree, t_nrNodes)

  USE INNERS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: t_n
  INTEGER, INTENT(IN) :: t_np
  REAL(dp), DIMENSION(1:t_n*t_np), INTENT(IN) :: t_x
  REAL(dp), DIMENSION(1:nt*t_n), INTENT(IN) :: t_pr
  INTEGER, DIMENSION(1:t_n), INTENT(IN) :: t_delta
  INTEGER, DIMENSION(1:t_n), INTENT(IN) :: t_delta_m
  INTEGER, INTENT(IN) :: t_mTry
  INTEGER, DIMENSION(1:t_np), INTENT(IN) :: t_nCat
  INTEGER, INTENT(IN) :: t_sampleSize
  INTEGER, INTENT(IN) :: t_nTree
  INTEGER, INTENT(IN) :: t_nrNodes

  INTEGER :: i
  LOGICAL :: are_equal
  

  ! PRINT *, "******************** setUpInners ********************"
  nAll = t_n
  np = t_np

  IF (isAllocated) THEN
    DEALLOCATE(xAll, prAll, deltaAll, deltaAll_m, nCat, forest%Func, forest%mean,  &
             & forest%Prob, trees)
  END IF

  ALLOCATE(xAll(1:nAll,1:np))
  ALLOCATE(prAll(1:nAll, 1:nt))
  ALLOCATE(deltaAll(1:nAll))
  ALLOCATE(deltaAll_m(1:nAll))
  ALLOCATE(nCat(1:np))

  isAllocated = .TRUE.

  xAll = reshape(t_x, (/nAll,np/))
  prAll = reshape(t_pr, (/nAll,nt/))
  deltaAll = t_delta
  deltaAll_m = t_delta_m

  nCat = t_nCat
  nLevs = max(maxval(nCat),1)

  mTry = t_mTry
  sampleSize = t_sampleSize

  ALLOCATE(forest%Func(1:nt, 1:nAll))
  ALLOCATE(forest%mean(1:nAll))
  ALLOCATE(forest%Prob(1:nAll))
  forest%Func = 0.d0
  forest%mean = 0.d0
  forest%Prob = 0.d0

  nTree = t_nTree

  ALLOCATE(trees(1:nTree))

  nrNodes = t_nrNodes

END SUBROUTINE setUpInners
! =================================================================================

! access function for calculating forest
SUBROUTINE survTree(tFunc, mean, Prob)
  USE INNERS
  IMPLICIT NONE

  REAL(dp), DIMENSION(1:nrNodes*nt), INTENT(OUT) :: tFunc
  REAL(dp), DIMENSION(1:nrNodes), INTENT(OUT) :: mean
  REAL(dp), DIMENSION(1:nrNodes), INTENT(OUT) :: Prob

  ! PRINT *, "******************** survTree ********************"
  CALL tsurvTree(tFunc, mean, Prob)

END SUBROUTINE survTree
! =================================================================================

! access function for calculating forest
SUBROUTINE cifTree(tFunc, mean, Prob)
  USE INNERS
  IMPLICIT NONE

  REAL(dp), DIMENSION(1:nrNodes*nt), INTENT(OUT) :: tFunc
  REAL(dp), DIMENSION(1:nrNodes), INTENT(OUT) :: mean
  REAL(dp), DIMENSION(1:nrNodes), INTENT(OUT) :: Prob

  ! PRINT *, "******************** cifTree ********************"
  CALL tsurvTree(tFunc, mean, Prob)

END SUBROUTINE cifTree
! =================================================================================


! retrieve dimensions of node matrix for the iTree-th tree
SUBROUTINE treeDim(iTree, nr, nc)
  USE INNERS
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iTree
  INTEGER, INTENT(OUT) :: nr
  INTEGER, INTENT(OUT) :: nc

  ! PRINT *, "******************** treeDim ********************"
  nr = size(trees(iTree)%matrix,1)
  nc = size(trees(iTree)%matrix,2)

END SUBROUTINE treeDim
! =================================================================================

! retrieve the nodes, survival/cif function, mean survival/cif, and survival/cif probability
! for the iTree-th tree
SUBROUTINE getTree(iTree, nr, nc, nodes, Func, mean, Prob)
  USE INNERS
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iTree
  INTEGER, INTENT(IN) :: nr
  INTEGER, INTENT(IN) :: nc
  REAL(dp), DIMENSION(1:nr*nc), INTENT(OUT) :: nodes
  REAL(dp), DIMENSION(1:nt*nr), INTENT(OUT) :: Func
  REAL(dp), DIMENSION(1:nr), INTENT(OUT) :: mean
  REAL(dp), DIMENSION(1:nr), INTENT(OUT) :: Prob

  INTEGER :: i, j

  ! PRINT *, "******************** getTree ********************"
  ! This change is to eliminate a strange lto warning from R
!  nodes = reshape(trees(iTree)%matrix, (/nr*nc/))
  j = 0
  DO i = 1, SIZE(trees(iTree)%matrix,2)
    nodes(j+1:j+SIZE(trees(iTree)%matrix,1)) = trees(iTree)%matrix(:,i)
    j = j + SIZE(trees(iTree)%matrix,1)
  END DO

  ! This change is to eliminate a strange lto warning from R
!  Func = reshape(trees(iTree)%Func, (/nt*nr/))
  j = 0
  DO i = 1, SIZE(trees(iTree)%Func,2)
    Func(j+1:j+SIZE(trees(iTree)%Func,1)) = trees(iTree)%Func(:,i)
    j = j + SIZE(trees(iTree)%Func,1)
  END DO

  mean = trees(iTree)%mean
  Prob = trees(iTree)%Prob

END SUBROUTINE getTree

! =================================================================================


!FUNCTION RND(a, b) RESULT(dRes)
!  implicit none

!  REAL(dp), INTENT(IN) :: a
!  REAL(dp), INTENT(IN) :: b
!  REAL(dp) :: dRes

!  CALL random_number(dRes)
!END FUNCTION
! =================================================================================
