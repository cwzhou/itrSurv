!itrSurv.f90 in testpackage/itrSurv/src
MODULE INNERS
  IMPLICIT NONE

  PUBLIC

  INTEGER, PARAMETER :: dp = selected_real_kind(15,307)
  INTEGER, PARAMETER :: ng_set2 = 2
  INTEGER, PARAMETER :: nst_set1 = 1

  INTEGER, SAVE :: ERT ! 0/1 1 = use extremely randomized tree
  INTEGER, SAVE :: minEvent ! minimum number of events in a node !Phase1: death, Phase2CR: PC death, Phase2RE: recurrent events
  INTEGER, SAVE :: minEventSurv ! minimum number of Survival events in a node ! this is the same as minEvent for Phase1. Not used for Phase2CR (because CR is subset)
  INTEGER, SAVE :: mTry ! maximum number of covariates to try
  INTEGER, SAVE :: n  ! number of sampled cases
  INTEGER, SAVE :: n_surv  ! number of sampled SUBJECTS (needed for Phase2RE bc cases = records)
  INTEGER, SAVE :: nAll ! number of cases (records for Phase2RE)
  INTEGER, SAVE :: nAll_surv ! number of subjects (same as nAll for Phase1 and Phase2CR) (needed for Phase2RE because records != subjects)
  INTEGER, SAVE :: nLevs ! maximum number of levels
  INTEGER, SAVE :: nodeSize ! minimum number of cases in a node
  INTEGER, SAVE :: nodeSizeSurv ! minimum number of SUBJECTS in a node (needed for Phase2RE because records != subjects)
  INTEGER, SAVE :: np ! number of covariates
  INTEGER, SAVE :: nrNodes ! maximum number of nodes in a tree
  INTEGER, SAVE :: nt ! number of time points
  INTEGER, SAVE :: nt_death ! number of time points for survival
  INTEGER, SAVE :: nTree ! number of trees
  INTEGER, SAVE :: replace
  INTEGER, SAVE :: rule ! logrank:1; truncated mean: 2; gray: 3; gray2 (CSH): 4; Q_LR RE:5
  INTEGER, SAVE :: sampleSize
  INTEGER, SAVE :: sIndex ! for survival probability, the index of nearest time point
  INTEGER, SAVE :: uniformSplit ! 0/1 1 = random cutoff comes from values

  ! censoring indicator 1 = event from any cause (not censored; 0 = censored)
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: deltaAll
  ! censoring indicator 1 = event from cause m (not censored; 0 = censored) for CR
  ! recurrent event indicator for RE
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: deltaAll_m
  ! event indicator 1 = event from any cause (not censored; 0 = censored)
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: delta
  ! COMPETING RISKS: event indicator 1 = event from cause m (not censored for cause m for CR)
  ! RECURRENT EVENTS: event indicator 1 = RECURRENT EVENT (0 = not a recurrent event (could be death or censoring))
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: delta_m, id_RE2
  ! number of categories in each np covariate
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: nCat

  !INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: hehe
    !if (isPhase1) then
    !  hehe = delta
    !  print *, "hehe is delta"
    !end if

    !if (isPhase2RE .OR. isPhase2CR) then
    !  print *, "what is hehe?"
    !  PRINT *, hehe
    !end if

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
  ! at risk for RE vector for sampled cases in recurrent event setting
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: pr2
 ! at risk for DEATH vector for sampled cases in recurrent event setting
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: pr2surv
  ! probability mass vector of survival function for sampled cases (time points = death in RE setting)
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: prsurv
  ! probability mass vector of survival function for all cases
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: prAll
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: pr2All
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: prsurvAll  
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: pr2survAll
  ! id list for Phase2RE; used in conjunction with pr2 for death-at-risk in Phase2RE
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: id_RE

  REAL, DIMENSION(:), ALLOCATABLE, SAVE :: ord_responseAll
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: ord_causeindAll

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
  LOGICAL, SAVE :: isPhase2RE = .FALSE.
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

      subroutine find_unique(val, num_unique) ! (val, final, num_unique)
        implicit none
        integer, intent(in) :: val(:)               ! Input array
        !integer, allocatable, intent(out) :: final(:) ! Output array of unique elements
        integer, intent(out) :: num_unique           ! Number of unique elements
        integer :: i, min_val, max_val
        integer, allocatable :: unique(:)

        ! Initialize counters and limits
        i = 0
        min_val = minval(val) - 1   ! Start just below the minimum value
        max_val = maxval(val)       ! Maximum value in the array

        ! Allocate memory for unique values (at most the size of val)
        allocate(unique(size(val)))

        ! Loop to find unique elements
        do while (min_val < max_val)
            i = i + 1
            min_val = minval(val, mask=val > min_val)
            unique(i) = min_val
        end do

        ! Number of unique elements found
        num_unique = i

        !! Allocate final array to store only the unique elements
        !allocate(final(num_unique))
        !final = unique(1:num_unique)
        !! Deallocate temporary array
        !deallocate(unique)
    end subroutine find_unique

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
! for RE: this is each record and nCases is number of records
! we need a nCases_person and 
SUBROUTINE tfindSplit(nCases, casesIn, nv, varsIn, &
                    & splitVar, cutoffBest, splitFound, casesOut, nCuts, lft)
  use, intrinsic :: ieee_arithmetic
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

  LOGICAL, DIMENSION(nCases) :: mk

  INTEGER :: cnt, cnt_m, i, cc, ikv, j, jj, k, kv, l, nUncensored, nUncensored_m, ptr, rightNode, rightNode_m
  INTEGER :: nnn, iii, jjj, count_nonzero, ix
  INTEGER :: rUnifSet, rUnifSet_m, set, splitLeft, splitLeftFinal, tieCovariate
  INTEGER :: splitLeft_m, splitLeftFinal_m, splitLeft_doloop, splitLeftFinal_doloop
  INTEGER :: tieValue, variablesTried
  INTEGER, DIMENSION(1:nCases) :: cases, dSorted, dSorted_m, tcases !dSorted_RE
  REAL(dp), DIMENSION(1:nCases) :: sorted_cases
  INTEGER, DIMENSION(1:nv) :: variables
  INTEGER, DIMENSION(:), ALLOCATABLE :: ind, ind_m, indSingles, indSingles_m, leftCases, leftCases_m, rightCases, rightCases_m
  INTEGER, DIMENSION(:), ALLOCATABLE :: uncensoredIndices, uncensoredIndices_m

  INTEGER, DIMENSION(1:nAll) :: group_cr
  INTEGER, DIMENSION(:), ALLOCATABLE :: group_cr2, cases2

  REAL(dp) :: cutoff, maxValueSplit, maxValueXm, rUnif, rUnif_m, valuej, tester3a, tester3b
  REAL(dp), DIMENSION(1:nt_death) :: atRiskLeft, atRiskRight, D, denJ, numJ
  REAL(dp), DIMENSION(1:nt_death) :: eventsLeft, eventsRight
  REAL(dp), DIMENSION(1:nt) :: pd1, pd2, pd21, pd22
  REAL(dp), DIMENSION(1:nt_death) :: Rcum ! length is number of failure time points
  REAL(dp), DIMENSION(1:nt) :: atRiskLeft_m, atRiskRight_m, D_m
  REAL(dp), DIMENSION(1:nt) :: eventsLeft_m, eventsRight_m, Rcum_m, pd1_m, pd2_m
  REAL(dp), DIMENSION(1:nCases) :: xSorted
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: prl, prr, prl_m, prr_m, pr2l, pr2r 
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: pr2survl, pr2survr, prsurvl, prsurvr

  LOGICAL :: randomSplit
  LOGICAL, DIMENSION(:), ALLOCATABLE :: singles, singles_m

  LOGICAL :: are_equal
  INTEGER :: index_delta
  LOGICAL :: notEqual
  REAL(dp) :: rnd, random1

  ! for crstm  ! rho_set0, nst_set1, and ng_set2 defined in modules inners
  REAL(dp) :: rho_set0
  REAL(dp), DIMENSION(:), ALLOCATABLE :: y_set
  INTEGER, DIMENSION(:), ALLOCATABLE :: ig_setSplit, m_set
  INTEGER :: ist_set1(nAll)
  REAL(dp) :: ys_set(nAll)
  INTEGER :: ms_set(nAll)
  INTEGER :: igs_set(nAll)
  REAL(dp) :: s_set(ng_set2-1)
  REAL(dp) :: vs_set(ng_set2-1, ng_set2-1)
  REAL(dp) :: v_set(ng_set2*(ng_set2-1)/2)
  REAL(dp) :: st_set(ng_set2-1)
  REAL(dp) :: vt_set(ng_set2*(ng_set2-1)/2)
  REAL(dp), DIMENSION(ng_set2*(4+3*ng_set2)) :: wk_set
  INTEGER, DIMENSION(4*ng_set2) :: iwk_set

  ! below is old code that doesn't initialize
  ! REAL(dp) :: s_set(ng_set2-1), vs_set(ng_set2-1, ng_set2-1)
  ! REAL(dp) :: v_set(ng_set2*(ng_set2-1)/2), st_set(ng_set2-1), vt_set(ng_set2*(ng_set2-1)/2)
  ! REAL(dp), DIMENSION(ng_set2*(4+3*ng_set2)) :: wk_set
  ! integer, dimension(4*ng_set2) :: iwk_set

  EXTERNAL :: rnd
  are_equal = .TRUE.

  IF (isPhase2RE) THEN
    WRITE(*,'(/,A)') '============================ tfindSplit ============================'
    PRINT *, "******************** tfindSplit ********************"
    PRINT *, "casesIn"
    PRINT *, casesIn
    PRINT *, "nCases:", nCases
    PRINT *, "nodeSize:", nodeSize
    PRINT *, "nodeSizeSurv:", nodeSizeSurv
  END IF

  ! determine if this is to be a random split
  randomSplit = rnd(0.d0, 1.d0) <= rs
  !PRINT *, "rs =", rs

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
 
 ! nCases is 
  ! Terminal node criteria
  IF (isPhase2RE) THEN
      IF (nCases < 2*nodeSizeSurv) RETURN
  ELSE
      ! Ensure nodeSizeSurv equals nodeSize in non-Phase2RE cases
      IF (nodeSizeSurv .NE. nodeSize) THEN 
          PRINT *, "Error: nodeSizeSurv and nodeSize mismatch"
          STOP
      END IF
      IF (nCases < 2*nodeSize) RETURN
      ! nodeSizeSurv = nodeSize for isPhase1 and isPhase2CR
  END IF

  IF (isPhase2RE) THEN
    PRINT *, "tcases"
    PRINT *, tcases
    PRINT *, "================================"
    STOP
  END IF

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
      !PRINT *, "nCat(kv) > 1 so using mean survival/cif time for each factor level"
      CALL getCovariate(nCases, casesIn, kv, xSorted)
    ELSE
      !PRINT *, "nCat(kv) is either 0 = continuous or 1 = ordered factors"
      xSorted = x(casesIn,kv)
    END IF

    !PRINT *, "old xSorted"
    !PRINT *, SIZE(xSorted)
    !PRINT *, xSorted
    !PRINT *, "first cases"
    !PRINT *, cases
    !PRINT *, "caseseIn"
    !PRINT *, casesIn
    !PRINT *, "nCases", nCases

    cases = casesIn

    !PRINT *, "second cases"
    !PRINT *, cases

    ! sort the covariate and track the indices
    CALL qsort4(xSorted, cases, 1, nCases)
!    CALL hpsort_eps_epw(nCases, xSorted, cases, 1d-8)

    !PRINT *, "sorted xSorted"
    !!PRINT *, SIZE(xSorted)
    !PRINT *, xSorted

    ! sort event indicator data accordingly
    ! Phase 1: overall survival (delta = event indicator from any cause)
    dSorted = delta(cases) 
    ! Phase 2: CR (delta_m = indicator for event from cause m)
    !          RE (delta_m = indicator for recurrent event)
    dSorted_m = delta_m(cases) 
             
    !dSorted_RE = delta_RE(cases) ! Phase 2 for recurrent events (delta_RE = indicator for recurrent events)

IF (isPhase2RE) THEN
    PRINT *, "cases"
    PRINT *, cases
    PRINT *, "tcases"
    PRINT *, tcases
    PRINT *, "dSorted"
    PRINT *, dSorted
END IF

    ! ******************** splitBoundaries ********************
    ! identify minimum cases for left and right splits 
    ! based on minimum uncensored cases, minimum node size, 
    ! and assurance that all equal valued cases are included 
    ! in the minimum nodes
    !PRINT *, "******************** splitBoundaries ********************"

    rUnif = 0.d0
    rUnif_m = 0.d0
    rUnifSet = -1
    rUnifSet_m = -1
    !PRINT *, "test1"

    ! below is for dSorted == Phase 1
    ! Do below for both isPhase1 and isPhase2CR
    ! cases that are not-censored
    uncensoredIndices = pack(tcases, dSorted .EQ. 1)
    nUncensored = size(uncensoredIndices)
        IF (isPhase2RE) THEN
          PRINT *, "uncensoredIndices: ", uncensoredIndices
          PRINT *, "Number of Uncensored: ", nUncensored
          PRINT *, "minimum Recurrent Events: ", minEvent
          PRINT *, "minimum Deaths: ", minEventSurv
        END IF

    ! if too few cases to meet minimum number of uncensored cases, CYCLE
    IF (nUncensored .LT. (minEventSurv * 2)) CYCLE

    !! able to split and satisfy minimum number of events in each node

    ! cases to left include all indices up to and including minEvent case
    ! must have at least nodeSize cases
    splitLeft = max(uncensoredIndices(minEvent), nodeSize)
    !PRINT *, "test3a"

    ! move splitLeft up to include cases with equivalent values of x
    splitLeft = count(xSorted .LE. (xSorted(splitLeft) + 1e-8))
    !PRINT *, "test4a"

    ! cases to right
    ! include all indices down to and including nUncensored - minEvent + 1 case
    ! must have at least nodeSize cases
    rightNode = min(uncensoredIndices(nUncensored - minEvent + 1), &
                  & nCases - nodeSize + 1)
    !PRINT *, "test4.5a"

    ! move rightNode down to include cases with equivalent values of x
    ! splitLeftFinal is the last possible case for the left node
    splitLeftFinal = count(xSorted .LT. xSorted(rightNode))

    ! if the splitLeft index is above the splitLeftFinal index cycle,
    ! split is not possible
    IF (splitLeft .GT. splitLeftFinal) CYCLE
    !PRINT *, "test6a"

    ! below is for Phase 2 (Endpoint: CR or RE)
    ! cases that are not-censored
    uncensoredIndices_m = pack(tcases, dSorted_m .EQ. 1)
    nUncensored_m = size(uncensoredIndices_m)
    
    ! if too few cases to meet minimum number of uncensored cases, CYCLE
    IF (nUncensored_m .LT. (minEvent * 2)) CYCLE

    !! able to split and satisfy minimum number of events in each node

    ! cases to left include all indices up to and including minEvent case
    ! must have at least nodeSize cases
    splitLeft_m = max(uncensoredIndices_m(minEvent), nodeSize)
    !PRINT *, "test3b"

    ! move splitLeft up to include cases with equivalent values of x
    splitLeft_m = count(xSorted .LE. (xSorted(splitLeft_m) + 1e-8))
    !PRINT *, "test4b"

    ! cases to right
    ! include all indices down to and including nUncensored - minEvent + 1 case
    ! must have at least nodeSize cases
    rightNode_m = min(uncensoredIndices_m(nUncensored_m - minEvent + 1), &
                  & nCases - nodeSize + 1)
    !PRINT *, "test5b"

    ! move rightNode down to include cases with equivalent values of x
    ! splitLeftFinal is the last possible case for the left node
    splitLeftFinal_m = count(xSorted .LT. xSorted(rightNode_m))

    ! if the splitLeft index is above the splitLeftFinal index cycle,
    ! split is not possible
    IF (splitLeft_m .GT. splitLeftFinal_m) CYCLE
    
    !PRINT *, "test7"
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
        !PRINT *, "TESTER_UNIF1"

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

        ! IF (isPhase2CR) THEN ! comment out bc true for isPhase2RE too 
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
        !PRINT *, "TESTER_UNIF2"
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
         !PRINT *, "TESTER CHECKING: splitLeft_m: ", splitLeft_m, "splitLeft: ", splitLeft
         !PRINT *, "tester3a:", tester3a, "tester3b", tester3b
         ! IF (tester3a /= tester3b) THEN
      !PRINT *, "ERROR: TESTER: SPLITLEFT NEQ SPLITLEFT_M"
      !PRINT *, tester3b
      !PRINT *, tester3a
    ! END IF

      END IF

    END IF

    !PRINT *, "TESTER_UNIF-END"

    ! -1 is returned if cannot satisfy minimum requirements for nodes
    ! cycle to next covariate
    !PRINT *, "rUnifSet:", rUnifSet
    IF (rUnifSet .EQ. -1) CYCLE

    ! increment the number of covariates that have been explored
    variablesTried = variablesTried + 1

    !write(*,'(A)') '********************************* maxValue *********************************'
    !***************** maxValue ***************

    ! set initial values for outputs
    set = 0
    maxValueXm = 0.d0
    cutOff = 0.d0

    !**********************************************************************************************************
    !**********************************************************************************************************
    ! SPLITTING IS DONE (ABOVE). NOW WE LOOK AT CASES (SUBJECTS) IN LEFT AND RIGHT DAUGHTER NODES
    IF (isPhase2RE) THEN
    PRINT *, "SPLITTING IS DONE"
    END IF
    !**********************************************************************************************************
    !**********************************************************************************************************
    leftCases = cases(1:(splitLeft-1))
    rightCases = cases(splitLeft:nCases)
    leftCases_m = cases(1:(splitLeft_m-1))
    rightCases_m = cases(splitLeft_m:nCases)
    !PRINT *, "TEST00"
    !PRINT *, "prsurv"
    !PRINT *, prsurv
    !PRINT *, "pr"
    !PRINT *, pr
    !PRINT *, "isPhase2RE", isPhase2RE
    if (isPhase1) THEN
      prl = pr(leftCases,:) ! status change no matter what the status is (censoring, failure, cause specific failure, etc)
      prr = pr(rightCases,:)
      prl_m = prl
      prr_m = prr
    else if (isPhase2CR) THEN
      prl = pr(leftCases_m,:) ! status change no matter what the status is (censoring, failure, cause specific failure, etc)
      prr = pr(rightCases_m,:)
      prl_m = pr(leftCases_m,:) ! the exact same for CR
      prr_m = pr(rightCases_m,:) ! should be the same for CR
    else if (isPhase2RE) THEN
      prl = prsurv(leftCases_m,:) ! survival times
      prr = prsurv(rightCases_m,:) 
      prl_m = pr(leftCases_m,:) ! recurrent event times
      prr_m = pr(rightCases_m,:)
    END IF
    !PRINT *, "TEST0"

    IF (isPhase1) THEN
    ! below: this gives # events for survival. 
    eventsLeft = sum(prl * &
                   & spread(dSorted(1:(splitLeft-1)), 2, nt), DIM = 1)
    eventsRight = sum(prr * &
                    & spread(dSorted(splitLeft:nCases), 2, nt), DIM = 1)
    ELSE IF (isPhase2CR) THEN
    ! CR: this gives # events for survival (using regular pr and failure delta, but splitLeft_m)
    ! splitLeft_m was determined using endpoint delta to determine min # censored in each node
    ! use splitLeft_m
      eventsLeft = sum(prl * &
                   & spread(dSorted(1:(splitLeft_m-1)), 2, nt), DIM = 1)
      eventsRight = sum(prr * &
                    & spread(dSorted(splitLeft_m:nCases), 2, nt), DIM = 1)
    ELSE IF (isPhase2RE) THEN
      ! for RE we need survival events to be using nt_death not nt
      eventsLeft = sum(prl * &
                   & spread(dSorted(1:(splitLeft_m-1)), 2, nt_death), DIM = 1)
      eventsRight = sum(prr * &
                    & spread(dSorted(splitLeft_m:nCases), 2, nt_death), DIM = 1)
    END IF
    !PRINT *, "TEST1"

    ! use dSorted_m, which is endpoint delta (CR: cause m, RE: recurrent event)
    eventsLeft_m = sum(prl_m * &
                  & spread(dSorted_m(1:(splitLeft_m-1)), 2, nt), DIM = 1)
    eventsRight_m = sum(prr_m * &
                  & spread(dSorted_m(splitLeft_m:nCases), 2, nt), DIM = 1)
    !PRINT *, "Test22"

    IF (isPhase2CR) THEN
      !PRINT *, "isPhase2CR:", isPhase2CR
      ! Initialize the group vector with zeros
      group_cr = 0
      ! Assign 1 to the indices in leftCases
      group_cr(leftCases_m) = 1
      ! Assign 2 to the indices in rightCases
      group_cr(rightCases_m) = 2
    END IF
    !PRINT *, "TEST2"

    ! Print the result
    !PRINT *, "leftCases"
    !PRINT *, leftCases
    !PRINT *, "rightCases"
    !PRINT *, rightCases
    !PRINT *, "nCases", nCases
    !PRINT *, "cases"
    !PRINT *, cases

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

    ! looking at "at risk" set now (we want overall failure for both Step1 and Step2CR)
    IF (isPhase1 .OR. isPhase2CR) THEN
      pd1 = sum(prl, DIM = 1) ! for group 1
      pd2 = sum(prr, DIM = 1) ! for group 2
      !PRINT *, "group1: pd1 = sum(prl, DIM = 1): ", pd1
      !PRINT *, "group2: pd2 = sum(prr, DIM = 1) ", pd2
      ! at risk is same for Phase1 and Phase2CR because same subjects and same timepoint
      !PRINT *, "tt3"
      IF (isPhase1) THEN
        atRiskLeft(1) = splitLeft - 1
        atRiskRight(1) = nCases - splitLeft + 1
      ELSE IF (isPhase2CR) THEN
        atRiskLeft(1) = splitLeft_m - 1
        atRiskRight(1) = nCases - splitLeft_m + 1
      END IF
      !PRINT *, "tt4"

      !IF (isPhase1 .OR. isPhase2CR) THEN
      DO j = 2, nt
        atRiskLeft(j) = atRiskLeft(j-1) - pd1(j-1)
        atRiskRight(j) = atRiskRight(j-1) - pd2(j-1)
      END DO
      !END IF
      !PRINT *, "tt5"

    ELSE IF (isPhase2RE) THEN
      ! because this is witihin phase2re, there are multiple records per person
      ! we need to identify 1) at risk at survival times (using colsums pr2surv)
      ! 2) at risk at recurrent event times (using colsums of pr2)
      ! Earlier, we already defined:
      ! 3) death events at survival times (using prsurv and delta)
      ! 4) recurrent events at RE times (using pr and delta_m)

      ! only important for RE
      pr2l = pr2(leftCases_m,:)
      pr2r = pr2(rightCases_m,:)
      pr2survl = pr2surv(leftCases_m,:)
      pr2survr = pr2surv(rightCases_m,:)

      pd1_m = sum(pr2l, DIM = 1) ! RE at risk for group 1
      pd2_m = sum(pr2r, DIM = 1) ! RE at risk for group 2
      atRiskLeft_m = pd1_m
      atRiskRight_m = pd2_m
      pd1 = sum(pr2survl, DIM = 1) ! death at risk in RE setting for group 1
      pd2 = sum(pr2survr, DIM = 1) ! death at risk in RE setting for group 2
      atRiskLeft = pd1
      atRiskRight = pd2

      PRINT *, "leftCases_m"
      PRINT *, leftCases_m
      PRINT *, "rightCases_m"
      PRINT *, rightCases_m
      PRINT *, "nt_death", nt_death
      PRINT *, "nt", nt
      PRINT *, "--------------------------"
      PRINT *, "RE: at risk left"
      PRINT *, atRiskLeft_m
      PRINT *, "Death in RE: at risk left"
      PRINT *, atRiskLeft
      PRINT *, "--------------------------"
      PRINT *, "splitLeft"
      PRINT *, splitLeft
      PRINT *, "splitLeftFinal"
      PRINT *, splitLeftFinal

    END IF

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

    cnt = 1
    cnt_m = 1

    if (isPhase1) then
      splitLeft_doloop = splitLeft
      splitLeftFinal_doloop = splitLeftFinal
    ELSE
      splitLeft_doloop = splitLeft_m
      splitLeftFinal_doloop = splitLeftFinal_m
    end if

    DO j = splitLeft_doloop, splitLeftFinal_doloop
      !PRINT *, "j: ", j, "/ncnt: ", cnt, "/ncnt_m: ", cnt_m

      ! status change indicators for jth case
      pd1 = prr(cnt,:) ! to get events and at risk
      pd1_m = prr_m(cnt_m,:) 

      cnt = cnt + 1
      cnt_m = cnt_m + 1

      ! number of survival events for jth case
      D = pd1*delta(cases(j))
      !PRINT *, "* survival events D: ", D
      ! number of endpoint events for jth case
      ! number of priority cause events for jth case
      ! number of recurrent events for jth case
      D_m = pd1*delta_m(cases(j))
      !PRINT *, "* endpoint events D_m: ", D_m
      !PRINT *, "! add the jth case to the left node"
      ! add the jth case to the left node
      eventsLeft = eventsLeft + D
      eventsLeft_m = eventsLeft_m + D_m
      !PRINT *, '! remove the jth case from the right node'
      ! remove the jth case from the right node
      eventsRight = eventsRight - D
      eventsRight_m = eventsRight_m - D_m
      
      ! calculating cumulative at risk!
      IF (isPhase1 .OR. isPhase2CR) THEN
      ! same for phase1 and phase2cr
        Rcum(1) = 0.0
        Rcum(2) = pd1(1)
        !PRINT *, "Rcum for timepoints 1 and 2: ", Rcum
        ! from time point 3 to nt
        DO k = 3, nt
          IF (pd1(k-1) .GT. 1d-8) THEN
            Rcum(k) = Rcum(k-1) + pd1(k-1)
          ELSE
            Rcum(k) = Rcum(k-1)
          END IF
        END DO
        !PRINT *, "Rcum:", Rcum
        Rcum = 1.d0 - Rcum
        ! number at risk
        !PRINT *, "Rcum = 1.d0 - Rcum:", Rcum

        ! add the jth case to the left node
        atRiskLeft = atRiskLeft + Rcum
        ! remove the jth case from the right node
        atRiskRight = atRiskRight - Rcum
      END IF

      IF (isPhase2RE) THEN
        ! for survival time points
        Rcum(1) = 0.0
        Rcum(2) = atRiskLeft(1)
        PRINT *, "Rcum for timepoints 1 and 2: "
        PRINT *, Rcum
        ! from time point 3 to nt_death (number of unique observed failure times + 0 and tau)
        DO k = 3, nt_death
          IF (atRiskLeft(k-1) .GT. 1d-8) THEN
            Rcum(k) = Rcum(k-1) + atRiskLeft(k-1)
          ELSE
            Rcum(k) = Rcum(k-1)
          END IF
        END DO
        PRINT *, "Rcum:", Rcum
        !removing below 1-rcum bc we are not using HC pr, we are using pr2
        ! Rcum = 1.d0 - Rcum
        ! number at risk
        PRINT *, "printing 1.d0 - Rcum (but we removed it for RE):"
        PRINT *, 1.d0 - Rcum

        ! add the jth case to the left node
        atRiskLeft = atRiskLeft + Rcum
        ! remove the jth case from the right node
        atRiskRight = atRiskRight - Rcum

        Rcum_m(1) = 0.0
        Rcum_m(2) = atRiskLeft_m(1)
        PRINT *, "Rcum_m for timepoints 1 and 2: "
        PRINT *, Rcum_m
        ! from time point 3 to nt
        DO k = 3, nt
          IF (atRiskLeft_m(k-1) .GT. 1d-8) THEN
            Rcum_m(k) = Rcum_m(k-1) + atRiskLeft_m(k-1)
          ELSE
            Rcum_m(k) = Rcum_m(k-1)
          END IF
        END DO
        PRINT *, "Rcum_m:"
        PRINT *, Rcum_m
        !removing below 1-rcum bc we are not using HC pr, we are using pr2
        ! Rcum = 1.d0 - Rcum
        ! number at risk
        PRINT *, "printing 1.d0 - Rcum_m (but we removed it for RE):"
        PRINT *, 1.d0 - Rcum_m

        ! add the jth case to the left node
        atRiskLeft_m = atRiskLeft_m + Rcum_m
        ! remove the jth case from the right node
        atRiskRight_m = atRiskRight_m - Rcum_m

        PRINT *, "atRiskLeft"
        PRINT *, atRiskLeft
        PRINT *, "atRiskLeft_m"
        PRINT *, atRiskLeft_m

      END IF 

      ! if the case is not the last case with this covariate value, cycle
      IF (xSorted(j) .GE. (xSorted(j+1) - 1d-8)) CYCLE


      !write(*,'(/,/,A)'), '######### Line 836 IN FORTRAN'
      ! number of events
      !PRINT *, "=======TEST_R1: eventsRight:"
      !PRINT *, eventsRight
      !PRINT *, "=======TEST_R1: eventsRight for endpoint:"
      !PRINT *, eventsRight_m
      !PRINT *, "=======TEST_L1: eventsLeft:" 
      !PRINT *, eventsLeft
      !PRINT *, "=======TEST_L1: eventsLeft for endpoint:"
      !PRINT *, eventsLeft_m
      !PRINT *, "D: "
      !PRINT *, D
      !PRINT *, "D_m: "
      !PRINT *, D_m
      !PRINT *, "Rcum"
      !PRINT *, Rcum
      !PRINT *, "atRiskLeft_m:"
      !PRINT *, atRiskLeft_m
      !PRINT *, "atRiskRight:", atRiskRight
      !PRINT *, "=======TEST_R2: eventsRight:", eventsRight
      !PRINT *, "=======TEST_R2: eventsRight for cause m:", eventsRight_m
      !PRINT *, "=======TEST_L2: eventsLeft:", eventsLeft
      !PRINT *, "=======TEST_L2: eventsLeft for cause m:", eventsLeft_m
      !PRINT *, "numJ:", numJ
      !PRINT *, "denJ:", denJ
      ! write(*,'(/,/)')

      !PRINT *, "original valuej: ", valuej
      ! calculate test statistic
      IF (rule == 1) THEN
        ! PRINT *, "~~~~~ SPLITTING TEST: PHASE 1: logrank test ~~~~~"
        CALL logrank(atRiskLeft, atRiskRight, eventsLeft, numJ, &
                   & denJ, valuej)
        ! PRINT *, "logrank test statistic valuej = ", valuej
      ELSE IF (rule == 2) THEN
        ! PRINT *, "~~~~~ SPLITTING TEST: PHASE 1: truncated mean test ~~~~~"
        CALL meanSplit(atRiskLeft, atRiskRight, eventsLeft, eventsRight, valuej)
        ! PRINT *, "mean test statistic valuej = ", valuej
      ELSE IF (rule == 3) THEN
        ! PRINT *, "~~~~~ SPLITTING TEST: PHASE 2 (CR): gray's test ~~~~~"
        ! set up for crstm
        IF (isPhase2CR) THEN
          !!!!!!!PRINT *, "CR: Gray's Test Set-Up to split nodes"
          ! below are placeholder/initial values
          rho_set0 = 0;
          ist_set1 = 1
          sorted_cases = cases
          !!!!!!!PRINT *, "cases"
          !!!!!!!PRINT *, cases
          !!!PRINT *, "old sorted_cases"
          !!!PRINT *, sorted_cases
          CALL qsort4(sorted_cases, cases, 1, size(cases))
          !!!!!!!PRINT *, "testinggg qsort: sorted sorted_cases"
          !!!!!!!PRINT *, sorted_cases
          y_set = ord_responseAll(sorted_cases) ! import observed ordered times
          m_set = ord_causeindAll(sorted_cases) ! ordered failure status (0 = censored, 1 = pc, 2 = other cause)
          ig_setSplit = group_cr(sorted_cases); ! put in splitting indicator aka group membership
          !!!PRINT *, "----------------"
          ! Print values to verify
          !!!!!!!print *, 'y_set:', y_set
          !!!!!!!print *, 'ig_setSplit:', ig_setSplit
          !!!!!!!print *, 'm_set:', m_set
          !!!PRINT *, ">>>>>>>>>>> DOES THIS WORK???????? <<<<<<<<<<<"
          !!!!!!!PRINT *, "size(sorted_cases)", size(sorted_cases)
          CALL crstm(y_set, m_set, ig_setSplit, ist_set1, size(sorted_cases), rho_set0, nst_set1, ng_set2, &
        & s_set, vs_set, ys_set, ms_set, igs_set, v_set, st_set, vt_set, wk_set, iwk_set, valuej)
          !!!PRINT *, "end crstm"
          !!!PRINT *, ">>>>>>>>>>> DOES THIS WORK???????? <<<<<<<<<<<"
        ELSE
          PRINT *, "THIS IS WERID: itrSurv.f90 LINE 786: why are we in gray's test if isPhase2CR is false?"  
        END IF
        ! PRINT *, "Gray's test statistic valuej = ", valuej
      ELSE IF (rule == 4) THEN
        PRINT *, "~~~~~ SPLITTING TEST: PHASE 2 (CR): gray's test (kaplan way - CSH) ~~~~~"
        CALL Gray_m(nt, nt, atRiskLeft, eventsLeft_m, &
        & atRiskRight, eventsRight_m, valuej)
        PRINT *, "CSH test statistic valuej = ", valuej
      ELSE IF (rule == 5) THEN
        PRINT *, "~~~~~ SPLITTING TEST: PHASE 2 (RE): Q_LR (extension of gray's for RE) ~~~~~"
        ! atRisk is same for RE and survival within RE test statistic Q_LR
        !CALL MeanFreqFunc(nt, nt_death, atRiskLeft_m, eventsLeft_m, atRiskLeft_m, eventsLeft, &
        !& survRE_Left, dRhat_Left, mu_Left, dmu_Left)
        !CALL MeanFreqFunc(nt, nt_death, atRiskRight_m, eventsRight_m, atRiskRight_m, eventsRight, &
        !& survRE_Right, dRhat_Right, mu_Right, dmu_Right)
        valuej = 99
        PRINT *, "Q_LR test statistic valuej = ", valuej
      END IF

      ! Make sure there are no NaN test statistics
      IF (ieee_is_nan(valuej)) then
        PRINT *, "The test statistic is NaN. Stopping the execution."
        STOP
      ELSE 
        !PRINT *, "Test Statistic: ", valuej
        !print *,""
        !print *,""
      END IF

      IF ((set .EQ. 0) .OR. (valuej .GT. maxValueXm)) THEN
        ! if first value or value > current max, save
        IF (rUnifSet .EQ. 1) THEN
          !PRINT *, "rUnifSet: ", rUnifSet
          !PRINT *, "rUnifSet_m: ", rUnifSet_m
          !PRINT *, "rUnif: ", rUnif
          !PRINT *, "rUnif_m: ", rUnif_m
          cutoff = rUnif
        ELSE
          cutoff = (xSorted(j) + xSorted(j+1))/2.d0
        END IF
        !PRINT *, "cutoff: ", cutoff
        maxValueXm = valuej
        tieValue = 1
        set = 1

      ELSE IF (valuej > (maxValueXm - 1d-8)) THEN
        !PRINT *, "~~~~~ TEST 6 ~~~~~"
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
      !PRINT *, "SuCCESSful."

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
      !PRINT *, "~~~~~ TEST 7 ~~~~~"

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
  !PRINT *, "splitFound: ", splitFound
  !PRINT *, "splitVar: ", splitVar
  !PRINT *, "cutoffBest: ", cutoffBest
  !PRINT *, "splitFound: ", splitFound
  !PRINT *, "casesOut: ", casesOut
  !PRINT *, "nCuts: ", nCuts
  !PRINT *, "lft: ", lft

  !PRINT *, "============================END OF tfindSplit============================"

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

  !PRINT *, "******************** kaplan ********************"
  num = nj - oj

  z(1) = num(1) / nj(1)
  !PRINT * , "z(1)", z(1)

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

  !PRINT *, "Kaplan Meier estimator for ns timepoints: ", z

END SUBROUTINE kaplan
! =================================================================================

 SUBROUTINE sort_and_subset(og, cases, subset_vector)
    implicit none
    integer, dimension(:), intent(in) :: og
    integer, dimension(:), intent(in) :: cases
    integer, dimension(:), intent(out) :: subset_vector
    integer :: i, j, temp

    ! Sort the cases array in ascending order
    !do i = 1, size(cases)
    !    do j = i + 1, size(cases)
    !        if (cases(i) > cases(j)) then
    !            temp = cases(i)
    !            cases(i) = cases(j)
    !            cases(j) = temp
    !        end if
    !    end do
    !end do

    !PRINT *, "og(cases)"
    !PRINT *, og(cases)


    ! Initialize subset_vector with values based on the sorted cases
    !do i = 1, size(cases)
    !    subset_vector(i) = og(cases(i))
    !end do
    subset_vector = 0
END SUBROUTINE sort_and_subset

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

  !PRINT *, "******************** meanSplit ********************"
  CALL kaplan(nt, N1j, O1j, E1)
  !PRINT *, "Kaplan Meier Estimate for Group 1: ", E1

  CALL kaplan(nt, N2j, O2j, E2)
  !PRINT *, "Kaplan Meier Estimate for Group 2: ", E2

  Z = sum((E1 - E2) * dt)
  Z = Z*Z
  !PRINT *, "Z statistic for truncated mean test: ", Z

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

  !PRINT *, "******************** LogRankSetup ********************"
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

  !PRINT *, "******************** logrank ********************"
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
  !PRINT *, "logrank test statistic Z = ", Z

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
  CIF(1) = 0 ! update 9/1/24 changing to 0 from 1/Nkj(1) because we go from 0 to tau so 0 at time point 0 !old: 1/Nkj(1) ! * Okj(1)
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

  !PRINT *, "******************** Gray_m (test statistic for m cause) ********************"

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
  CALL CIF_mk(ns1, N1j, O1j, CIFm1, JumpCIFm1) ! cause m group k=1
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

  !PRINT *, "CIFm1: m-cause CIF for group 1: ", CIFm1
  !PRINT *, "CIFm2: m-cause CIF for group 2: ", CIFm2
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
  !PRINT *, "Gray CIF statistic: ", Gray_CIF

END SUBROUTINE Gray_m
! # comparing group 1 vs group 2 for node splitting
! CALL Gray_m(ns1, ns2, delta_11j, N1j, O1j, delta_12j, N2j, O2j, Gray_CIF1) ! for cause 1
! CALL Gray_m(ns1, ns2, delta_21j, N1j, O1j, delta_22j, N2j, O2j, Gray_CIF2) ! for cause 2
! =================================================================================

! =================================================================================
! USING A MODIFIED VERSION OF CMPRSK CRSTM AND CRST
! Copyright (C) 2000 Robert Gray
! Distributed under the terms of the GNU public license
SUBROUTINE crstm(y, m, ig, ist, no, rho, nst, ng, s, vs, ys, ms, igs, v, st, vt, wk, iwk,z)
use, intrinsic :: ieee_arithmetic
IMPLICIT NONE

! subroutine to calculate score and variance matrix for comparing
! cumulative incidence curves for a specific cause among groups.
!  test statistic given by s' inv(vs) s, dist approx chi-square(ng-1)
!
!  everything starting with i-n is integer, all others double precision
!
!  On input:
!    y is the failure times (sorted in increasing order)
!    m is coded 0 if censored, 1 if failed from the cause of interest
!       2 if failed from some other cause.
!    ig denotes group membership, must be coded 1,2,...,ng (ie
!       consecutive integers from 1 to ng), where ng is the number of groups
!    ist denotes strata membership, must be coded 1,2,...,nst, where nst
!       is the number of strata (code all 1's if only 1 strata)
!    no is the # of observations (length of y,m,ig,ist)
!    rho is the power used in the weight function in the test statistic
!    nst and ng are the # strata and # groups
!  Length of ys, ms, and igs must be at least as long as the size of the
!    largest strata
!  Length of s and st must be at least ng-1
!  Length of v and vt must be at least ng*(ng-1)/2
!  Length of vs must be at least (ng-1)^2
!  Length of wk must be at least ng*(4+3*ng)
!  Length of iwk must be at least 4*ng
!
!  On output:
!    s gives the scores for the first ng-1 groups, and
!    vs the estimated variance covariance matrix of these scores

! double precision: a,b,c,d,e,f,g,h, o,p,q,r,s,t,u,v,w,x,y,z
! integer: i,j,k,l,m,n

    ! Input variables
    integer, intent(in) :: ig(no), ist(no), m(no), no, nst, ng
    double precision, intent(in) :: y(no), rho
    double precision, intent(inout) :: ys(no)
    integer, intent(inout) :: igs(no), ms(no)

    ! Output variables
    double precision, intent(out) :: s(ng-1), vs(ng-1, ng-1)
    double precision, intent(out) :: v(ng*(ng-1)/2), st(ng-1), vt(ng*(ng-1)/2)
    double precision, intent(out) :: z

    ! Temporary arrays and variables
    double precision, dimension(ng*(4+3*ng)) :: wk
    integer, dimension(4*ng) :: iwk
    integer :: i, j, l, ks, ng1, ng2, n
    integer :: ng3, ng4

    !!!!!!!PRINT *, "-------- starting crstm -------------"
   !!!PRINT *, "CRSTM initial inputs:"
   !!!PRINT *, "y"
   !!!PRINT *, y
   !!!PRINT *, "m"
   !!!PRINT *, m
   !!!PRINT *, "ig"
   !!!PRINT *, ig
    
   !!!PRINT *, "CRSTM0"
    ng1 = ng - 1
    ng2 = ng * ng1 / 2
    l = 0

    ! Initialize s and v
    s = 0.0d0
    v = 0.0d0

    DO i = 1, ng1
      s(i) = 0
      DO j = 1, i
        l = l + 1
        v(l) = 0
      END DO
    END DO

!PRINT *, "CRSTM1"
    ! Loop over strata
    do ks = 1, nst !do 20
      n = 0
      do i = 1, no !do 21
        if (ist(i).ne.ks) EXIT ! go to 21
        n = n + 1
        ys(n) = y(i)
        ms(n) = m(i)
        igs(n) = ig(i)
      end do !21
      ng3 = 4 * ng + 1
      ng4 = ng * ng
      !PRINT *, "CRSTM2"
      !PRINT *, "----------------------------"
      !PRINT *, "ys"
      !PRINT *, ys
      !PRINT *, "ys(1)"
      !PRINT *, ys(1)
      !PRINT *, "----------------------------"
      !PRINT *, "ms"
      !PRINT *, ms
      !PRINT *, "ms(1)"
      !PRINT *, ms(1)
      !PRINT *, "----------------------------"
      !PRINT *, "igs"
      !PRINT *, igs
      !PRINT *, "igs(1)"
      !PRINT *, igs(1)
      !PRINT *, "----------------------------"
      !PRINT *, "n", n
      !PRINT *, "ng", ng
      !PRINT *, "rho", rho
      !PRINT *, "st", st
      !PRINT *, "vt", vt
      !!!!!PRINT *, "ng1", ng1
      !PRINT *, "ng2", ng2
      !PRINT *, "----------------------------"
      !PRINT *, "wk"
      !PRINT *, wk
      !PRINT *, "wk(1)", wk(1)
      !PRINT *, "----------------------------"
      !PRINT *, "iwk"
      !PRINT *, iwk
      !PRINT *, "iwk(1)", iwk(1)
      ! Call subroutine crst
      !PRINT *, "--- calling crst ----"
      call crst(ys(1), ms(1), igs(1), n, ng, rho, st, vt, ng1, ng2, &
                & wk(1), wk(ng+1), wk(2*ng+1), wk(3*ng+1), &
                & wk(ng3), wk(ng3+ng4), wk(ng3+2*ng4), &
                & wk(ng3+2*ng4+ng), iwk(1), iwk(ng+1))

!PRINT *, "***************************************************************"
!PRINT *, "* PRINT: CRST OUTPUT after CRST is over*"
!PRINT *, "ys:"
!PRINT *, ys(1:n)

!PRINT *, "ms:"
!PRINT *, ms(1:n)

!PRINT *, "igs:"
!PRINT *, igs(1:n)

!!!!!PRINT *, "st:"
!!!!!PRINT *, st

!PRINT *, "wk(1:ng):"
!PRINT *, wk(1:ng)

!PRINT *, "wk(ng+1:2*ng):"
!PRINT *, wk(ng+1:2*ng)

!PRINT *, "wk(2*ng+1:3*ng):"
!PRINT *, wk(2*ng+1:3*ng)

!PRINT *, "wk(3*ng+1:4*ng):"
!PRINT *, wk(3*ng+1:4*ng)

!PRINT *, "wk(ng3:ng3+ng4-1):"
!PRINT *, wk(ng3:ng3+ng4-1)

!PRINT *, "wk(ng3+ng4:ng3+2*ng4-1):"
!PRINT *, wk(ng3+ng4:ng3+2*ng4-1)

!PRINT *, "wk(ng3+2*ng4:ng3+2*ng4+ng-1):"
!PRINT *, wk(ng3+2*ng4:ng3+2*ng4+ng-1)

!PRINT *, "vt:"
!PRINT *, vt

!PRINT *, "iwk(1:ng):"
!PRINT *, iwk(1:ng)

!PRINT *, "iwk(ng+1:2*ng):"
!PRINT *, iwk(ng+1:2*ng)
!!!!!PRINT *, "* END OF CRST OUTPUT *"
!!!!!PRINT *, "***************************************************************"

!!!!!PRINT *, "s(1)", s(1)
    !!!!!PRINT *, "vs(1,1)", vs(1,1)
    !!!!!PRINT *, "ng1", ng1
    !!!!!PRINT *, "st(1)", st(1)

      l = 0
      ! Update s and v
      do i = 1, ng1 !do 23
        s(i) = s(i) + st(i)
        do j = 1, i !do 24
          l = l + 1
          v(l) = v(l) + vt(l)
        end do !24
        !!!!!PRINT *, "end do 24: s(1)", s(1)
      end do !23
      !!!!!PRINT *, "end do 23: s(1)", s(1)
    end do !20
    !!!!!PRINT *, "CRSTM3"
    !!!!!PRINT *, "s(1)", s(1)
    !!!!!PRINT *, "vs(1,1)", vs(1,1)
    !!!!!PRINT *, "v(1)", v(1)

    ! Populate vs matrix
    l = 0
    do i = 1, ng1 !do 31
        do j = 1, i !do 332
            l = l + 1
            vs(i, j) = v(l)
            vs(j, i) = vs(i, j)
        end do !332
        !PRINT *, "end do 332: vs(1,1)", vs(1,1)
    end do !31
    !!!!!PRINT *, "CRSTM4"
    !!!!!PRINT *, "s(1)", s(1)
    !!!!!PRINT *, "vs(1,1)", vs(1,1)
    z = s(1)*s(1)/vs(1,1) ! chi-sq test statistic for 2 groups (priority cause is first cause)
    !!!!!PRINT *, "z", z
    ! return
    !!!!!PRINT *, "========= end of CRSTM ========="

END SUBROUTINE crstm
! =================================================================================

 !           call crst(y = ys(1), m = ms(1), ig =igs(1), n = n, ng = ng, rho = rho, s = st, v = vt, ng1 = ng1, nv = ng2, &
 !                f1m = wk(1), f1 = wk(ng+1), skmm = wk(2*ng+1), skm = wk(3*ng+1), &
 !                 c= wk(ng3), a= wk(ng3+ng4), v3 = wk(ng3+2*ng4), &
 !                 v2 = wk(ng3+2*ng4+ng), rs=iwk(1), d=iwk(ng+1))
SUBROUTINE crst(y, m, ig, n, ng, rho, s, v, ng1, nv, f1m, f1, skmm, skm, c, a, v3, v2, rs, d)
    IMPLICIT NONE

    ! Declare variables with explicit types and dimensions
    integer, intent(in) :: n, ng, ng1, nv
    real(dp), intent(in) :: rho
    real(dp), dimension(n), intent(inout) :: y
    real(dp), dimension(ng1), intent(inout) :: s
    real(dp), dimension(ng), intent(inout) :: f1m, f1, skmm, skm, v3
    real(dp), dimension(ng, ng), intent(inout) :: c, a
    real(dp), dimension(nv), intent(inout) :: v
    real(dp), dimension(ng1, ng), intent(inout) :: v2
    integer, dimension(n), intent(inout) :: m, ig
    integer, dimension(ng), intent(inout) :: rs
    integer, dimension(0:2, ng), intent(inout) :: d

    integer :: i, j, k, l, ll, lu, nd1, nd2
    real(dp) :: fm, f, tr, tq, td, t1, t2, t3, t4, t5, t6, fb

    !PRINT *, "***************************************************************"
    !PRINT *, "********************* STARTING CRST *********************"
    !PRINT *, "***************************************************************"

!PRINT *, "y", y
!PRINT *, "m", m
!PRINT *, "ig", ig
!PRINT *, "n", n
!PRINT *, "ng", ng
!PRINT *, "rho", rho
!!!!!PRINT *, "within crstm, st = s"
!!!!!PRINT *, s
!PRINT *, "v", v
!PRINT *, "ng1", ng1
!PRINT *, "nv", nv
!PRINT *, "f1m", f1m
!PRINT *, "f1", f1
!!!PRINT *, "input skmm", skmm
!PRINT *, "skm", skm
!PRINT *, "c", c
!PRINT *, "a", a
!PRINT *, "v3", v3
!PRINT *, "v2", v2
!PRINT *, "rs", rs
!PRINT *, "d", d

    ! Initialize rs array
    rs(:) = 0

    ! Populate rs array based on ig
    do i = 1, n
        j = ig(i)
        rs(j) = rs(j) + 1
    end do

    ! Initialize arrays
    s(:) = 0.0
    f1m(:) = 0.0
    f1(:) = 0.0
    skmm(:) = 1.0
    !!!PRINT *, "skmm initialized"
    !!!PRINT *, skmm
    skm(:) = 1.0
    v3(:) = 0.0
    v(:) = 0.0
    v2(:, :) = 0.0
    c(:, :) = 0.0
    a(:, :) = 0.0

    ! Begin looping over unique times
    fm = 0.0
    f = 0.0
    ll = 1
    lu = ll
    
    !PRINT *, "STARTING CRST=====1"
50  do
        lu = lu + 1
        if (lu > n) exit
        if (y(lu) > y(ll)) exit
    end do
    !PRINT *, "STARTING CRST=====2"
    
    lu = lu - 1
    nd1 = 0
    nd2 = 0
    ! d will contain the # in each group censored, failed from
    ! cause 1, and failing from cause 2, at this time
    d(:, :) = 0
    !PRINT *, "STARTING CRST=====3"

    do i = ll, lu
        j = ig(i)
        k = m(i)
        d(k, j) = d(k, j) + 1
    end do

    !PRINT *, "STARTING CRST=====4"
    nd1 = sum(d(1, :))
    nd2 = sum(d(2, :))

    !if (nd1 == 0 .and. nd2 == 0) PRINT *, "LINE 1539: goto 90"
    if (nd1 == 0 .and. nd2 == 0) goto 90
    !PRINT *, "STARTING CRST=====5"

    tr = 0.0
    tq = 0.0

    do i = 1, ng
        if (rs(i) <= 0) cycle
        td = d(1, i) + d(2, i)
        !!!PRINT *, "group ", i, "rs: ", rs(i)
        !!!PRINT *, "td: ", td
        !!!PRINT *, "(rs(i) - td) / rs(i)", (rs(i) - td) / rs(i)
        !!!PRINT *, "group ", i, "skmm: ", skmm(i)
        ! skmm is left continuous, and skm right continuous, km est.
        skm(i) = skmm(i) * (rs(i) - td) / rs(i)
        !!!PRINT *, "For group ", i, "we have skm(i) to be:", skm(i)
        ! f1m is left continuous, and f1 right continuous, cuminc est.
        f1(i) = f1m(i) + (skmm(i) * d(1, i)) / rs(i)
        ! in notation of the Gray 1988 paper, tr is \sum_r\hat{h}_r, and tq is \sum_r R_r
        tr = tr + rs(i) / skmm(i)
        tq = tq + rs(i) * (1.0 - f1m(i)) / skmm(i)
    end do
    !PRINT *, "STARTING CRST=====6"

    f = fm + nd1 / tr
    fb = (1.0 - fm)**rho

    a(:, :) = 0.0

!PRINT *, "STARTING CRST=====7"
    do i = 1, ng
        if (rs(i) <= 0) cycle
        t1 = rs(i) / skmm(i)
        a(i, i) = fb * t1 * (1.0 - t1 / tr)
        if (a(i, i) /= 0.0) c(i, i) = c(i, i) + a(i, i) * nd1 / (tr * (1.0 - fm))

        do j = i + 1, ng
            if (rs(j) <= 0) cycle
            a(i, j) = -fb * t1 * rs(j) / (skmm(j) * tr)
            if (a(i, j) /= 0.0) c(i, j) = c(i, j) + a(i, j) * nd1 / (tr * (1.0 - fm))
        end do
    end do

    !PRINT *, "STARTING CRST=====8"
    ! Make a symmetric matrix
    do i = 2, ng
        k = i - 1
        a(i, 1:k) = a(1:k, i)
        c(i, 1:k) = c(1:k, i)
    end do

    !!!PRINT *, "Before Loop: s(1):", s(1)
    !!!PRINT *, "fb", fb
    !!!PRINT *, "d(1, 1)", d(1, 1)
    !!!PRINT *, "nd1", nd1 
    !!!PRINT *, "rs(1)", rs(1)
    !!!PRINT *, "f1m(1)", f1m(1) 
    !!!PRINT *, "skm"
    !!!PRINT *, skm
    !!!PRINT *, "tq", tq
    !!!PRINT *, "d(1, 1) - nd1 * rs(1) * (1.0 - f1m(1)) / (skmm(1) * tq)", d(1, 1) - nd1 * rs(1) * (1.0 - f1m(1)) / (skmm(1) * tq)
    !!!PRINT *, "---"
    do i = 1, ng1
        !!!print *, "DO LOOP: i is: ", i, " and ng1 is: ", ng1
        !!!print *, "rs(i):", rs(i)
        if (rs(i) <= 0) THEN 
          !!!PRINT *, "cycling"
          cycle
        else
          !!!PRINT *, "NOT CYCLING"
        end if
        s(i) = s(i) + fb * (d(1, i) - nd1 * rs(i) * (1.0 - f1m(i)) / (skmm(i) * tq))
        !!!if (s(1) == 0) THEN
        !!!   PRINT *, "HELPPPPPP: THERE IS AN issue: s(1) =", s(1)
        !!!ELSE
        !!!   PRINT *, "no issue :). s(1) = ", s(1)
        !!!END IF
    end do

    !!!!!PRINT *, "within crst, s = st pt2:", s(1)

    if (nd1 > 0) then
        do k = 1, ng
            if (rs(k) <= 0) cycle
            t4 = 1.0
            if (skm(k) > 0.0) t4 = 1.0 - (1.0 - f) / skm(k)
            t5 = 1.0
            if (nd1 > 1) t5 = 1.0 - (nd1 - 1) / (tr * skmm(k) - 1.0)
            t3 = t5 * skmm(k) * nd1 / (tr * rs(k))
            v3(k) = v3(k) + t4 * t4 * t3

            do i = 1, ng1
                t1 = a(i, k) - t4 * c(i, k)
                v2(i, k) = v2(i, k) + t1 * t4 * t3

                do j = 1, i
                    l = i * (i - 1) / 2 + j
                    t2 = a(j, k) - t4 * c(j, k)
                    v(l) = v(l) + t1 * t2 * t3
                end do
            end do
        end do
    end if

!PRINT *, "STARTING CRST=====4"
    if (nd2 > 0) then
        do k = 1, ng
            if (skm(k) <= 0.0 .or. d(2, k) <= 0) cycle
            t4 = (1.0 - f) / skm(k)
            t5 = 1.0
            if (d(2, k) > 1) t5 = 1.0 - (d(2, k) - 1.0) / (rs(k) - 1.0)
            t6 = rs(k)
            t3 = t5 * ((skmm(k)**2) * d(2, k)) / (t6**2)
            v3(k) = v3(k) + t4 * t4 * t3

            do i = 1, ng1
                t1 = t4 * c(i, k)
                v2(i, k) = v2(i, k) - t1 * t4 * t3

                do j = 1, i
                    l = i * (i - 1) / 2 + j
                    t2 = t4 * c(j, k)
                    v(l) = v(l) + t1 * t2 * t3
                end do
            end do
        end do
    end if

90  if (lu >= n) goto 30

    do i = ll, lu
        j = ig(i)
        rs(j) = rs(j) - 1
    end do

    fm = f
    f1m(:) = f1(:)
    skmm(:) = skm(:)

    ll = lu + 1
    lu = ll
    goto 50

30  l = 0
    do i = 1, ng1
        do j = 1, i
            l = l + 1
            do k = 1, ng
                v(l) = v(l) + c(i, k) * c(j, k) * v3(k)
                v(l) = v(l) + c(i, k) * v2(j, k)
                v(l) = v(l) + c(j, k) * v2(i, k)
            end do
        end do
    end do
!PRINT *, "***************************************************************"
!PRINT *, "========= end of CRST====="
!PRINT *, "...!PRINTing outputs..."

!PRINT *, "y:"
!PRINT *, y

!PRINT *, "m:"
!PRINT *, m

!PRINT *, "ig:"
!PRINT *, ig

!!!!!PRINT *, "s:"
!!!!!PRINT *, s

!PRINT *, "f1m:"
!PRINT *, f1m

!PRINT *, "f1:"
!PRINT *, f1

!PRINT *, "skmm:"
!PRINT *, skmm

!PRINT *, "skm:"
!PRINT *, skm

!PRINT *, "v:"
!PRINT *, v

!PRINT *, "c:"
!DO i = 1, ng
    !PRINT *, c(i, :)
!END DO

!PRINT *, "a:"
!DO i = 1, ng
!    PRINT *, a(i, :)
!END DO

!PRINT *, "v3:"
!PRINT *, v3

!PRINT *, "v2:"
!DO i = 1, ng1
!    PRINT *, v2(i, :)
!END DO

!PRINT *, "rs:"
!PRINT *, rs

!PRINT *, "d:"
!DO i = 0, 2
    !PRINT *, d(i, :)
!END DO

!!!PRINT *, "************************* END OF CRST **************************************"
    
    return
END SUBROUTINE crst
! =================================================================================

! =================================================================================
! Recurrent Event Phase 2: mu/dmu (Ghosh and Lin 2000)
! assuming no ties (aka unique observed RE times)
! Nj_endpoint: real(:), at risk for RE (using RE time points) (calculated using pr2 which is the same for RE and death events)
! Oj_endpoint: real(:), events for RE (calcualted using pr and delta_RE for RE)
! Nj_survival: real(:), at risk for death (using death time points) (calculated using pr2 which is the same for RE and death events)
! Oj_survival: real(:), events for death (calculated using pr and delta for death)
! survRE: real(:), survival for RE times
! dRhat: real(:), nelson aalen estimator for recurrence data
! mu: real(:), mean frequency function
! dmu: real(:), derivative of mff
! tp_survival is unique observed death times (plus 0 and tau)
! tp_endpoint is unique observed RE times (plus 0 and tau)

! nt_endpoint is the number of tp_endpoint (RE)
! nt_death is the number of tp_survival (death)




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

  !PRINT *, "******************** calcValueSingle ********************"
  !PRINT *, "estimate the survival/cif function and mean survival/cif time"
  Func = 0.d0
  mean = 0.d0

  ! We now calculate the number of at risk cases at each time point, using Rb
  ! {nt}
  Rb = sum(pr(casesIn,:), DIM = 1) ! DIM = 1 is columns in fortran, this is time
  !PRINT *, "sum(pr(casesIn,:), DIM = 1): ", Rb

  Nj(1) = nCases !number at risk at first time point is everyone
  DO i = 2, nt
    ! at each time point after first tp, we take # at risk at previous tp and subtract the number of events that happened
    Nj(i) = Nj(i-1) - Rb(i-1)
  END DO
  !PRINT *, "Number at Risk Cases at Each Time Point: Nj", Nj

  ! We now calculate the number of events at each time point
  ! {nt}
  DO i = 1, nt
    ! here: delta depends on if you are doing survival or CR
    Oj_m(i) = sum(pr(casesIn, i)*delta_m(casesIn))
    Oj(i) = sum(pr(casesIn, i)*delta(casesIn)) ! HOW MANY FAILURES AT EACH TIME POINT, * IS ELEMENT-WISE MULTIPLICATION.
    !PRINT *, "((((((((((((( CHECKING OJ(i) ))))))))))))): ", Oj(i)
    !PRINT *, "((((((((((((( CHECKING OJ_m(i) ))))))))))))): ", Oj_m(i)

    IF (isPhase1) THEN
      IF (Oj_m(i) /= Oj(i)) THEN
        PRINT *, "ERROR: LINE 1146"
        PRINT *, "time: ", i
        PRINT *, "pr(casesIn, i): ", pr(casesIn, i)
        PRINT *, "pr(casesIn, i)*delta_m(casesIn): ", pr(casesIn, i)*delta_m(casesIn)
        PRINT *, "pr(casesIn, i)*delta(casesIn): ", pr(casesIn, i)*delta(casesIn)
        PRINT *, "casesIn", casesIn
        PRINT *, "delta_m(casesIn): ", delta_m(casesIn)
        PRINT *, "delta(casesIn): ", delta(casesIn)
        PRINT *, "Oj_m(i) = ", Oj_m(i)
        PRINT *, "Oj(i) = ", Oj(i)
      END IF
  END IF

  END DO

  IF (isPhase1) THEN
    !PRINT *, "STEP1 estimate KM SURVIVAL FUNCTION"
    ! Kaplan-Meier estimate survival function
    ! {nt}
    !PRINT *, "Number of Events at Each Time Point: Oj", Oj
    !PRINT *, "Number of Events at Each Time Point for Cause M: Oj_m", Oj_m
    !PRINT *, "Number of time points: nt", nt
    !PRINT *, "Number at Risk: Nj", Nj

    CALL kaplan(nt, Nj, Oj, Func)
    !PRINT *, "Kaplan Meier Estimate Survival Function: ", Func
    ! mean survival time
    mean = sum(Func * dt)
    !PRINT *, "mean survial time: ", mean
    IF (mean < 0.0) PRINT *, "ERROR: mean < 0"
  END IF

  IF (isPhase2CR) THEN
    !PRINT *, "STEP2CR estimate CIF"
    ! AJ estimate of CIF at each time point {nt}
    CALL CIF_mk(nt, Nj, Oj_m, Func, jumpCIF)
    !PRINT *, "end of STEP2CR estimate CIF"

    ! mean CIF time
    mean = sum(Func * dt)
    !PRINT *, "mean cumulative incidence time: ", mean

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

  !PRINT *, "******************** getCovariate ********************"
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
    !PRINT *, "getCovariate: CALL calcValueSingle"
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

  !PRINT *, "******************** tcalculateValue ********************"
  Func = 0.d0
  mean = 0.d0
  Prob = 0.d0

  !PRINT *, "tcalculateValue: CALL calcValueSingle"
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

  INTEGER :: i, iTree, j, k, lft, m, nc, ncur, splitFound, splitVar, i_surv
  INTEGER, DIMENSION(1:sampleSize) :: indices, jdex, xrand
  !INTEGER, DIMENSION(1:sampleSize_surv) :: indices_surv
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

  integer :: n_surv    ! Declare the integer that stores the number of unique elements
  integer, allocatable :: val(:), final(:)  ! Declare allocatable arrays


  are_equal = .TRUE.

  !PRINT *, "******************** tsurvTree ********************"
  !PRINT *, "nAll: ", nAll
  !PRINT *, "sampleSize", sampleSize
  tforestSurvFunc = 0.d0
  forestMean = 0.d0
  forestSurvProb = 0.d0

  DO iTree = 1, nTree
    !PRINT *, "iTree: ", iTree
    !PRINT *, "################# Line 1319"
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
      n = sampleSize ! the number of cases to sample for each tree: people for Phase1,2CR; records (RE) for Phase2RE
      id_RE2 = id_RE(xrand)
      x = xAll(xrand,:)
      pr = prAll(xrand,:)
      pr2 = pr2All(xrand,:)
      prsurv = prsurvAll(xrand,:)
      pr2surv = pr2survAll(xrand,:)
      delta = deltaAll(xrand)
      delta_m = deltaAll_m(xrand)
    ELSE IF (nAll .NE. sampleSize) THEN
      PRINT *, "sample without replacement"
      xrand = sampleWithoutReplace(nAll, sampleSize)
      n = sampleSize
      id_RE2 = id_RE(xrand)
      x = xAll(xrand,:)
      pr = prAll(xrand,:)
      pr2 = pr2All(xrand,:)
      prsurv = prsurvAll(xrand,:)
      pr2surv = pr2survAll(xrand,:)
      delta = deltaAll(xrand)
      delta_m = deltaAll_m(xrand)
    ELSE
      !PRINT *, "replace \neq 1 and nAll = sampleSize"
      n = sampleSize
      id_RE2 = id_RE
      x = xAll
      pr = prAll
      pr2 = pr2All
      prsurv = prsurvAll
      pr2surv = pr2survAll
      delta = deltaAll
      delta_m = deltaAll_m
    END IF

    IF (isPhase2RE) THEN
      PRINT *, "id_RE2:"
      PRINT *, id_RE2
      !PRINT *, "sampleSize", sampleSize
      PRINT *, "Number of REs:", n
      ! Call the subroutine to find unique elements
      call find_unique(id_RE2, n_surv) ! (id_RE2, num_unique)
      ! Print the result
      !print *, 'Unique values: ', final
      print *, 'Number of subjects: ', n_surv
    END IF

    ! cutoff for identifying covariates to be explored
    srs = stratifiedSplit / REAL(np)

    ! indices for all cases
    indices = (/(i,i=1,n)/)
    PRINT *, "indices for all cases"
    PRINT *, indices
    !indices_surv = (/(i_surv, i_surv = 1, n_surv)/)
    jdex = indices

    ! indices for all covariates
    pindices = (/(i,i=1,np)/)
    PRINT *, "indices for all covariates: pindices"
    PRINT *, pindices
    PRINT *, "Phase1?", isPhase1
    PRINT *, "Phase2? (RE)", isPhase2RE
    PRINT *, "nt and nt_death"
    PRINT *, nt
    PRINT *, nt_death

    !! initialize first node

    ! calculate survival function and mean survival of the first node
    ! IF (isPhase1) PRINT *, "calculate survival function and mean survival of the first node"
    ! IF (isPhase2CR) PRINT *, "calculate CIF function and mean CIF of the first node"
    !PRINT *, "tsurvTree: CALL calcValueSingle"
    CALL calcValueSingle(n, indices, Func(:,1), mean(1))
    IF (isPhase2RE) THEN
      PRINT *, "--------end of calcValueSingle: OUTPUT mean(1): "
      PRINT *, mean(1)
    END IF

    IF (isSurvival) THEN
      ! estimate survival/cif probability at SurvivalTime/CIFTime
      Prob(1) = Func(sIndex,1) * (1.d0 - sFraction) + &
                  & Func(sIndex+1,1) * sFraction
      IF (Prob(1) .LT. 1d-8) Prob(1) = 0.d0
    END IF

    IF (isPhase2RE) THEN
      PRINT *, "delta_m(indices) and sum"
      PRINT *, delta_m(indices)
      PRINT *, "Number of current RE events here:", sum(delta_m(indices))
      PRINT *, "Number of current DEATH here:", sum(delta(indices))
      PRINT *, "n", n
      PRINT *, "nodeSize", nodeSize
      PRINT *, "n_surv", n_surv
      PRINT *, "nodeSizeSurv", nodeSizeSurv
    END IF
    ! For endpoint = CR: Phase 1: delta = delta_m
    ! For endpoint = RE: Phase 1: delta_m is just delta
    ! nMatrix: DIMENSION(1:nrNodes, 1:(5+nLevs))
    ! determine if the node can split based on basic minimum requirements
    if (n .LE. nodeSize .OR. sum(delta_m(indices)) .LE. 1) THEN
      if (isPhase2RE) THEN
        if (n_surv .LE. nodeSizeSurv .OR. sum(delta(indices)) .LE. 1) THEN
          nMatrix(1,1) = -1
        else 
          nMatrix(1,1) = -2
        end if
      ELSE 
        nMatrix(1,1) = -1
      end if
    ELSE
      nMatrix(1,1) = -2
    END IF

    !If (isPhase2RE) THEN
      !PRINT *, "nMatrix is"
      !PRINT *, nMatrix
      !STOP
    !END IF

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

    if (isPhase2RE) PRINT *, "nrNodes: ", nrNodes
    DO k = 1, nrNodes

      if (isPhase2RE) THEN 
        if (k .GE. 2) STOP
      END IF

      ! if k is beyond current node count or
      ! current node count at limit, break from loop
      IF (k .GT. ncur .OR. ncur .GT. (nrNodes - 2)) EXIT
      !PRINT *, "TEST7"

      ! if node is not to be split, cycle to next node
      IF (nint(nMatrix(k,1)) .EQ. -1) CYCLE
      !PRINT *, "TEST8"

      if (isPhase2RE) THEN
        PRINT *, "stm"
        PRINT *, stm
        PRINT *, "stm(k,1); k = 1."
        PRINT *, stm(k,1)
        PRINT *, "jdex(stm(k,1):stm(k,2))"
        PRINT *, jdex(stm(k,1):stm(k,2))
      END IF
      ! indices for cases contained in node
      ind = jdex(stm(k,1):stm(k,2))
      !PRINT *, "ind"
      !PRINT *, ind

      ! if there are deficient variables, use only these variables
      cand = cstat(:,k) .LT. floor(srs * sum(cStat(:,k)))
      pind = pack(pindices, cand)
      IF (size(pind) .EQ. 0) pind = pindices
      !PRINT *, "TEST9"

      ! split cases
      indOut = ind

      if (isPhase2RE) THEN
        PRINT *, size(ind)
        PRINT *, "ind"
        PRINT *, ind
        print *, "pind"
        PRINT *, pind
      END IF 

      !PRINT *, "################# Line 1441"
      CALL tfindSplit(size(ind), ind, size(pind), pind, splitVar, cutoffBest, &
                    & splitFound, indOut, nc, lft)
      !PRINT *, "TEST10"

      IF (splitFound .EQ. 0 ) THEN
        ! if no split available, set node k as terminal node
        nMatrix(k,1) = -1
        CYCLE
      END IF
      !PRINT *, "TEST11"

      ! set node k to be interior (i.e. has split)
      nMatrix(k,1) = -3
      !PRINT *, "TEST12"

      ! add split information to node
      nMatrix(k,4) = pindices(splitVar)
      nMatrix(k,5) = nc
      nMatrix(k,6:(6+nc-1)) = cutoffBest(1:nc)
      nMatrix(k,2) = ncur + 1
      nMatrix(k,3) = ncur + 2

      !PRINT *, "increment the times the variable was used in a split."
      ! increment the times the variable was used in a split
      newstat = cStat(:,k)
      newstat(nint(nMatrix(k,4))) = newstat(nint(nMatrix(k,4))) + 1

      ! store new case order in jdex
      jdex(stm(k,1):stm(k,2)) = indOut

      !! left node
      !PRINT *, "!!!!! left node !!!!!"

      ncur = ncur + 1

      ! index boundaries for cases in left node
      stm(ncur,1) = stm(k,1)
      stm(ncur,2) = stm(k,1) + lft - 1

      leftCases = jdex(stm(ncur,1):stm(ncur,2))
      !PRINT *, "leftCases: ", leftCases

      ! get basic node information for left daughter
      ! IF (isPhase1) PRINT *, "get basic node information for left daughter."
      !PRINT *, "tsurvTree: Line 1494: CALL calcValueSingle"
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

      !PRINT *, "################# Line 1524"
      !! right node
      !PRINT *, "!!!!! right node !!!!!"

      ncur = ncur + 1

      ! index boundaries for cases in right node
      stm(ncur,1) = stm(k,1) + lft
      stm(ncur,2) = stm(k,2)

      ! retrieve left and right cases
      rightCases = jdex(stm(ncur,1):stm(ncur,2))

      ! calculate survival function and mean survival time for right node
      ! IF (isPhase1) PRINT *, "calculate survival function and mean survival time for right node."
      !PRINT *, "tsurvTree: Line 1531: CALL calcValueSingle"
      CALL calcValueSingle(size(rightCases), rightCases, Func(:,ncur), &
                         & mean(ncur))

      !PRINT *, "################# Line 1545"
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
        !PRINT *, "if a numeric variable, use the cutoff value to identify if individual i goes left"
        tst = xm .LE. nMatrix(k,6)
      ELSE
        ! if an unordered factor, use category to identify if individual
        ! i goes left
        !PRINT *, "if an unordered factor, use category to identify if individual i goes left"
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

  !PRINT *, "forestSurvFunc: ", forestSurvFunc
  !PRINT *, "forestMean: ", forestMean
  !PRINT *, "forestSurvProb: ", forestSurvProb

  !PRINT *, "END OF SUBROUTINE TSURVTREE"

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

 !!!PRINT *, "******************** predict ********************"
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
  !PRINT *, "hi"
  !PRINT *, "forest%Prob: ", forest%Prob
  !PRINT *, "stat:", stat
  !PRINT *, "trees(iTree)%Prob(stat): ", trees(iTree)%Prob(stat)
  forest%Prob = forest%Prob + trees(iTree)%Prob(stat)
  !PRINT *, "end: forest%Prob: ", forest%Prob

  RETURN

END SUBROUTINE predict
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

  !PRINT *, "******************** predictSurvTree ********************"
  !PRINT *, "size of (mean survival of each node) is: ", SIZE(mean)
  !PRINT *, mean
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
  !PRINT *, "TESTING0 stat: ", stat
  !PRINT *, "nNodes:", nNodes

  DO i = 1, nNodes
    !PRINT *, "NODE #: ", i
    !PRINT *, "TESTING1 stat: ", stat

    ! if terminal cycle
    IF (nint(nodes(i,1)) .EQ. -1) CYCLE

    ! identify individuals in this node
    sti = stat .EQ. i
    !PRINT *, "indiv in this node:", sti
    !PRINT *, "TESTING2 stat: ", stat

    ! retrieve the variable on which data is split
    m = nint(nodes(i,4))
    !PRINT *, "retrieve the variable on which data is split: m = ", m

    ! retrieve the covariate
    xm = x(:,m)

    !PRINT *, "nCat(m): ", nCat(m)

    IF (nCat(m) .LE. 1) THEN
      ! if a numeric variable, use the cutoff value to
      ! identify if individual i goes left
      !PRINT *, "if a numeric variable, use the cutoff value to identify if inidv i goes left"
      tst = xm .LE. nodes(i,6)
    ELSE
      ! if an unordered factor, use category to identify if individual
      ! i goes left
      !PRINT *, "if an unordered factor, use category to identify if indiv i goes left"
      tst = .FALSE.
      DO j = 1, nint(nodes(i,5))
        tst = tst .OR. nint(xm) .EQ. nint(nodes(i,j+5))
      END DO
    END IF
    !PRINT *, "TESTING3 stat: ", stat

    !PRINT *,"nodes(i,2)", nodes(i,2)

    WHERE (tst .AND. sti) stat = nint(nodes(i,2))
    WHERE ( (.NOT. tst) .AND. sti) stat = nint(nodes(i,3))
    !PRINT *, "TESTING4 stat: ", stat

  END DO

  ! retrieve appropriate values based on terminal node
  ! stat contains the terminal node each person is in
  tsurv = Func(:,stat)
  ! print *, "tsurv (first ", 5, " subjects)"
    ! Print the first N subjects of tsurv
    ! do i = 1, 5
       !PRINT *, tsurv(:,i)
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

  !PRINT *, "TESTING100 stat: ", stat
  !PRINT *, "size(stat):", size(stat)
  !PRINT *, "mean(stat) which is terminal node mean surv time for each person:"
  !PRINT *, mean(stat)
  !PRINT *, "size(mean(stat)):", size(mean(stat))
  ! mean = mean survival of each node (size is node size)
  ! stat = terminal node index (size is sample size)
  ! mean(stat) takes the the mean value for each index for each person because stat is 1:n
  ! This means if person one is in 8th terminal node, then predMean has the 8th value of mean for its first index.
  predMean = mean(stat) ! selecting index of terminal node in mean array
  predProb = Prob(stat) ! selecting index of terminal node in Prob array
  !PRINT *, "End of predictSurvTree"

END SUBROUTINE predictSurvTree
! =================================================================================

! set up basic information for the module
! t_nt, integer, the number of observed unique event time points including 0 and tau; death for isPhase1/isPhase2CR and recurrent events for isPhase2RE
! t_nt_death, integer, the number of observed unique death time points including 0 and tau; needed for isPhase2RE
! t_dt, real(:), the time differences between time points
! t_rs, real, the probability for a random split
! t_ERT, integer, the indicator of extremely randomized trees
! t_uniformSplit, integer, the indicator of method for determining cut-off
!   when using ERT
! t_nodeSize, integer, the minimum number of cases in each node (Phase1/2CR: cases = subjects = people; Phase2RE: cases = recurrent events/records)
! t_nodeSizeSurv, integer, the minimum number of SUBJECTS in each node
! t_minEvent, integer, the minimum number of events in each node
! t_rule, integer, 1 = logrank, 2 = mean, 3 = CR gray (Gray 1988), 4 = CR csh, 5 = RE gray Q_LR (Ghosh & Lin 2000)
! t_sIndex, integer, the indices of time points that is closest to the
!   requested survival time
! t_sFraction, real, the fractional distance between time points the the
!   requested survival time
! t_stratifiedSplit, real, the coefficient for determining stratification
! t_replace, integer, indicator of sampling with replacement
SUBROUTINE setUpBasics(t_nt, t_nt_death, t_dt, t_rs, t_ERT, t_uniformSplit, &
                     & t_nodeSize, t_nodeSizeSurv, &
                     & t_minEvent, t_minEventSurv, &
                     & t_rule, t_sIndex, t_sFraction, &
                     & t_stratifiedSplit, t_replace)

  USE INNERS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: t_nt
  INTEGER, INTENT(IN) :: t_nt_death
  REAL(dp), DIMENSION(1:t_nt), INTENT(IN) :: t_dt
  REAL(dp), INTENT(IN) :: t_rs
  INTEGER, INTENT(IN) :: t_ERT
  INTEGER, INTENT(IN) :: t_uniformSplit
  INTEGER, INTENT(IN) :: t_nodeSize
  INTEGER, INTENT(IN) :: t_nodeSizeSurv
  INTEGER, INTENT(IN) :: t_minEvent
  INTEGER, INTENT(IN) :: t_minEventSurv
  INTEGER, INTENT(IN) :: t_rule
  INTEGER, INTENT(IN) :: t_sIndex
  REAL(dp), INTENT(IN) :: t_sFraction
  REAL(dp), INTENT(IN) :: t_stratifiedSplit
  INTEGER, INTENT(IN) :: t_replace

  !PRINT *, "******************** setUpBasics ********************"
  !PRINT *, "t_rule", t_rule
  !PRINT *, "isPhase1:", isPhase1
  !PRINT *, "isPhase2:", isPhase2
  !PRINT *, "isPhase2CR:", isPhase2CR
  !PRINT *, "isPhase2RE:", isPhase2RE
  isPhase1 = .FALSE.
  isPhase2 = .FALSE.
  isPhase2CR = .FALSE.
  isPhase2RE = .FALSE.
  !PRINT *, "t_rule", t_rule
  !PRINT *, "isPhase1:", isPhase1
  !PRINT *, "isPhase2:", isPhase2
  !PRINT *, "isPhase2CR:", isPhase2CR
  !PRINT *, "isPhase2RE:", isPhase2RE
  
  nt = t_nt
  nt_death = t_nt_death

  sIndex = t_sIndex
  sFraction = t_sFraction

  isSurvival = sIndex > 0

  ! 1: logrank; 2: mean; 3: cr gray; 4: cr csh; 5: re q_lr
  
  isPhase1 = t_rule < 3
  isPhase2 = t_rule > 2
  isPhase2CR = (t_rule == 3 .OR. t_rule == 4)
  isPhase2RE = t_rule == 5

  !PRINT *, "t_rule", t_rule
  !PRINT *, "isPhase1:", isPhase1
  !PRINT *, "isPhase2:", isPhase2
  !PRINT *, "isPhase2CR:", isPhase2CR
  !PRINT *, "isPhase2RE:", isPhase2RE

IF (isPhase1 .OR. isPhase2CR) THEN
   IF (nt /= nt_death) THEN
      PRINT *, "nt"
      PRINT *, nt
      PRINT *, "nt_death"
      PRINT *, nt_death
      PRINT *, "itrSurv.f90: Line 2888: nt should equal nt_death for Phase1 and Phase2CR"
      STOP
   END IF
END IF


  IF (dtAllocated) DEALLOCATE(dt)

  ALLOCATE(dt(1:nt))

  dtAllocated = .TRUE.
  ! write(*, *) 'Check2 Fortran'

  dt = t_dt

  rs = t_rs
  ERT = t_ERT
  uniformSplit = t_uniformSplit
  nodeSize = t_nodeSize
  nodeSizeSurv = t_nodeSizeSurv
  minEvent = t_minEvent
  minEventSurv = t_minEventSurv
  rule = t_rule
  stratifiedSplit = t_stratifiedSplit
  replace = t_replace
  ! write(*,*) 'End of setUpBasics Fortran'

END SUBROUTINE setUpBasics
! =================================================================================

! set up basic information for the module that is step dependent
! t_n, integer, the number of cases under consideration
! t_idvec, integer(:), the id label for Phase2 dataset (needed for MeanFreqFunc for RE endpoint to calculate at risk for death in RE setting)
! t_np, integer, the number of covariates
! t_x, real(:), the covariates
! t_pr, real(:), the probability mass vector of survival function
! t_pr2, real(:), the at-risk vector for subjects in recurrent event setting (ignore for Phase 1 OS, or Phase 2 CR)
! t_pr2surv, real(:) for RE: death at risk during phase2 (records x nt_death)
! t_prsurv, real(:) for RE: death pr during phase 2 (records x nt_death)
! t_ord_causeind, integer(:), status indicator for subjects (ordered by response) (vector of 0 for RE)
! t_ord_response, real(:), ordered response for subjects (vector of 0 for RE)
! t_delta, integer(:), the indicator of censoring
! t_delta_m, integer(:), the indicator of censoring for cause m (for CR endpoint)
! t_mTry, integer, the maximum number of covariates to try for splitting
! t_nCat, integer(:), the number of categories in each covariate
! t_sampleSize, integer, the number of cases to sample for each tree ! this is number of people for Phase1 and Phase2CR, but number of records for Phase2RE
! t_ntree, integer, the number of trees in the forest
! t_nrNodes, integer, the maximum number of nodes
SUBROUTINE setUpInners(t_n, t_n_surv, t_idvec, t_np, t_x, &
                      & t_pr, t_pr2, t_pr2surv, t_prsurv, t_ord_causeind, t_ord_response, &
                      & t_delta, t_delta_m, t_mTry, t_nCat, &
                      & t_sampleSize, t_nTree, t_nrNodes)

  USE INNERS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: t_n, t_n_surv
  INTEGER, DIMENSION(1:t_n), INTENT(IN) :: t_idvec
  INTEGER, INTENT(IN) :: t_np
  REAL(dp), DIMENSION(1:t_n*t_np), INTENT(IN) :: t_x
  REAL(dp), DIMENSION(1:nt*t_n), INTENT(IN) :: t_pr
  REAL(dp), DIMENSION(1:nt*t_n), INTENT(IN) :: t_pr2 ! for RE only
  REAL(dp), DIMENSION(1:nt_death*t_n), INTENT(IN) :: t_pr2surv ! for RE only
  REAL(dp), DIMENSION(1:nt_death*t_n), INTENT(IN) :: t_prsurv ! for RE only
  INTEGER, DIMENSION(1:t_n), INTENT(IN) :: t_ord_causeind ! for CR gray's test only
  REAL(dp), DIMENSION(1:t_n), INTENT(IN) :: t_ord_response ! for CR gray's test only
  INTEGER, DIMENSION(1:t_n), INTENT(IN) :: t_delta
  INTEGER, DIMENSION(1:t_n), INTENT(IN) :: t_delta_m ! endpoint delta
  INTEGER, INTENT(IN) :: t_mTry
  INTEGER, DIMENSION(1:t_np), INTENT(IN) :: t_nCat
  INTEGER, INTENT(IN) :: t_sampleSize
  INTEGER, INTENT(IN) :: t_nTree
  INTEGER, INTENT(IN) :: t_nrNodes

  INTEGER :: i
  LOGICAL :: are_equal


 !PRINT *, "******************** setUpInners ********************"
  nAll = t_n
  nAll_surv = t_n_surv
  !print *, "nAll is:", nAll
  ord_causeindAll = t_ord_causeind
  id_RE = t_idvec
  ord_responseAll = t_ord_response
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
  pr2All = reshape(t_pr2, (/nAll,nt/))
  prsurvAll = reshape(t_pr2, (/nAll,nt_death/))
  pr2survAll = reshape(t_pr2, (/nAll,nt_death/))
  deltaAll = t_delta
  deltaAll_m = t_delta_m

  nCat = t_nCat
  nLevs = max(maxval(nCat),1)

  mTry = t_mTry
  sampleSize = t_sampleSize !the number of cases to sample for each tree

  ALLOCATE(forest%Func(1:nt, 1:nAll))
  ALLOCATE(forest%mean(1:nAll))
  ALLOCATE(forest%Prob(1:nAll))
  forest%Func = 0.d0
  forest%mean = 0.d0
  forest%Prob = 0.d0

  nTree = t_nTree

  ALLOCATE(trees(1:nTree))

  nrNodes = t_nrNodes

  IF (isPhase2RE) THEN
    PRINT *, "******************** setUpInners ********************"
    PRINT *, "Number of cases under consideration: nAll:", nAll
    PRINT *, "Number of cases survival:", nAll_surv
    PRINT *, "id_RE"
    PRINT *, id_RE
    PRINT *, "Number of Covariates np: ", np
    PRINT *, "t_delta"
    PRINT *, t_delta
    PRINT *, "t_delta_m"
    PRINT *, t_delta_m
    PRINT *, "******************************************************"
  END IF

  !PRINT *, "End of SetupInners"

END SUBROUTINE setUpInners
! =================================================================================

! access function for calculating forest
SUBROUTINE survTree(tFunc, mean, Prob)
  USE INNERS
  IMPLICIT NONE

  REAL(dp), DIMENSION(1:nrNodes*nt), INTENT(OUT) :: tFunc
  REAL(dp), DIMENSION(1:nrNodes), INTENT(OUT) :: mean
  REAL(dp), DIMENSION(1:nrNodes), INTENT(OUT) :: Prob

  !PRINT *, "******************** survTree ********************"
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

  !PRINT *, "******************** cifTree ********************"
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

  !PRINT *, "******************** treeDim ********************"
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

  !PRINT *, "******************** getTree ********************"
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

! delete this later
SUBROUTINE temporary(s,vs,z)
USE INNERS
IMPLICIT NONE

double precision, intent(in) :: s, vs
double precision, intent(out) :: z

! CALL crstm()

z = s/sqrt(vs)
z = z*z

END SUBROUTINE temporary


! =================================================================================
! USING A MODIFIED VERSION OF CMPRSK CRSTM AND CRST
! Copyright (C) 2000 Robert Gray
! Distributed under the terms of the GNU public license
SUBROUTINE crstm(y, m, ig, ist, no, rho, nst, ng, s, vs, ys, ms, igs, v, st, vt, wk, iwk,z)
IMPLICIT NONE

! subroutine to calculate score and variance matrix for comparing
! cumulative incidence curves for a specific cause among groups.
!  test statistic given by s' inv(vs) s, dist approx chi-square(ng-1)
!
!  everything starting with i-n is integer, all others double precision
!
!  On input:
!    y is the failure times (sorted in increasing order)
!    m is coded 0 if censored, 1 if failed from the cause of interest
!       2 if failed from some other cause.
!    ig denotes group membership, must be coded 1,2,...,ng (ie
!       consecutive integers from 1 to ng), where ng is the number of groups
!    ist denotes strata membership, must be coded 1,2,...,nst, where nst
!       is the number of strata (code all 1's if only 1 strata)
!    no is the # of observations (length of y,m,ig,ist)
!    rho is the power used in the weight function in the test statistic
!    nst and ng are the # strata and # groups
!  Length of ys, ms, and igs must be at least as long as the size of the
!    largest strata
!  Length of s and st must be at least ng-1
!  Length of v and vt must be at least ng*(ng-1)/2
!  Length of vs must be at least (ng-1)^2
!  Length of wk must be at least ng*(4+3*ng)
!  Length of iwk must be at least 4*ng
!
!  On output:
!    s gives the scores for the first ng-1 groups, and
!    vs the estimated variance covariance matrix of these scores

! double precision: a,b,c,d,e,f,g,h, o,p,q,r,s,t,u,v,w,x,y,z
! integer: i,j,k,l,m,n

    ! Input variables
    integer, intent(in) :: ig(no), ist(no), m(no), no, nst, ng
    double precision, intent(in) :: y(no), rho
    double precision, intent(inout) :: ys(no)
    integer, intent(inout) :: igs(no), ms(no)

    ! Output variables
    double precision, intent(out) :: s(ng-1), vs(ng-1, ng-1)
    double precision, intent(out) :: v(ng*(ng-1)/2), st(ng-1), vt(ng*(ng-1)/2)
    double precision, intent(out) :: z

    ! Temporary arrays and variables
    double precision, dimension(ng*(4+3*ng)) :: wk
    integer, dimension(4*ng) :: iwk
    integer :: i, j, l, ks, ng1, ng2, n
    integer :: ng3, ng4

   !!!PRINT *, "CRSTM0"
    ng1 = ng - 1
    ng2 = ng * ng1 / 2
    l = 0

    ! Initialize s and v
    s = 0.0d0
    v = 0.0d0

    DO i = 1, ng1
      s(i) = 0
      DO j = 1, i
        l = l + 1
        v(l) = 0
      END DO
    END DO

!PRINT *, "CRSTM1"
    ! Loop over strata
    do ks = 1, nst !do 20
      n = 0
      do i = 1, no !do 21
        if (ist(i).ne.ks) EXIT ! go to 21
        n = n + 1
        ys(n) = y(i)
        ms(n) = m(i)
        igs(n) = ig(i)
      end do !21
      ng3 = 4 * ng + 1
      ng4 = ng * ng
     !!!PRINT *, "CRSTM2"
     !PRINT *, "----------------------------"
     !PRINT *, "ys"
     !PRINT *, ys
     !!!PRINT *, "ys(1)"
     !!!PRINT *, ys(1)
     !!!PRINT *, "----------------------------"
     !PRINT *, "ms"
     !PRINT *, ms
     !!!PRINT *, "ms(1)"
     !!!PRINT *, ms(1)
     !!!PRINT *, "----------------------------"
     !PRINT *, "igs"
     !PRINT *, igs
     !!!PRINT *, "igs(1)"
     !!!PRINT *, igs(1)
     !PRINT *, "----------------------------"
     !!!PRINT *, "n", n
     !!!PRINT *, "ng", ng
     !!!PRINT *, "rho", rho
     !!!PRINT *, "st", st
     !!!PRINT *, "vt", vt
     !PRINT *, "ng1", ng1
     !!!PRINT *, "ng2", ng2
     !!!PRINT *, "----------------------------"
     !!!PRINT *, "wk"
     !!!PRINT *, wk
     !!!PRINT *, "wk(1)", wk(1)
     !!!PRINT *, "----------------------------"
     !!!PRINT *, "iwk"
     !!!PRINT *, iwk
     !!!PRINT *, "iwk(1)", iwk(1)
      ! Call subroutine crst
      call crst(ys(1), ms(1), igs(1), n, ng, rho, st, vt, ng1, ng2, &
                & wk(1), wk(ng+1), wk(2*ng+1), wk(3*ng+1), &
                & wk(ng3), wk(ng3+ng4), wk(ng3+2*ng4), &
                & wk(ng3+2*ng4+ng), iwk(1), iwk(ng+1))

PRINT *, "***************************************************************"
PRINT *, "*!!!PRINT: CRST OUTPUT after CRST is over*"
PRINT *, "ys:"
PRINT *, ys(1:n)

PRINT *, "ms:"
PRINT *, ms(1:n)

PRINT *, "igs:"
PRINT *, igs(1:n)

PRINT *, "st:"
PRINT *, st

PRINT *, "wk(1:ng):"
PRINT *, wk(1:ng)

PRINT *, "wk(ng+1:2*ng):"
PRINT *, wk(ng+1:2*ng)

PRINT *, "wk(2*ng+1:3*ng):"
PRINT *, wk(2*ng+1:3*ng)

PRINT *, "wk(3*ng+1:4*ng):"
PRINT *, wk(3*ng+1:4*ng)

PRINT *, "wk(ng3:ng3+ng4-1):"
PRINT *, wk(ng3:ng3+ng4-1)

PRINT *, "wk(ng3+ng4:ng3+2*ng4-1):"
PRINT *, wk(ng3+ng4:ng3+2*ng4-1)

PRINT *, "wk(ng3+2*ng4:ng3+2*ng4+ng-1):"
PRINT *, wk(ng3+2*ng4:ng3+2*ng4+ng-1)

PRINT *, "vt:"
PRINT *, vt

PRINT *, "iwk(1:ng):"
PRINT *, iwk(1:ng)

PRINT *, "iwk(ng+1:2*ng):"
PRINT *, iwk(ng+1:2*ng)
PRINT *, "* END OF CRST OUTPUT *"
PRINT *, "***************************************************************"

      l = 0
      ! Update s and v
      do i = 1, ng1 !do 23
        s(i) = s(i) + st(i)
        do j = 1, i !do 24
          l = l + 1
          v(l) = v(l) + vt(l)
        end do !24
      end do !23
    end do !20
   !!!PRINT *, "CRSTM3"

    ! Populate vs matrix
    l = 0
    do i = 1, ng1 !do 31
        do j = 1, i !do 332
            l = l + 1
            vs(i, j) = v(l)
            vs(j, i) = vs(i, j)
        end do !332
    end do !31
   PRINT *, "CRSTM4"

    z = s(1)*s(1)/vs(1,1) ! chi-sq test statistic for 2 groups (priority cause is first cause)
   PRINT *, "z", z
    ! return
   PRINT *, "========= end of CRSTM ========="

END SUBROUTINE crstm
! =================================================================================

 !           call crst(y = ys(1), m = ms(1), ig =igs(1), n = n, ng = ng, rho = rho, s = st, v = vt, ng1 = ng1, nv = ng2, &
 !                f1m = wk(1), f1 = wk(ng+1), skmm = wk(2*ng+1), skm = wk(3*ng+1), &
 !                 c= wk(ng3), a= wk(ng3+ng4), v3 = wk(ng3+2*ng4), &
 !                 v2 = wk(ng3+2*ng4+ng), rs=iwk(1), d=iwk(ng+1))
SUBROUTINE crst(y, m, ig, n, ng, rho, s, v, ng1, nv, f1m, f1, skmm, skm, c, a, v3, v2, rs, d)
    IMPLICIT NONE

    ! Declare variables with explicit types and dimensions
    integer, intent(in) :: n, ng, ng1, nv
    double precision, intent(in) :: rho
    double precision, dimension(n), intent(inout) :: y
    double precision, dimension(ng1), intent(inout) :: s
    double precision, dimension(ng), intent(inout) :: f1m, f1, skmm, skm, v3
    double precision, dimension(ng, ng), intent(inout) :: c, a
    double precision, dimension(nv), intent(inout) :: v
    double precision, dimension(ng1, ng), intent(inout) :: v2
    integer, dimension(n), intent(inout) :: m, ig
    integer, dimension(ng), intent(inout) :: rs
    integer, dimension(0:2, ng), intent(inout) :: d

    integer :: i, j, k, l, ll, lu, nd1, nd2
    double precision :: fm, f, tr, tq, td, t1, t2, t3, t4, t5, t6, fb

   !!!PRINT *, "***************************************************************"
   !!!PRINT *, "********************* STARTING CRST *********************"
   !!!PRINT *, "***************************************************************"

PRINT *, "y", y
PRINT *, "m", m
PRINT *, "ig", ig
PRINT *, "n", n
PRINT *, "ng", ng
PRINT *, "rho", rho
PRINT *, "s", s
PRINT *, "v", v
PRINT *, "ng1", ng1
PRINT *, "nv", nv
PRINT *, "f1m", f1m
PRINT *, "f1", f1
PRINT *, "skmm", skmm
PRINT *, "skm", skm
PRINT *, "c", c
PRINT *, "a", a
PRINT *, "v3", v3
PRINT *, "v2", v2
PRINT *, "rs", rs
PRINT *, "d", d

    ! Initialize rs array
    rs(:) = 0

    ! Populate rs array based on ig
    do i = 1, n
        j = ig(i)
        rs(j) = rs(j) + 1
    end do

    ! Initialize arrays
    s(:) = 0.0
    f1m(:) = 0.0
    f1(:) = 0.0
    skmm(:) = 1.0
    skm(:) = 1.0
    v3(:) = 0.0
    v(:) = 0.0
    v2(:, :) = 0.0
    c(:, :) = 0.0
    a(:, :) = 0.0

    ! Begin looping over unique times
    fm = 0.0
    f = 0.0
    ll = 1
    lu = ll
    
   !!!PRINT *, "STARTING CRST=====1"
50  do
        lu = lu + 1
        if (lu > n) exit
        if (y(lu) > y(ll)) exit
    end do
   !!!PRINT *, "STARTING CRST=====2"
    
    lu = lu - 1
    nd1 = 0
    nd2 = 0
    ! d will contain the # in each group censored, failed from
    ! cause 1, and failing from cause 2, at this time
    d(:, :) = 0
   !!!PRINT *, "STARTING CRST=====3"

    do i = ll, lu
        j = ig(i)
        k = m(i)
        d(k, j) = d(k, j) + 1
    end do

   !!!PRINT *, "STARTING CRST=====4"
    nd1 = sum(d(1, :))
    nd2 = sum(d(2, :))

    !if (nd1 == 0 .and. nd2 == 0) PRINT *, "LINE 1539: goto 90"
    if (nd1 == 0 .and. nd2 == 0) goto 90
   !!!PRINT *, "STARTING CRST=====5"

    tr = 0.0
    tq = 0.0

    do i = 1, ng
        if (rs(i) <= 0) cycle
        td = d(1, i) + d(2, i)
        ! skmm is left continuous, and skm right continuous, km est.
        skm(i) = skmm(i) * (rs(i) - td) / rs(i)
        ! f1m is left continuous, and f1 right continuous, cuminc est.
        f1(i) = f1m(i) + (skmm(i) * d(1, i)) / rs(i)
        ! in notation of the Gray 1988 paper, tr is \sum_r\hat{h}_r, and tq is \sum_r R_r
        tr = tr + rs(i) / skmm(i)
        tq = tq + rs(i) * (1.0 - f1m(i)) / skmm(i)
    end do
   !!!PRINT *, "STARTING CRST=====6"

    f = fm + nd1 / tr
    fb = (1.0 - fm)**rho

    a(:, :) = 0.0

!!!PRINT *, "STARTING CRST=====7"
    do i = 1, ng
        if (rs(i) <= 0) cycle
        t1 = rs(i) / skmm(i)
        a(i, i) = fb * t1 * (1.0 - t1 / tr)
        if (a(i, i) /= 0.0) c(i, i) = c(i, i) + a(i, i) * nd1 / (tr * (1.0 - fm))

        do j = i + 1, ng
            if (rs(j) <= 0) cycle
            a(i, j) = -fb * t1 * rs(j) / (skmm(j) * tr)
            if (a(i, j) /= 0.0) c(i, j) = c(i, j) + a(i, j) * nd1 / (tr * (1.0 - fm))
        end do
    end do

   !!!PRINT *, "STARTING CRST=====8"
    ! Make a symmetric matrix
    do i = 2, ng
        k = i - 1
        a(i, 1:k) = a(1:k, i)
        c(i, 1:k) = c(1:k, i)
    end do

    do i = 1, ng1
        if (rs(i) <= 0) cycle
        s(i) = s(i) + fb * (d(1, i) - nd1 * rs(i) * (1.0 - f1m(i)) / (skmm(i) * tq))
    end do

    if (nd1 > 0) then
        do k = 1, ng
            if (rs(k) <= 0) cycle
            t4 = 1.0
            if (skm(k) > 0.0) t4 = 1.0 - (1.0 - f) / skm(k)
            t5 = 1.0
            if (nd1 > 1) t5 = 1.0 - (nd1 - 1) / (tr * skmm(k) - 1.0)
            t3 = t5 * skmm(k) * nd1 / (tr * rs(k))
            v3(k) = v3(k) + t4 * t4 * t3

            do i = 1, ng1
                t1 = a(i, k) - t4 * c(i, k)
                v2(i, k) = v2(i, k) + t1 * t4 * t3

                do j = 1, i
                    l = i * (i - 1) / 2 + j
                    t2 = a(j, k) - t4 * c(j, k)
                    v(l) = v(l) + t1 * t2 * t3
                end do
            end do
        end do
    end if

!!!PRINT *, "STARTING CRST=====4"
    if (nd2 > 0) then
        do k = 1, ng
            if (skm(k) <= 0.0 .or. d(2, k) <= 0) cycle
            t4 = (1.0 - f) / skm(k)
            t5 = 1.0
            if (d(2, k) > 1) t5 = 1.0 - (d(2, k) - 1.0) / (rs(k) - 1.0)
            t6 = rs(k)
            t3 = t5 * ((skmm(k)**2) * d(2, k)) / (t6**2)
            v3(k) = v3(k) + t4 * t4 * t3

            do i = 1, ng1
                t1 = t4 * c(i, k)
                v2(i, k) = v2(i, k) - t1 * t4 * t3

                do j = 1, i
                    l = i * (i - 1) / 2 + j
                    t2 = t4 * c(j, k)
                    v(l) = v(l) + t1 * t2 * t3
                end do
            end do
        end do
    end if

90  if (lu >= n) goto 30

    do i = ll, lu
        j = ig(i)
        rs(j) = rs(j) - 1
    end do

    fm = f
    f1m(:) = f1(:)
    skmm(:) = skm(:)

    ll = lu + 1
    lu = ll
    goto 50

30  l = 0
    do i = 1, ng1
        do j = 1, i
            l = l + 1
            do k = 1, ng
                v(l) = v(l) + c(i, k) * c(j, k) * v3(k)
                v(l) = v(l) + c(i, k) * v2(j, k)
                v(l) = v(l) + c(j, k) * v2(i, k)
            end do
        end do
    end do
!!!PRINT *, "***************************************************************"
!!!PRINT *, "========= end of CRST====="
!!!PRINT *, "...printing outputs..."

!!!PRINT *, "y:"
!!!PRINT *, y

!!!PRINT *, "m:"
!!!PRINT *, m

!!!PRINT *, "ig:"
!!!PRINT *, ig

!!!PRINT *, "s:"
!!!PRINT *, s

!!!PRINT *, "f1m:"
!!!PRINT *, f1m

!!!PRINT *, "f1:"
!!!PRINT *, f1

!!!PRINT *, "skmm:"
!!!PRINT *, skmm

!!!PRINT *, "skm:"
!!!PRINT *, skm

!!!PRINT *, "v:"
!!!PRINT *, v

!!!PRINT *, "c:"
DO i = 1, ng
    !!!PRINT *, c(i, :)
END DO

!!!PRINT *, "a:"
DO i = 1, ng
    !!!PRINT *, a(i, :)
END DO

!!!PRINT *, "v3:"
!!!PRINT *, v3

!!!PRINT *, "v2:"
DO i = 1, ng1
    !!!PRINT *, v2(i, :)
END DO

!!!PRINT *, "rs:"
!!!PRINT *, rs

!!!PRINT *, "d:"
DO i = 0, 2
    !!!PRINT *, d(i, :)
END DO

!!!PRINT *, "***************************************************************"

    return
END SUBROUTINE crst
! =================================================================================


subroutine example_sort(n)
  !use stdlib_sorting, only: sort
  implicit none

  integer, intent(in) :: n
  integer,  allocatable :: array(:)

  array = [5, 4, 3, 1, 10, 4, 9]
  print *, "array"
  print *, array
  !call sort(array)
  print *, "sorted array"
  print *, array   !print [1, 3, 4, 4, 5, 9, 10]
  print *, "input n", n

  ! qsort(arr, n, sizeof(int), compare);
  ! CALL qsort(xSorted, cases, 1, nCases)

! base − It represents pointer to the first element of the array to be sorted.
! nitems − It represents number of element in the array.
! size − It represents size of each element in the array.
! compare − It represent a function pointer to a comparison function that compares two elements.
  !CALL qsort4(array, 7, 1, 7)

end subroutine example_sort