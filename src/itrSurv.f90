!itrSurv.f90 in testpackage/itrSurv/src
MODULE INNERS
!use stdlib_sorting
  IMPLICIT NONE

  PUBLIC

  INTEGER, PARAMETER :: dp = selected_real_kind(15,307)
  INTEGER, PARAMETER :: ng_set2 = 2
  INTEGER, PARAMETER :: nst_set1 = 1

  INTEGER, SAVE :: ERT ! 0/1 1 = use extremely randomized tree ! Nov2024: this is not supported for Phase2RE.
  INTEGER, SAVE :: minEventEnd ! Phase 2.
  !CR: minimum number of people with priority cause events in a node
  !RE: minimum number of people with at least one recurrent event in a node
  INTEGER, SAVE :: minEventSurv ! minimum number of people with a survival event in a node 
  INTEGER, SAVE :: mTry ! maximum number of covariates to try
  INTEGER, SAVE :: n  ! number of sampled cases
  INTEGER, SAVE :: n_surv  ! number of sampled SUBJECTS (needed for Phase2RE bc cases = records)
  INTEGER, SAVE :: nAll ! number of cases (records for Phase2RE)
  INTEGER, SAVE :: nAll_surv ! number of subjects (same as nAll for Phase1 and Phase2CR) (needed for Phase2RE because records != subjects)
  INTEGER, SAVE :: nLevs ! maximum number of levels
  ! to do: delete nodesizeend and change nodesizesurv to just node size. update rest of script
  INTEGER, SAVE :: nodeSizeEnd
  INTEGER, SAVE :: nodeSizeSurv ! minimum number of SUBJECTS in a node (not cases because in isPhase2RE, cases = records, not subjects)
  INTEGER, SAVE :: np ! number of covariates
  INTEGER, SAVE :: nrNodes ! maximum number of nodes in a tree
  INTEGER, SAVE :: nt ! number of time points
  INTEGER, SAVE :: nt_death ! number of time points for survival
  INTEGER, SAVE :: nTree ! number of trees
  INTEGER, SAVE :: replace
  INTEGER, SAVE :: rule ! logrank:1; truncated mean: 2; gray: 3; gray2 (CSH): 4; Q_LR RE:5
  INTEGER, SAVE :: sampleSize
!  INTEGER, SAVE :: sampleSize_surv
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
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: delta_m, id_RE2, row_RE2
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
  REAL(dp), DIMENSION(:), ALLOCATABLE, SAVE :: dt_death
  REAL(dp), DIMENSION(:), ALLOCATABLE, SAVE :: end_tp
  REAL(dp), DIMENSION(:), ALLOCATABLE, SAVE :: surv_tp

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
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: row_RE
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: person_ind

  REAL, DIMENSION(:), ALLOCATABLE, SAVE :: ord_responseAll
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: ord_causeindAll

  ! covariates to be considered for split for sampled cases
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: x
  ! covariates to be considered for split for all cases
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: xAll

  ! TRUE = time differences vector has been allocated
  LOGICAL, SAVE :: dtAllocated = .FALSE.
  LOGICAL, SAVE :: isAlloc_list = .FALSE.
  LOGICAL, SAVE :: isAlloc_ind = .FALSE.
  LOGICAL, SAVE :: isAlloc_sampledArray = .FALSE.

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

      subroutine find_unique(val, final, num_unique) ! (val, final, num_unique)
        implicit none
        integer, intent(in) :: val(:)               ! Input array
        integer, allocatable, intent(out) :: final(:) ! Output array of unique elements
        integer, intent(out) :: num_unique           ! Number of unique elements
        integer :: i, min_val, max_val
        integer, allocatable :: unique(:)

        ! Initialize counters and limits
        i = 0
        min_val = minval(val) - 1   ! Start just below the minimum value
        max_val = maxval(val)       ! Maximum value in the array

        ! Allocate memory for unique values (at most the size of val)
        if (allocated(unique)) deallocate(unique)
        allocate(unique(size(val)))

        ! Loop to find unique elements
        do while (min_val < max_val)
            i = i + 1
            min_val = minval(val, mask=val > min_val)
            unique(i) = min_val
        end do

        ! Number of unique elements found
        num_unique = i

        ! Allocate final array to store only the unique elements
        if (allocated(final)) deallocate(final)
        allocate(final(num_unique))
        final = unique(1:num_unique)
    end subroutine find_unique

subroutine find_unique_subjects(id, unique_id, num_unique)
    implicit none
    integer, intent(in) :: id(:)               ! Input array
    integer, intent(out) :: unique_id(:) ! Output array of unique elements
    integer, intent(out) :: num_unique         ! Number of unique elements
    integer, allocatable :: diff(:), id1(:)    ! Temporary arrays
    integer :: n

    ! Get the size of the input array
    n = size(id)
    !PRINT *, "starting find unique subjects"
    ! Allocate temporary arrays
    allocate(id1(n))
    allocate(diff(n))
    !print *, "test1"
    ! Shift the elements of id to compare adjacent elements
    id1(1:n-1) = id(2:n)
    id1(n) = id(n) + 1   ! Set a value that ensures the last element is different
    !print *, "test2"

    ! Calculate the difference between adjacent elements
    diff = id - id1
    !print *, "test2.1"

    ! Pack the unique elements into unique_id array
    !allocate(unique_id(count(diff /= 0)))  ! Allocate unique_id to the correct size
    unique_id = pack(id, diff /= 0)
    !print *, "test3"
    !PRINT *, unique_id

    ! Get the number of unique elements
    num_unique = size(unique_id)
    !print *, "test3.5"

    ! Deallocate temporary arrays (optional in modern Fortran, but good practice)
    !deallocate(id1)
    !print *, "test4.0"
    !deallocate(diff)
    !print *, "test4"

end subroutine find_unique_subjects


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


! =================================================================================
SUBROUTINE sampleExpandWithoutReplacement(nCases, nSubj, sampleSize, id_RE, nRecords, sampledArray, sampledArray_index)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nCases            ! Total number of records
  INTEGER, INTENT(IN) :: nSubj             ! Total number of unique people (distinct IDs)
  INTEGER, INTENT(IN) :: sampleSize        ! Number of people to sample
  INTEGER, INTENT(IN) :: id_RE(nRecords)   ! Array of IDs (record-level mapping)
  INTEGER, INTENT(IN) :: nRecords          ! Total number of records
  INTEGER, ALLOCATABLE, INTENT(OUT) :: sampledArray(:), sampledArray_index(:)

  INTEGER :: xrand(sampleSize)             ! Sampled indices of unique IDs
  INTEGER, ALLOCATABLE :: uniqueIDs(:)     ! Array of unique IDs from id_RE
  INTEGER :: i, j, k, recCount

  LOGICAL, ALLOCATABLE :: isSampled(:)

  ! Find unique IDs from id_RE
  if (allocated(uniqueIDs)) deallocate(uniqueIDs)
  ALLOCATE(uniqueIDs(nSubj))
  uniqueIDs = 0
  k = 1
  DO i = 1, nRecords
    IF (ALL(uniqueIDs(1:k-1) /= id_RE(i))) THEN
      uniqueIDs(k) = id_RE(i)
      k = k + 1
    END IF
    IF (k > nSubj) EXIT
  END DO
  !PRINT *, "uniqueIDs with size ", size(uniqueIDs)
  !PRINT *, uniqueIDs

  ! Sample indices of unique IDs without replacement
  xrand = sampleWithoutReplace(nSubj, sampleSize)
  !PRINT *, "xrand with size", size(xrand)
  !PRINT *, xrand

  ! Allocate logical array to mark sampled records
  if (allocated(isSampled)) deallocate(isSampled)
  ALLOCATE(isSampled(nRecords))
  isSampled = .FALSE.

  ! Mark records corresponding to sampled unique IDs
  DO i = 1, nRecords
    DO j = 1, sampleSize
      IF (id_RE(i) == uniqueIDs(xrand(j))) THEN
        isSampled(i) = .TRUE.
        EXIT
      END IF
    END DO
  END DO

  ! Count total sampled records
  recCount = COUNT(isSampled)

  ! Allocate output arrays
  if (allocated(sampledArray)) deallocate(sampledArray)
  ALLOCATE(sampledArray(recCount))
  if (allocated(sampledArray_index)) deallocate(sampledArray_index)
  ALLOCATE(sampledArray_index(recCount))

  ! Fill sampledArray and sampledArray_index
  k = 1
  DO i = 1, nRecords
    IF (isSampled(i)) THEN
      sampledArray(k) = id_RE(i)
      sampledArray_index(k) = i
      k = k + 1
    END IF
  END DO

END SUBROUTINE sampleExpandWithoutReplacement
! =================================================================================

! =================================================================================
SUBROUTINE sampleCasesWithoutReplacement(nCases, n, isPhase2RE, personIndex, nRecords, sampledArray)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nCases       ! Total cases (e.g., subjects)
  INTEGER, INTENT(IN) :: n            ! Sample size
  LOGICAL, INTENT(IN) :: isPhase2RE   ! Phase 2 flag
  INTEGER, INTENT(IN) :: personIndex(nRecords) ! Array mapping records to people
  INTEGER, INTENT(IN) :: nRecords     ! Total number of records
  INTEGER, ALLOCATABLE, INTENT(OUT) :: sampledArray(:) ! Output array for sampled indices

  INTEGER :: xrand(n)                 ! Temporary array for sampled people
  INTEGER :: i, j, recCount
  LOGICAL, ALLOCATABLE :: isSampled(:)

  ! Sample people
  xrand = sampleWithoutReplace(nCases, n)
  PRINT *, "xrand with size ", size(xrand)
  PRINT *, xrand

  IF (isPhase2RE) THEN
    ! Allocate logical array for marking sampled records
    ALLOCATE(isSampled(nRecords))
    isSampled = .FALSE.

    PRINT *, "personIndex"
    PRINT *, personIndex
    PRINT *
    ! Mark records corresponding to sampled people
    DO i = 1, n
      PRINT *
      PRINT *, "i:", i
      PRINT *, "xrand(i)", xrand(i)
      isSampled = isSampled .OR. (personIndex == xrand(i))
      PRINT *, "isSampled:", isSampled
      PRINT *
    END DO
    !DO i = 1, n
      !DO j = 1, nRecords
        !IF (personIndex(j) == xrand(i)) THEN
        !  isSampled(j) = .TRUE.
        !END IF
      !END DO
    !END DO

    PRINT *, "isSampled"
    PRINT *, isSampled

    ! Count how many records were marked
    recCount = COUNT(isSampled)
    PRINT *, "recCount #1:", recCount

    ! Allocate sampledArray to hold all selected records
    ALLOCATE(sampledArray(recCount))
    recCount = 1
    DO j = 1, nRecords
      IF (isSampled(j)) THEN
        sampledArray(recCount) = j
        recCount = recCount + 1
      END IF
    END DO

    DEALLOCATE(isSampled)
    PRINT *, "recCount:", recCount

  ELSE
    ! If not Phase 2, simply use xrand
    ALLOCATE(sampledArray(n))
    sampledArray = xrand
  END IF

END SUBROUTINE sampleCasesWithoutReplacement
! =================================================================================


! Identify the optimal split
!   nCases : integer, the number of elements in input casesIn
!   casesIn : integer(:), the indices of the cases in this node
!   nSubj : integer, the number of elements in input subjIn; this is the same as nCases for Phase1/2CR
!   subjIn : integer(:), the indices of the people in this node; this is the same as casesIn for Phase1/2CR
!   nv : integer, the number of covariates to include in search
!   varsIn : integer(:), covariates to include in search
!   splitVar : integer, the index of the selected variable for splitting
!   cutoffBest : real(:), the cutoff (<= go to 'left' node)
!   splitFound : integer, 0 = no split; 1 = found a split
!   casesOut : integer(:), elements of casesIn that go left; subject level index if yes, 0
!   casesOutRE : integer(:), records of casesIn that go left; record level index if yes, 0
!     otherwise
!   casesOut_people_RE : integer(:), records for original subject ID labels that go left.
!   nCuts : integer, the number of cutoff values returned
!   lft : integer, the number of cases in the left node
! CALL tfindSplit(size(ind), ind, indRE, size(pind), pind, splitVar, cutoffBest, splitFound, indOut, indOutRE, nc, lft)
      !CALL tfindSplit(k, size(ind), ind, indRE, &
      !              & size_unique_ind_people_RE, ind_people_RE, &
      !              & sampledArray_index_new, &
      !              & record_ind, record_ind_trt, record_surv_RE, &
      !              & size(pind), pind, splitVar, cutoffBest, &
      !              & splitFound, indOut, indOutRE, &
      !              & indOut_people_RE, indOut_record_RE, nc, lft)

! for RE: this is each record and nCases is number of records
! we need a nCases_person and 
SUBROUTINE tfindSplit(node1, nCases, casesIn, casesInRE, &
                    nSubj, subjIn, &
                    sampledArray_index_new1, &
                    record_ind, record_ind_trt, record_people_RE, &
                    nv, varsIn, splitVar, cutoffBest, splitFound, &
                    casesOut, casesOutRE, casesOut_people_RE, &
                    casesOut_record_RE, nCuts, lft)
  use, intrinsic :: ieee_arithmetic
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: node1
  INTEGER, INTENT(IN) :: nCases
  INTEGER, INTENT(IN) :: nSubj
  INTEGER, DIMENSION(1:nCases), INTENT(IN) :: casesIn
  INTEGER, DIMENSION(1:nCases), INTENT(IN) :: casesInRE
  INTEGER, DIMENSION(1:nCases), INTENT(IN) :: subjIn
  INTEGER, DIMENSION(1:nCases), INTENT(IN) :: sampledArray_index_new1 ! dimension(:)
  INTEGER, DIMENSION(1:nCases), INTENT(IN) :: record_ind
  INTEGER, DIMENSION(1:nCases), INTENT(IN) :: record_ind_trt
  INTEGER, DIMENSION(1:nCases), INTENT(IN) :: record_people_RE
  INTEGER, INTENT(IN) :: nv
  INTEGER, DIMENSION(1:nv), INTENT(IN) :: varsIn
  INTEGER, INTENT(OUT) :: splitVar
  REAL(dp), DIMENSION(1:nLevs), INTENT(OUT) :: cutoffBest
  INTEGER, INTENT(OUT) :: splitFound
  INTEGER, DIMENSION(1:nCases), INTENT(OUT) :: casesOut
  INTEGER, DIMENSION(1:nCases), INTENT(OUT) :: casesOutRE
  INTEGER, DIMENSION(1:nCases), INTENT(OUT) :: casesOut_people_RE
  INTEGER, DIMENSION(1:nCases), INTENT(OUT) :: casesOut_record_RE
  INTEGER, INTENT(OUT) :: nCuts
  INTEGER, INTENT(OUT) :: lft

  LOGICAL, DIMENSION(nCases) :: mk

  INTEGER :: cnt, cnt_m, i, iperson, cc, ikv, j, jj, k, kv, l, nUncensored, nUncensored_m, ptr, rightNode, rightNode_m
  INTEGER :: nnn, iii, jjj, count_nonzero, ix
  INTEGER :: rUnifSet, set, splitLeft, splitLeftFinal, tieCovariate
  INTEGER :: splitLeft_m, splitLeftFinal_m, splitLeft_doloop, splitLeftFinal_doloop
  INTEGER :: splitLeft_loop
  INTEGER :: tieValue, variablesTried
  INTEGER, DIMENSION(1:nCases) :: delta_m_sub, delta_sub
  INTEGER, DIMENSION(1:nCases) :: cases, dSorted, dSorted_m, tcases, subjects!, tsubjects
  INTEGER, DIMENSION(1:nCases) :: recordID_og, recordID_trt, recordID, recordID_new
  INTEGER, DIMENSION(1:nCases) :: personID_og, personID, person_ind_sorted, personID_new !tsubjind
  INTEGER, DIMENSION(:), ALLOCATABLE :: uniqueID, firstIndex, lastIndex
  INTEGER :: uniqueCount, doi, doj, doitmp, dojtmp
  INTEGER :: dostart, doend, numCases
  INTEGER, DIMENSION(1:nCases) :: tcasesSorted, subjSorted, sorted_cases_index
  REAL(dp), DIMENSION(1:nCases) :: sorted_cases
  INTEGER, DIMENSION(1:nv) :: variables
  INTEGER, DIMENSION(:), ALLOCATABLE :: ind, ind_m, indSingles, indSingles_m
  INTEGER, DIMENSION(:), ALLOCATABLE :: leftCases, leftCases_loop, rightCases, rightCases_loop
  INTEGER :: caseR
  INTEGER, DIMENSION(:), ALLOCATABLE :: leftPeople_loop, unique_leftPeople_loop, rightPeople_loop, unique_rightPeople_loop
  INTEGER, DIMENSION(:), ALLOCATABLE :: rightPeople_loop_og, leftPeople_loop_og
  INTEGER :: nleftPeople_loop, nrightPeople_loop
  INTEGER, DIMENSION(:), ALLOCATABLE :: uncensoredIndices, uncensoredIndices_m

  INTEGER, DIMENSION(1:nAll) :: group_cr
  INTEGER, DIMENSION(:), ALLOCATABLE :: group_cr2, cases2

  REAL(dp) :: cutoff, maxValueSplit, maxValueXm, rUnif, valuej, tester3a, tester3b
  REAL(dp) :: valuej_num, valuej_denom
  REAL(dp), DIMENSION(1:nt_death) :: atRiskLeft, atRiskRight, D, denJ, numJ
  REAL(dp), DIMENSION(1:nt_death) :: atRiskLeft_loop, atRiskRight_loop
  REAL(dp), DIMENSION(1:nt_death) :: eventsLeft, eventsRight, eventsLeft_loop, eventsRight_loop
  REAL(dp), DIMENSION(:), ALLOCATABLE :: tempSum
  REAL(dp), DIMENSION(1:nt_death) :: pd1, pd2, pd21, pd22, pd1_loop, pd2_loop 
  REAL(dp), DIMENSION(1:nt_death) :: Rcum ! length is number of failure time points
  REAL(dp), DIMENSION(1:nt) :: atRiskLeft_m, atRiskRight_m, D_m
  REAL(dp), DIMENSION(1:nt) :: atRiskLeft_m_loop, atRiskRight_m_loop
  REAL(dp), DIMENSION(1:nt) :: eventsLeft_m, eventsRight_m,  eventsLeft_m_loop, eventsRight_m_loop
  REAL(dp), DIMENSION(1:nt) :: Rcum_m
  REAL(dp), DIMENSION(1:nt) :: pd1_m, pd2_m, pd1_m_loop, pd2_m_loop
  REAL(dp), DIMENSION(1:nCases) :: xSorted, covar_sorted_RE
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: prl, prr, prl_m, prr_m, pr2l, pr2r 
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: pr2survl, pr2survr, prsurvl, prsurvr
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: pr2surv_sub, pr2_sub, prsurv_sub, pr_sub
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: prl_loop, prr_loop, prl_m_loop, prr_m_loop
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: pr2l_loop, pr2r_loop, pr2survl_loop, pr2survr_loop

  LOGICAL :: randomSplit
  LOGICAL, DIMENSION(:), ALLOCATABLE :: singles, singles_m
  INTEGER, DIMENSION(:), ALLOCATABLE :: unique_uncensoredIndices_m
  REAL(dp), DIMENSION(:), ALLOCATABLE :: x_splitLeft_m_vec, x_splitRight_m_vec
  REAL(dp) :: x_splitLeft_m, x_splitRight_m

  REAL(dp), DIMENSION(1:nt) :: survRE_right, survRE_left
  REAL(dp), DIMENSION(1:nt) :: dRhat_right, dRhat_left
  REAL(dp), DIMENSION(1:nt) :: mu_right, mu_left
  REAL(dp), DIMENSION(1:nt) :: dmu_right, dmu_left

  REAL(dp), ALLOCATABLE, DIMENSION(:, :) :: dPsi_left, dPsi_right
  
  ! Declare a temporary array to hold the result of the multiplication
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: tempArray

  LOGICAL :: are_equal, print_check, all_delta_same
  INTEGER :: index_delta
  LOGICAL :: notEqual
  REAL(dp) :: rnd, random1
  INTEGER ::  testi

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

  !REAL(dp), DIMENSION(1:nCases, 6) :: combArray
  REAL(dp), ALLOCATABLE :: inputVectors(:,:), combArray(:,:)
  
  ! below is for calculating dMi and dMiD
  REAL(dp), DIMENSION(nt) :: tmp_events1, tmp_events2, dlam 
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: dNi, YidR, dMi, dNiD, YidLam, dMiD

  real(dp), dimension(13) :: TESTINGpd1
  INTEGER :: iiii, tmp_i, record_ind_size, irec

  ! below is old code that doesn't initialize
  ! REAL(dp) :: s_set(ng_set2-1), vs_set(ng_set2-1, ng_set2-1)
  ! REAL(dp) :: v_set(ng_set2*(ng_set2-1)/2), st_set(ng_set2-1), vt_set(ng_set2*(ng_set2-1)/2)
  ! REAL(dp), DIMENSION(ng_set2*(4+3*ng_set2)) :: wk_set
  ! integer, dimension(4*ng_set2) :: iwk_set

  EXTERNAL :: rnd
  are_equal = .TRUE.
  print_check = .FALSE.

  !IF (isPhase2RE) WRITE(*,'(/,A)') '============================ tfindSplit ============================'

  IF (isPhase2RE) THEN
  IF (print_check) THEN
    PRINT *, "subjIn with size", size(subjIn)
    PRINT *, subjIn
    PRINT *, "casesIn with size", size(casesIn)
    PRINT *, casesIn
    PRINT *, "casesInRE with size", size(casesInRE)
    PRINT *, casesInRE
    PRINT *, "record_ind with size", size(record_ind)
    PRINT *, record_ind
    PRINT *, "record_ind_trt with size", size(record_ind_trt)
    PRINT *, record_ind_trt
    PRINT *, "record_people_RE with size", size(record_people_RE)
    PRINT *, record_people_RE
    PRINT *, "new Sampled Array Indices with size", size(sampledArray_index_new1)
    PRINT *, sampledArray_index_new1
    WRITE(*,'(/,A)') '============================ tfindSplit ============================'
    !PRINT *, "******************** tfindSplit ********************"
    PRINT *, "nCases:", nCases, "and nSubj:", nSubj
    PRINT *, "nodeSizeEnd:", nodeSizeEnd
    PRINT *, "nodeSizeSurv:", nodeSizeSurv
  END IF
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
  casesOutRE = 0
  casesOut_people_RE = 0
  casesOut_record_RE = 0

  ! set number of cutoffs to 0
  nCuts = 0

  ! indices of variables to be explored
  variables = varsIn

  ! location of last available parameter in variables
  ptr = nv

  ! tracks the number of variables tried
  variablesTried = 0

  tcases = (/(i,i=1,nCases)/)

  ! Terminal node criteria     
  IF (isPhase2RE) THEN
    IF (nSubj < 2*nodeSizeSurv .OR. nCases < 2*nodeSizeEnd) RETURN ! if number of subjects is less than 2*min node size for people, then its terminal node
  ELSE
    ! Ensure nodeSizeSurv equals nodeSizeEnd in non-Phase2RE cases
    IF (nodeSizeSurv .NE. nodeSizeEnd) THEN 
      PRINT *, "Error: nodeSizeSurv and nodeSizeEnd mismatch so stopping"
      STOP
    ELSE
      ! nodeSizeSurv = nodeSizeEnd for isPhase1 and isPhase2CR
      IF (nCases < 2*nodeSizeEnd) RETURN ! if number of cases is less than 2*min node size for cases then this is terminal node
    END IF
  END IF

  !IF (isPhase2RE) PRINT *, "*** Starting covariate do-loop LINE 670 ***"
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
      IF (.NOT. isPhase2RE) THEN
        CALL getCovariate(nCases, casesIn, kv, xSorted)
      ELSE
        !PRINT *, "RE: ================================================================"
        !PRINT *, "get covar"
        PRINT *, "WARNING: This may not be coded up yet for RE"
        CALL getCovariate(nCases, casesIn, kv, covar_sorted_RE)
        !PRINT *, "shape x:", shape(x)
        !PRINT *, "sahpe covar_sorted_RE:", shape(covar_sorted_RE)
      END IF
    ELSE
      !PRINT *, "nCat(kv) is either 0 = continuous or 1 = ordered factors"
      IF (.NOT. isPhase2RE) THEN
        xSorted = x(casesIn,kv)
      ELSE ! Phase2RE- use record_ind
        !PRINT *, "RE: ================================================================"
        !PRINT *, "not get covar"
        covar_sorted_RE = x(casesInRE,kv) !record_ind is old
        !PRINT *, "shape x:", shape(x)
        !PRINT *, "shape covar_sorted_RE:", shape(covar_sorted_RE)
        !PRINT *, "casesInRE with size:", size(casesInRE)
        !PRINT *, casesInRE
        delta_sub = delta(casesInRE)
        delta_m_sub = delta_m(casesInRE)
        IF (print_check) THEN
          PRINT *, "shape delta:", shape(delta)
          PRINT *, "shape delta_sub:", shape(delta_sub)
          PRINT *, "kv=", kv
          PRINT *, "x with size", size(x)
          PRINT *, x
          PRINT *, size(covar_sorted_RE)
          PRINT *, covar_sorted_RE
          PRINT *
        END IF
      END IF
    END IF

    !IF (isPhase2RE) THEN
      ! this is for pr2, pr, prsurv, pr2surv
      !PRINT *, "```````````````````````````"
      !PRINT *, "shape of pr2", shape(pr2)
      !PRINT *, "shape of pr", shape(pr)
      !PRINT *, "record_ind_trt with size:", size(record_ind_trt)
      !PRINT *, record_ind_trt
      !pr_sub = pr!(record_ind_trt,:)
      !prsurv_sub = prsurv!(record_ind_trt,:)
      !pr2_sub = pr2!(record_ind_trt,:)
      !pr2surv_sub = pr2surv!(record_ind_trt,:)
      !PRINT *, "shape of delta", size(delta)
      !delta_sub = delta
      !delta_m_sub = delta_m
      !PRINT *, "shape of pr2_sub", shape(pr2_sub)
      !PRINT *, "shape of pr_sub", shape(pr_sub)
      !PRINT *, "```````````````````````````"
    !END IF

    !================================================================
    !================================================================
    IF (isPhase1 .OR. isPhase2CR) THEN
      cases = casesIn
      !IF (isPhase1) PRINT *, "first casesIn with size", size(casesIn)
      !IF (isPhase1) PRINT *, casesIn
      !PRINT *, "pre-sorting cases"
      !PRINT *, cases
      ! sort the covariate and track the indices
      CALL qsort4(xSorted, cases, 1, nCases)
      !PRINT *, "Phase1/2CR: post-sorting cases with size", size(cases)
      !PRINT *, cases
      
      !PRINT *, "Combined Details:"
      !PRINT *, "------------------"
      !PRINT *, size(xSorted)
      !PRINT *, size(cases)
      !PRINT *, " xSorted     cases "
      !PRINT *, "------------------"
      !DO tmp_i = 1, SIZE(cases)
      !PRINT '(I12, 3X, F6.3)', cases(tmp_i), xSorted(tmp_i)
      !END DO
      !PRINT *, "------------------"
      !PRINT *

      ! sort event indicator data accordingly
      ! Phase 1: overall survival (delta = event indicator from any cause)
      dSorted = delta(cases) 
      IF (isPhase2CR) THEN
        ! Phase 2: CR (delta_m = indicator for event from cause m)
        dSorted_m = delta_m(cases) 
      END IF
    !================================================================
    !================================================================
    ELSE IF (isPhase2RE) THEN
      !PRINT *, "i thought casesIn was:", size(casesIn)
      !PRINT *, casesIn
      tcasesSorted = sampledArray_index_new1 ! tcasesSorted = tcases
      !PRINT *, "now setting tcasesSorted (before sorting) with inital tcasesSorted size:  ", SIZE(tcasesSorted)
      !PRINT *, tcasesSorted
      !PRINT *, "pre-sort covar:", size(covar_sorted_RE)
      !PRINT *, covar_sorted_RE(1:5)
      CALL qsort4(covar_sorted_RE, tcasesSorted, 1, nCases)
      !PRINT *, "post-sort covar:", size(covar_sorted_RE)
      !PRINT *, covar_sorted_RE(1:5)
      !PRINT *, "sorted tcasesSorted size:  ", SIZE(tcasesSorted)
      !PRINT *, tcasesSorted
      
      IF (print_check) THEN
        ! Print sizes before executing the main block of code
        PRINT *, "Initial array sizes:"
        PRINT *, "pr2 has size:       ", SIZE(pr2, 1), " x ", SIZE(pr2, 2)
        PRINT *, "pr2surv has size:   ", SIZE(pr2surv, 1), " x ", SIZE(pr2surv, 2)
        PRINT *, "delta has size:     ", SIZE(delta)
        PRINT *, "delta_m has size:   ", SIZE(delta_m)
        PRINT *, "tcasesSorted size:  ", SIZE(tcasesSorted)
        PRINT *, "subjIn has size:", size(subjIn)
        PRINT *, "casesIn has size:", size(casesIn)
        PRINT *, "delta(tcasesSorted) size:", SIZE(delta(tcasesSorted))
        PRINT *, "covar_sorted_RE with size", size(covar_sorted_RE)
        PRINT *, "---------------------------------------------"
      END IF

      ! below is old - delete
      ! i think what i need to do is edit casesIn(), delta(), and delta_m() to be soemthing not tcasesSorted
            !CALL group_and_sort(covar_sorted_RE, subjIn(tcasesSorted), &
            !              tcasesSorted, casesIn(tcasesSorted), &
            !              delta(tcasesSorted), delta_m(tcasesSorted), &
            !              size(covar_sorted_RE), combArray)
            !xSorted = combArray(:, 1) ! xSorted
            !personID_og = combArray(:, 2) ! personID_og
            !recordID = combArray(:, 3) ! recordID
            !personID = combArray(:, 4) ! personID
            !dSorted = combArray(:, 5) ! dSorted
            !dSorted_m = combArray(:, 6) ! dSorted_m
            !CALL create_new_vector(personID_og, size(personID_og), personID_new)

        if (allocated(inputVectors)) then
            deallocate(inputVectors)
        end if
        ! Allocate the inputVectors array
        ALLOCATE(inputVectors(size(covar_sorted_RE), 7))
        ! Combine vectors into the inputVectors array
        inputVectors(:, 1) = covar_sorted_RE
        inputVectors(:, 2) = REAL(subjIn(tcasesSorted), dp) ! Cast integers to REAL(dp)
        inputVectors(:, 3) = REAL(record_people_RE(tcasesSorted), dp) 
        inputVectors(:, 4) = REAL(casesIn(tcasesSorted), dp)
        inputVectors(:, 5) = REAL(casesInRE(tcasesSorted), dp) !REAL(sampledArray_index_new1(tcasesSorted), dp)
        inputVectors(:, 6) = REAL(delta_sub(tcasesSorted), dp)
        inputVectors(:, 7) = REAL(delta_m_sub(tcasesSorted), dp)

        if (allocated(combArray)) then
          deallocate(combArray)
        end if
        ! Allocate the combinedArray (optional; can also be handled in subroutine)
        ALLOCATE(combArray(size(covar_sorted_RE), 7))
        CALL group_and_sort(inputVectors, size(covar_sorted_RE), combArray)

          
    !PRINT *, "Size of delta_sub: ", SIZE(delta_sub)
    !PRINT *, delta_sub
    !PRINT *, "Size of delta: ", SIZE(delta)
    !PRINT *, delta
    !PRINT *, "Size of delta_sub(tcasesSorted): ", SIZE(delta_sub(tcasesSorted))
    !PRINT *, delta_sub(tcasesSorted)
    !PRINT *, "Size of delta_m_sub: ", SIZE(delta_m_sub)
    !PRINT *, delta_m_sub
    !PRINT *, "Size of delta_m: ", SIZE(delta_m)
    !PRINT *, delta_m
    !PRINT *, "Size of delta_m_sub(tcasesSorted): ", SIZE(delta_m_sub(tcasesSorted))
    !PRINT *, delta_m_sub(tcasesSorted)

    !PRINT *, "FINAL SORTED ARRAY:"
    !PRINT *, "xSorted    personID_og     recordID_og     personID     recordID     dSorted      dSorted_m"
    !DO tmp_i = 1, 3 !SIZE(combArray, 1)
    !    WRITE(*, '(1X, F10.4, *(1X, I10))') combArray(tmp_i, 1), &
    !        (INT(combArray(tmp_i, j)), j = 2, SIZE(combArray, 2))
    !END DO

    xSorted = combArray(:, 1) ! xSorted
    personID_og = combArray(:, 2) ! personID_og
    recordID_og = combArray(:, 3) ! recordID_og
    personID = combArray(:, 4) ! 
    recordID = combArray(:, 5) ! personID within this sample
    dSorted = combArray(:, 6) ! dSorted
    dSorted_m = combArray(:, 7) ! dSorted_m
    CALL create_new_vector(personID_og, size(personID_og), personID_new)
    CALL create_new_vector(recordID, size(recordID), recordID_new)

      !PRINT *, "dimensions of combArray:", shape(combArray)

      !PRINT *, "delta has size:     ", SIZE(delta)
      !PRINT *, delta

      !PRINT *, "delta_m has size:   ", SIZE(delta_m)
      !PRINT *, delta_m

      !PRINT *, "delta_sub has size:     ", SIZE(delta_sub)
      !PRINT *, delta_sub

      !PRINT *, "delta_m_sub has size:   ", SIZE(delta_m_sub)
      !PRINT *, delta_m_sub

      !PRINT *, "dSorted size:", SIZE(dSorted)
      !PRINT *, dSorted

      !PRINT *, "dSorted_m size:", SIZE(dSorted_m)
      !PRINT *, dSorted_m

        IF (print_check) THEN
          PRINT *, "sampledArray_index_new1 with size:", size(sampledArray_index_new1)
          PRINT *, sampledArray_index_new1
          PRINT *, "sampledArray_index_new1(tcasesSorted) with size:", size(sampledArray_index_new1(tcasesSorted))
          PRINT *, sampledArray_index_new1(tcasesSorted)
          PRINT *, "tcasesSorted with size:", size(tcasesSorted)
          PRINT *, size(tcasesSorted)
          PRINT *, "Combined Details:"
          PRINT *, "------------------"
          PRINT *, "covar_sorted_RE has size: ", size(covar_sorted_RE)
          PRINT *, "subjIn has size: ", size(subjIn)
          PRINT *, "tcasesSorted has size: ", size(tcasesSorted)
          PRINT *, "delta has size:     ", SIZE(delta)
          PRINT *, "delta_m has size:   ", SIZE(delta_m)
          PRINT *, "delta_sub has size:     ", SIZE(delta_sub)
          PRINT *, "delta_m_sub has size:   ", SIZE(delta_m_sub)
          PRINT *, "personID_og has size: ", size(personID_og)
          PRINT *, "recordID has size: ", size(recordID)
          PRINT *, "personID_new has size: ", size(personID_new)
          PRINT *, "pr2 has size: ", shape(pr2)
          PRINT *, "pr2_sub has size: ", shape(pr2_sub)
          PRINT *, "pr has size: ", shape(pr)
          !PRINT *, "pr_sub has size: ", shape(pr_sub)
          PRINT *
          PRINT *, "xSorted  personID_og  recordID_og  recordID_trt  personID_new  recordID_new  dSorted  dSorted_m"
          PRINT *, "----------------------------------------------------------------------------------------------------------"
          DO tmp_i = 1, 5 !SIZE(personID_og)
          PRINT '(F6.3, I12, 1X, I12, 1X, I12, 1X, I12, 2X, I12, I12, I12)', &
            & xSorted(tmp_i), personID_og(tmp_i), &
            & recordID_og(tmp_i), recordID_trt(tmp_i), &
            & personID_new(tmp_i), recordID_new(tmp_i), &
            dSorted(tmp_i), dSorted_m(tmp_i)
          END DO
          PRINT *, "----------------------------------------------------------------------------------------------------------"
          PRINT *
          !PRINT *, "delta sub:", size(delta_sub)
          !PRINT *, delta_sub
          !PRINT *, "delta sub endpoint:", size(delta_m_sub)
          !PRINT *, delta_m_sub
          !PRINT *, "delta_sub(tcasesSorted):", size(delta_sub(tcasesSorted))
          !PRINT *, delta_sub(tcasesSorted)
          !PRINT *, "delta_m_sub(tcasesSorted):", size(delta_m_sub(tcasesSorted))
          !PRINT *, delta_m_sub(tcasesSorted)
          PRINT *, "delta has size:     ", SIZE(delta)
          PRINT *, "delta_m has size:   ", SIZE(delta_m)
          PRINT *, "delta_sub has size:     ", SIZE(delta_sub)
          PRINT *, "delta_m_sub has size:   ", SIZE(delta_m_sub)
          PRINT *, "dSorted size:", size(dSorted)
          PRINT *, "dSorted_m size:", size(dSorted_m)
          PRINT *, "subjIn with size:", size(subjIn)
          PRINT *, subjIn
          PRINT *, "row      delta      delta_m"
          PRINT *, "--------------------------------------"
          DO tmp_i = 1, SIZE(delta)
          PRINT '(I4, 1X, I4, 1X, I4)', &
            & tmp_i, delta(tmp_i), &
            & delta_m(tmp_i)
          END DO
          PRINT *, "--------------------------------------"
          PRINT *
          PRINT *, "row    personID_og    tcasesSorted    dSorted    dSorted_m"
          PRINT *, "--------------------------------------------------------------"
          DO tmp_i = 1, SIZE(personID_og)
          PRINT '(I4, 1X, I10, 1X, I7, 1X, I4, 1X, I4)', &
            & tmp_i, personID_og(tmp_i), &
            & tcasesSorted(tmp_i), dSorted(tmp_i), dSorted_m(tmp_i)
          END DO
          PRINT *, "--------------------------------------------------------------"
          PRINT *
        END IF

        IF (print_check) THEN

            !IF (size(prsurv_sub,1) .NE. size(prsurv,1)) THEN
              ! Print sizes after executing the main block of code
              !PRINT *, "Derived array sizes:"
              !PRINT *, "xSorted size:       ", SIZE(xSorted)
              !PRINT *, "personID_og size:   ", SIZE(personID_og)
              !PRINT *, "recordID size:      ", SIZE(recordID)
              !PRINT *, "dSorted size:       ", SIZE(dSorted)
              !PRINT *, "dSorted_m size:     ", SIZE(dSorted_m)
              !PRINT *, "prsurv_sub size:    ", SIZE(prsurv_sub, 1), " x ", SIZE(prsurv_sub, 2)
              !PRINT *, "pr_sub size:        ", SIZE(pr_sub, 1), " x ", SIZE(pr_sub, 2)
              !PRINT *, "pr2surv_sub size:   ", SIZE(pr2surv_sub, 1), " x ", SIZE(pr2surv_sub, 2)
              !PRINT *, "pr2_sub size:       ", SIZE(pr2_sub, 1), " x ", SIZE(pr2_sub, 2)
              PRINT *, "prsurv size:    ", SIZE(prsurv, 1), " x ", SIZE(prsurv, 2)
              PRINT *, "pr size:        ", SIZE(pr, 1), " x ", SIZE(pr, 2)
              PRINT *, "pr2surv size:   ", SIZE(pr2surv, 1), " x ", SIZE(pr2surv, 2)
              PRINT *, "pr2 size:       ", SIZE(pr2, 1), " x ", SIZE(pr2, 2)
              !PRINT *, "---------------------------------------------"
            !END IF

        ! Assuming all arrays have the same size
        PRINT *, " personID_og   recordID   personID_new"
        PRINT *, "------------------------------------"
        DO tmp_i = 1, SIZE(personID_og)
          PRINT '(I12, I12, I12)', personID_og(tmp_i), recordID(tmp_i), personID_new(tmp_i)
          END DO
        PRINT *, "------------------------------------"
      END IF

    !================================================================
    !================================================================
    ELSE 
      PRINT *, "ERROR: Phase not recognized. Check isPhase1, isPhase2CR, or isPhase2RE, so stopping."
      STOP
    END IF
    !================================================================
    !================================================================

    ! ******************** splitBoundaries ********************
    ! identify minimum cases for left and right splits 
    ! based on minimum uncensored cases, minimum node size, 
    ! and assurance that all equal valued cases are included 
    ! in the minimum nodes
    !PRINT *, "******************** splitBoundaries ********************"
    rUnif = 0.d0
    rUnifSet = -1 !-1 means doesn't satisify requirements, move to next covariate. this is default and change to 0 if satisfied.

    ! tcases and tsubjects, and delta and delta_m are og in order indices
    ! cases and subjects, and dSorted and dSorted_m have been sorted with xSorted
    ! dSorted == Phase 1 and is for SURVIVAL for Phase2CR and Phase2RE.

    ! Do below for both isPhase1 and isPhase2CR
    ! SURVIVAL cases that are not-censored
    uncensoredIndices = pack(tcases, dSorted .EQ. 1) 
    nUncensored = size(uncensoredIndices)

    IF (isPhase2CR) THEN ! Phase1/2CR
    ! ENDPOINT cases that are not-censored ! not used for Phase1
    ! because one per person
      uncensoredIndices_m = pack(tcases, dSorted_m .EQ. 1) 
      nUncensored_m = size(uncensoredIndices_m)
    END IF

    IF (isPhase2RE) THEN
      ! because one per record
      uncensoredIndices_m = pack(personID_new, dSorted_m .EQ. 1)
      !IF (ALLOCATED(unique_uncensoredIndices_m)) THEN
        !PRINT *, "deallocating unique_uncensoredIndices_m"
        !DEALLOCATE(unique_uncensoredIndices_m)
      !END IF
      call find_unique(uncensoredIndices_m, &
      unique_uncensoredIndices_m, nUncensored_m)
    END IF

    ! gives cases indices for those who have delta = 1 (event)
    ! tcases is 1,2,...,nCases (so its index for current version)
    ! RE: uncensored means recurrent event (so censored = either death or censored)
    IF (print_check) THEN
        PRINT *, "Death Event Indices: ", uncensoredIndices
        PRINT *, "Number of people with Deaths: ", nUncensored
        PRINT *, "nSubj:", nSubj
        PRINT *, "nCases:", nCases
            IF (isPhase2RE) THEN
        PRINT *, "Recurrent Event Indices: ", uncensoredIndices_m
        PRINT *, "Number of people with Recurrent Events: ", nUncensored_m
        PRINT *, "Minimum # people with Recurrent Events required per node: ", minEventEnd
        PRINT *, "Minimum # people with Deaths required per node: ", minEventSurv
        PRINT *, "Minimum number of people required per node (nodeSize): ", nodeSizeSurv
      END IF
    END IF


    IF (isPhase1 .OR. isPhase2CR) THEN
        IF (nSubj .NE. nCases) THEN
            PRINT *, "ERROR: The number of subjects does not match the number of cases."
            PRINT *, "nSubj:", nSubj, "and nCases:", nCases
            STOP
        END IF
    END IF

    ! FOR CR: since PC is a subset of overall survival, we don't really need extra specification.
    !! Able to split and satisfy minimum number of events in each node
    IF (isPhase1) THEN
      ! If too few uncensored cases to meet minimum number of people with survival events, CYCLE
      IF (nUncensored .LT. (minEventSurv * 2)) CYCLE
    ELSE ! Phase 2
      ! If too few uncensored cases to meet minimum number of people with ENDPOINT events (PC event or at least one RE), CYCLE
      IF (nUncensored_m .LT. (minEventEnd * 2)) CYCLE
    END IF

    ! ============================================
    ! Below is for Phase 1
    ! ============================================
    IF (isPhase1) THEN
      ! we want to do below for all Phases 
      ! cases to left include all indices up to and including minEventSurv case
      ! must have at least nodeSizeSurv cases
      splitLeft = max(uncensoredIndices(minEventSurv), nodeSizeSurv)
      ! move splitLeft up to include cases with equivalent values of x
      splitLeft = count(xSorted .LE. (xSorted(splitLeft) + 1e-8))

      ! cases to right
      ! include all indices down to and including nUncensored - minEventSurv + 1 case
      ! must have at least nodeSizeSurv cases
      rightNode = min(uncensoredIndices(nUncensored - minEventSurv + 1), &
                    & nSubj - nodeSizeSurv + 1)
      ! nUncensord - minEventSurv + 1 = whatever is left from people who died minus min people required for death

      ! move rightNode down to include cases with equivalent values of x
      ! splitLeftFinal is the last possible case for the left node
      splitLeftFinal = count(xSorted .LT. xSorted(rightNode))
      ! if the splitLeft index is above the splitLeftFinal index cycle,
      ! split is not possible
      IF (splitLeft .GT. splitLeftFinal) CYCLE
    END IF

    ! ============================================
    ! Below is for Phase 2 (Endpoint: CR)
    ! ============================================
    IF (isPhase2CR) THEN
      ! cases to left include all indices up to and including minEventSurv case
      ! must have at least nodeSizeSurv cases
      splitLeft_m = max(uncensoredIndices_m(minEventEnd), nodeSizeEnd)
      ! move splitLeft up to include cases with equivalent values of x
      splitLeft_m = count(xSorted .LE. (xSorted(splitLeft_m) + 1e-8))
    ! cases to right
    ! include all indices down to and including nUncensored - minEvent + 1 case
    ! must have at least nodeSize cases
      rightNode_m = min(uncensoredIndices_m(nUncensored_m - minEventEnd + 1), &
                    & nSubj - nodeSizeEnd + 1)
      ! move rightNode down to include cases with equivalent values of x
      ! splitLeftFinal is the last possible case for the left node
      splitLeftFinal_m = count(xSorted .LT. xSorted(rightNode_m))
      ! if the splitLeft index is above the splitLeftFinal index cycle,
      ! split is not possible
      IF (splitLeft_m .GT. splitLeftFinal_m) CYCLE
    END IF

    ! below this section is problematic

    ! ============================================
    ! Below is for Phase 2 (Endpoint: RE)
    ! ============================================
    IF (isPhase2RE) THEN
      ! cases to left include all indices up to and including minEventSurv case
      ! must have at least nodeSizeSurv cases
      splitLeft_m = max(unique_uncensoredIndices_m(minEventEnd), nodeSizeEnd)
      ! move splitLeft up to include cases with equivalent values of x

      ! original: !splitLeft_m = count(xSorted .LE. (xSorted(splitLeft_m) + 1e-8)) ! original
      x_splitLeft_m_vec = pack(xSorted, personID_new .EQ. splitLeft_m)!we want xSorted(personID_new = splitLeft_m)

      x_splitLeft_m = x_splitLeft_m_vec(1) ! take one of them
      splitLeft_m = count(xSorted .LE. (x_splitLeft_m + 1e-8)) ! 
    
      ! cases to right
      ! include all indices down to and including nUncensored - minEvent + 1 case
      ! must have at least nodeSize cases
      rightNode_m = min(unique_uncensoredIndices_m(nUncensored_m - minEventEnd + 1), &
      & nSubj - nodeSizeEnd + 1)
      ! move rightNode down to include cases with equivalent values of x
      ! splitLeftFinal is the last possible case for the left node
      x_splitRight_m_vec = pack(xSorted, personID_new .EQ. rightNode_m)
      x_splitRight_m = x_splitRight_m_vec(1)
      splitLeftFinal_m = count(xSorted .LT. x_splitRight_m)

      ! if the splitLeft index is above the splitLeftFinal index cycle,
      ! split is not possible
      IF (splitLeft_m .GT. splitLeftFinal_m) CYCLE

    END IF

    ! above this section is problematic
    IF (isPhase2) THEN
      splitLeft = splitLeft_m
      splitLeftFinal = splitLeftFinal_m
    END IF

    rUnifSet = 0

    ! ============================================
    ! Below is for Phase 1 or Phase2CR
    ! ============================================
    IF (isPhase1 .OR. isPhase2CR) THEN
      IF ((.NOT. randomSplit) .AND. (ERT .EQ. 1)) THEN
        !PRINT *, "randomSplit:", randomSplit
        !PRINT *, "ERT:", ERT
      
        !************* getUniformSplit *****************
        !PRINT *, "************* getUniformSplit *****************"

        !! split based on extremely randomized trees and uniform split inputs

        ! if ERT splitLeft = splitLeftFinal, which is the splitting point
        ! and cutoff is set to random value or mid-point depending on uniERT

        !! extremely randomized tree methods used

        IF (uniformSplit .EQ. 0) THEN
          !PRINT *, "uniformSplit = 0"

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

        ELSE IF (uniformSplit .EQ. 1) THEN
          !PRINT *, "uniformSplit = 1"
          ! randomly select a value in the range of values that satisfy the allowed cases in the left/right nodes
          random1 = rnd(0.d0, 1.d0)
          rUnif = random1 * (xSorted(splitLeftFinal+1) - &
                & xSorted(splitLeft)) + xSorted(splitLeft)
          rUnifSet = 1
          
          ! identify the first case that splits to the right of this value
          ! the preceding case is the last case to the left node
          splitLeftFinal = nSubj - count(xSorted > rUnif)

          ! the last required case in the left node is now splitLeftFinal
          splitLeft = splitLeftFinal

        END IF

      END IF ! end of IF ((.NOT. randomSplit) .AND. (ERT .EQ. 1))
    END IF ! end of isPhase1/2CR

    IF (isPhase2RE) THEN
      IF (ERT .EQ. 1) THEN
        PRINT *, "Currently, ERT is not supported for Recurrent Events Endpoint. Please select ERT = 0."
        STOP
      END IF
    END IF

    IF (splitLeft <= 0) THEN
      PRINT *, "Line 718 ERROR: splitLeft <= 0"
      PRINT *, "splitLeft:", splitLeft
      STOP
    END IF

    ! -1 is returned if cannot satisfy minimum requirements for nodes
    ! cycle to next covariate
    !PRINT *, "rUnifSet:", rUnifSet
    IF (rUnifSet .EQ. -1) CYCLE

    ! increment the number of covariates that have been explored
    variablesTried = variablesTried + 1

    !***************** maxValue ***************
    ! set initial values for outputs
    set = 0
    maxValueXm = 0.d0
    cutOff = 0.d0

    !*******************************************************************
    !*******************************************************************
    !*******************************************************************
    ! SPLITTING IS DONE (ABOVE). NOW WE LOOK AT CASES (SUBJECTS) IN LEFT AND RIGHT DAUGHTER NODES
    !*******************************************************************
    !*******************************************************************
    !*******************************************************************
    !IF (isPhase2RE) PRINT *, "LINE 1264: SPLITTING IS DONE"
    if (splitLeft <= 0) then
      PRINT *, "ERROR: splitLeft is 0 or negative:", splitLeft
      STOP
    end if

    ! ============================================
    ! Below is for Phase 1 or Phase2CR
    ! ============================================
    if (isPhase1 .OR. isPhase2CR) THEN
      leftCases = cases(1:(splitLeft-1))
      rightCases = cases(splitLeft:nSubj)
    END IF

    !IF (isPhase2RE .AND. node1 .EQ. 2 .AND. i .EQ. 1) THEN
        !PRINT *, "------------------------------------------------------------"
        !PRINT *, "Testing malloc error: i =", i, ", node1 =", node1
        !PRINT *, "------------------------------------------------------------"
        !STOP "Program terminated due to testing malloc error."
    !END IF
    
    ! ============================================
    ! Below is for Phase2RE
    ! ============================================
    IF (isPhase2RE) THEN
      !PRINT *, "LINE 1346 splitLeft: ", splitLeft
      leftCases = recordID_new(1:(splitLeft)) ! 12/3/24: make this splitLeft not splitLeft-1 (diff from Phase1/2CR)
      rightCases = recordID_new(splitLeft+1:nCases) ! recordID(splitLeft+1:nCases)
      !PRINT *, "LINE 1329"
      !leftPeople = personID_new(1:(splitLeft))
      !rightPeople = personID_new(splitLeft+1:nCases)
    END IF

    IF (isPhase2RE) THEN
    IF (print_check) THEN
      ! Verify the sizes of dSorted_m and nt
      PRINT *, "Checking variable sizes for validation:"
      PRINT *, "Size of dSorted_m:", SIZE(dSorted_m)
      PRINT *, "Value of splitLeft:", splitLeft
      IF (splitLeft > SIZE(dSorted_m)) THEN
          PRINT *, "Error: splitLeft exceeds the size of dSorted_m."
          STOP "Invalid size for dSorted_m."
      END IF
      PRINT *, "Value of nt:", nt
      IF (nt <= 0) THEN
          PRINT *, "Error: nt must be greater than 0."
          STOP "Invalid value for nt."
      END IF
      PRINT *, "Size checks passed."

      ! Validate slicing for leftCases and rightCases
      PRINT *, "Validating indices for leftCases and rightCases:"
      ! Check that splitLeft is within the valid range for recordID_new
      IF (splitLeft < 1 .OR. splitLeft > SIZE(recordID_new)) THEN
          PRINT *, "Error: splitLeft =", splitLeft, "is out of bounds for recordID_new, size =", SIZE(recordID_new)
          STOP "Invalid splitLeft index for recordID_new."
      END IF
      ! Check that nCases is within the valid range for recordID_new
      IF (nCases < splitLeft + 1 .OR. nCases > SIZE(recordID_new)) THEN
          PRINT *, "Error: nCases =", nCases, "is out of bounds or invalid for recordID_new, size =", SIZE(recordID_new)
          STOP "Invalid nCases index for recordID_new."
      END IF
      ! Check the size of leftCases
      PRINT *, "Size of leftCases will be:", splitLeft
      IF (splitLeft > SIZE(recordID_new)) THEN
          PRINT *, "Error: Calculated size of leftCases exceeds bounds of recordID_new."
          STOP "Invalid size for leftCases."
      END IF
      ! Check the size of rightCases
      PRINT *, "Size of rightCases will be:", nCases - splitLeft
      IF (splitLeft + 1 > nCases) THEN
          PRINT *, "Error: splitLeft+1 exceeds nCases. Cannot compute rightCases."
          STOP "Invalid size for rightCases."
      END IF
      PRINT *, "Index validation for leftCases and rightCases passed."

      PRINT *, "nSubj: ", nSubj
      PRINT *, "nCases: ", nCases
      PRINT *, "# death timepoints:", nt_death
      PRINT *, "# RE timepoints:", nt
      ! Print the first 5 elements of the prsurv_sub array in a readable way
      PRINT *, "prsurv_sub array for the first 10 elements of leftCases:"
      DO testi = 1, MIN(10, SIZE(leftCases))  ! Limit to the first 10 elements or the size of leftCases
          PRINT *, "Person ID ", leftCases(testi)," with delta = ", real(dSorted(testi)), ": "
          DO j = 1, SIZE(prsurv_sub, 2)  ! Loop over each timepoint
              ! Break the print statement into two parts to avoid line truncation
              PRINT '(A, I3, A, F6.1, A, F6.1)', "prsurv_sub for timepoint ", j, ": ", real(prsurv_sub(leftCases(testi), j))
          END DO
      END DO

      ! Print the first 5 elements of the pr array in a readable way
      PRINT *, "prsurv array for the first 10 elements of leftCases:"
      DO testi = 1, MIN(10, SIZE(leftCases))  ! Limit to the first 10 elements or the size of leftCases
          PRINT *, "Person ID ", leftCases(testi)," with delta = ", real(dSorted_m(testi)), ": "
          DO j = 1, SIZE(prsurv, 2)  ! Loop over each timepoint
              ! Break the print statement into two parts to avoid line truncation
              PRINT '(A, I3, A, F6.1, A, F6.1)', "prsurv for timepoint ", j, ": ", real(prsurv(leftCases(testi), j))
          END DO
      END DO

      end if ! end of print_check
    END IF ! end of isPhase2RE
    
    ! ============================================
    ! Below is for Phase 1 or Phase2CR
    ! ============================================
    IF (isPhase1 .OR. isPhase2CR) THEN
      prl = pr(leftCases,:) ! status change no matter what the status is (censoring, failure, cause specific failure, etc)
      prr = pr(rightCases,:)
      IF (isPhase2CR) THEN
        prl_m = prl ! the exact same for CR
        prr_m = prr ! should be the same for CR
      END IF
    END IF

    ! ============================================
    ! Below is for Phase2RE
    ! ============================================
    IF (isPhase2RE) THEN
      !PRINT *, "LINE 1399: prl, prr, prl_m, prr_m"
      prl = prsurv(leftCases,:) ! survival times ! prsurv_sub
      prr = prsurv(rightCases,:) 
      prl_m = pr(leftCases,:) ! recurrent event times
      prr_m = pr(rightCases,:)! pr_sub
      ! LATER TO DO: JTH CASE do-loop: eventsRight; eventsLeft
      ! prl_loop = prsurv(leftCases_loop,:)
      ! prr_loop = prsurv(rightCases_loop,:)
      ! prl_m_loop = pr(leftCases_loop,:)
      ! prr_m_loop = pr(rightCases_loop,:)
      ! splitLeft_loop = size(leftCases_loop, dim = 1) to give how many records there are
      ! eventsLeft_loop = sum(prl_loop * &
      ! & spread(dSorted(1:splitLeft_loop), 2, nt_death), DIM = 1)
      ! eventsRight_loop = sum(prr_loop * &
      ! & spread(dSorted(splitLeft_loop+1:nCases), 2, nt_death), DIM = 1)
    END IF

      IF (print_check) THEN
          PRINT *, "pr2"
          ! Print the last 5 elements of the 'pr' array in a readable way
          PRINT *, "pr2 array for the last 5 elements of leftCases:"
          DO testi = MAX(1, SIZE(leftCases) - 4), SIZE(leftCases)  ! Limit to the last 5 elements or the size of leftCases
              PRINT *, "Person ID: ", leftCases(testi), " with delta = ", REAL(dSorted_m(testi)), ":"
              PRINT '(F6.1)', (pr(testi, j), j = 1, SIZE(pr2, 2))  ! Print each element with 1 decimal
          END DO

          PRINT *, "leftCases array with size", SIZE(leftCases), ":", leftCases
          PRINT *, "splitLeft:", splitLeft

          ! Debugging or notification message
          PRINT *, "HHHHHHHHEEEEEEEEEEEEEEEEEEEEYYYYYYYYYYYYYY LOOOOOOOOOOOOOK @ MEEEEEEEE"
          PRINT *, "HHHHHHHHEEEEEEEEEEEEEEEEEEEEYYYYYYYYYYYYYY LOOOOOOOOOOOOOK @ MEEEEEEEE"
          PRINT *, "HHHHHHHHEEEEEEEEEEEEEEEEEEEEYYYYYYYYYYYYYY LOOOOOOOOOOOOOK @ MEEEEEEEE"
          PRINT *, "HHHHHHHHEEEEEEEEEEEEEEEEEEEEYYYYYYYYYYYYYY LOOOOOOOOOOOOOK @ MEEEEEEEE"

          PRINT *  ! Blank line for readability

          ! Debug information for personID_og if needed
          ! PRINT *, "Use below to double check against atrisk/events using checking_atrisk_events.R"
          ! PRINT *, "left:", SIZE(personID_og(1:splitLeft))
          ! PRINT *, personID_og(1:splitLeft)
          ! PRINT *, "right:", SIZE(personID_og(splitLeft + 1:nCases))
          ! PRINT *, personID_og(splitLeft + 1:nCases)

          PRINT *, "size(prsurv_sub)", SIZE(prsurv_sub, 1), " x ", SIZE(prsurv_sub, 2)
          PRINT *, "size(pr2surv_sub)", SIZE(pr2surv_sub, 1), " x ", SIZE(pr2surv_sub, 2)
          PRINT *, "size(pr_sub)", SIZE(pr_sub, 1), " x ", SIZE(pr_sub, 2)
          PRINT *, "size(pr2_sub)", SIZE(pr2_sub, 1), " x ", SIZE(pr2_sub, 2)
          PRINT *
          PRINT *, "prl_m array for all records:"
          DO testi = 1, SIZE(leftCases)  ! Iterate through all records
              PRINT *, "Person ", leftCases(testi), " with delta = ", REAL(dSorted_m(testi)), ":"
              DO j = 1, SIZE(prl_m, 2)  ! Print all timepoints for each person
                  PRINT '(A, I3, A, F6.1)', "  Timepoint ", j, ": ", REAL(prl_m(testi, j))
              END DO
              PRINT *  ! New line for better readability
          END DO
      END IF


      if (splitLeft <= 0) then
        PRINT *, "Warning: splitLeft is 0 or negative"
        STOP
      end if

    ! ============================================
    ! survival events (death)
    ! ============================================
    !IF (isPhase2RE) PRINT *, "LINE 1467: survival events"
    eventsLeft = sum(prl * &
        & spread(dSorted(1:(splitLeft-1)), 2, nt_death), DIM = 1)
    eventsRight = sum(prr * &
        & spread(dSorted(splitLeft:nCases), 2, nt_death), DIM = 1)

    ! ============================================
    ! endpoint events (CR or RE)
    ! ============================================
    IF (.NOT. isPhase1) THEN
      ! endpoint events (overall for CR; RE for RE)
      ! use dSorted_m, which is endpoint delta (CR: cause m, RE: recurrent event)
      eventsLeft_m = sum(prl_m * &
                    & spread(dSorted_m(1:(splitLeft-1)), 2, nt), DIM = 1)
      eventsRight_m = sum(prr_m * &
                    & spread(dSorted_m(splitLeft:nCases), 2, nt), DIM = 1)
    END IF

    IF (print_check) THEN
        ! Print the array in a visually organized format
        PRINT *, "Array dSorted_m (1:(splitLeft-1)):"
        PRINT *, "splitLeft:", splitLeft

        ! Loop through each person up to splitLeft-1 and display relevant details
        DO testi = 1, splitLeft - 1
            PRINT *, "Person", testi, "dSorted_m:", dSorted_m(testi)
            PRINT *, "Status change over time points (prl_m):"
            DO j = 1, nt  ! Loop over time points
                PRINT '(F8.4)', prl_m(testi, j)
            END DO
            PRINT *  ! Newline for readability
        END DO

        ! Display summary information
        PRINT *, "dSorted_m(1:(splitLeft-1)):", dSorted_m(1:(splitLeft-1))
        PRINT *, "spread(dSorted_m(1:(splitLeft-1)), 2, nt):"
        PRINT *, spread(dSorted_m(1:(splitLeft-1)), 2, nt)

        ! Compute and display transformed arrays and results
        PRINT *, "prl_m * spread(dSorted_m(1:(splitLeft-1)), 2, nt):"
        PRINT *, prl_m * spread(dSorted_m(1:(splitLeft-1)), 2, nt)

        PRINT *, "Sum of eventsLeft_m along dimension 1:"
        PRINT *, sum(prl_m * spread(dSorted_m(1:(splitLeft-1)), 2, nt), DIM = 1)

        ! Display dimensions and content of prl_m
        PRINT *, "prl_m with dimensions rows =", SIZE(prl_m, 1), ", cols =", SIZE(prl_m, 2)
        PRINT *, prl_m
        PRINT *
    END IF

    ! ============================================
    ! Below is for Phase 2 (Endpoint: CR)
    ! ============================================
    IF (isPhase2CR) THEN
      ! Initialize the group vector with zeros
      group_cr = 0
      ! Assign 1 to the indices in leftCases
      group_cr(leftCases) = 1
      ! Assign 2 to the indices in rightCases
      group_cr(rightCases) = 2
    END IF

    ! ============================================
    ! Below is for Phase 1 or Phase2CR
    ! ============================================
    ! looking at "at risk" set now (we want overall failure for both Step1 and Step2CR)
    IF (isPhase1 .OR. isPhase2CR) THEN
      pd1 = sum(prl, DIM = 1) ! for group 1
      pd2 = sum(prr, DIM = 1) ! for group 2
      ! at risk is same for Phase1 and Phase2CR because same subjects and same timepoint
      atRiskLeft(1) = splitLeft - 1
      atRiskRight(1) = nSubj - splitLeft + 1

      if (nt .NE. nt_death) then
        PRINT *, "Phase1/2CR Error: nt (", nt, ") is not equal to nt_death (", nt_death, "). Program will stop."
        if (isPhase2CR) THEN
          PRINT *, "Since CR points are subset of overall survival, we just use overall survival unique failure time points."
        END IF
        STOP
      end if

      !PRINT *, "tt3"
      !IF (isPhase1) THEN
      !  atRiskLeft(1) = splitLeft - 1
      !  atRiskRight(1) = nCases - splitLeft + 1
      !ELSE IF (isPhase2CR) THEN
      !  atRiskLeft(1) = splitLeft_m - 1
      !  atRiskRight(1) = nCases - splitLeft_m + 1
      !END IF
      !PRINT *, "tt4"

      DO j = 2, nt
        atRiskLeft(j) = atRiskLeft(j-1) - pd1(j-1)
        atRiskRight(j) = atRiskRight(j-1) - pd2(j-1)
      END DO
      !PRINT *, "tt5"
      
    ! ============================================
    ! Below is for Phase 2 (Endpoint: RE)
    ! ============================================
    ELSE IF (isPhase2RE) THEN
      !PRINT *, "LINE 1569: pr2l, pr2r for atRisk_m"
      ! because this is witihin phase2re, there are multiple records per person
      ! we need to identify 
      ! 1) at risk at survival times (using colsums pr2surv)
      ! 2) at risk at recurrent event times (using colsums of pr2)
      ! Earlier, we already defined:
      ! 3) death events at survival times (using prsurv and delta)
      ! 4) recurrent events at RE times (using pr and delta_m)

      ! pr2 is status change for recurrent event
      ! pr2surv is status change for terminal event
      ! only important for RE
      ! LATER TO DO: JTH CASE do-loop at-risk:
      ! leftCases_loop = leftCases + add each jth case
      ! rightCases_loop = rightCases - rm jth case 
      ! then: pr2l_loop = pr2(leftCases_loop,:) 
      ! pr2r_loop = pr2(rightCases_loop,:) 
      ! pd1_loop = sum(pr2l_loop, DIM = 1)
      ! pd2_loop = sum(pr2r_loop, DIM = 1)
      ! atRiskLeft_loop = pd1_loop
      pr2l = pr2(leftCases,:) 
      pr2r = pr2(rightCases,:)
      pd1_m = sum(pr2l, DIM = 1) ! RE at risk for group 1
      pd2_m = sum(pr2r, DIM = 1) ! RE at risk for group 2
      atRiskLeft_m = pd1_m
      atRiskRight_m = pd2_m
      !PRINT *, "-------------------"
      !print *, 'Shape of pr2:', shape(pr2)
      !print *, 'Shape of pr2l:', shape(pr2l) 
      !print *, 'Shape of pr2r:', shape(pr2r)     
      !PRINT *, "ROWS pr2:", size(pr2,1)
      !PRINT *, "COLS pr2:", size(pr2,2)
      !PRINT *, "ROWS pr2l:", size(pr2l,1)
      !PRINT *, "COLS pr2l:", size(pr2l,2)
      !PRINT *, "ROWS pr2r:", size(pr2r,1)
      !PRINT *, "COLS pr2r:", size(pr2r,2)
      !PRINT *, "LENGTH pd1_m:", size(pd1_m,1)
      !PRINT *, "LENGTH pd2_m:", size(pd2_m,1)
      !PRINT *, "-------------------"
      !PRINT *, "Row 1 of pr2 (1 person with 23 timepoints): "
      !PRINT *, pr2(1, :)
      !PRINT *, "Col 1 of pr2 (55 people at first timepoint): "
      !PRINT *, pr2(:, 1)
      !PRINT *, "-------------------"
      !PRINT *, "Row 1 of pr2surv (1 person with 13 timepoints): "
      !PRINT *, pr2surv(1, :)
      !PRINT *, "Col 1 of pr2surv (55 people at first timepoint): "
      !PRINT *, pr2surv(:, 1)
      !PRINT *, "-------------------"

      pr2survl = pr2surv(leftCases,:)
      pr2survr = pr2surv(rightCases,:)
      pd1 = sum(pr2survl, DIM = 1) ! death at risk in RE setting for group 1
      pd2 = sum(pr2survr, DIM = 1) ! death at risk in RE setting for group 2
      atRiskLeft = pd1
      atRiskRight = pd2
      !print *, 'Shape of pr2surv:', shape(pr2surv)
      !print *, 'Size of pr2surv:', size(pr2surv)
      !PRINT *, "ROWS pr2surv:", size(pr2surv,1)
      !PRINT *, "COLS pr2surv:", size(pr2surv,2)
      !PRINT *, "ROWS pr2survl:", size(pr2survl,1)
      !PRINT *, "COLS pr2survl:", size(pr2survl,2)
      !PRINT *, "ROWS pr2survr:", size(pr2survr,1)
      !PRINT *, "COLS pr2survr:", size(pr2survr,2)
      !PRINT *, "LENGTH pd1:", size(pd1,1)
      !PRINT *, "LENGTH pd2:", size(pd2,1)
      !PRINT *, "-------------------"

      !PRINT *, "---- LEFT NODE: survival ----"
      !! Check the shape of pr2survl
      !print *, 'Shape of pr2survl:', shape(pr2survl)
      !print *, 'Size of pr2survl:', size(pr2survl)
      !PRINT *, "ROWS pr2survl:", size(pr2survl,1)
      !PRINT *, "COLS pr2survl:", size(pr2survl,2)
      !PRINT *, "Death event: at risk for group 1"
      !PRINT *, "LENGTH (vector) pd1 = sum(pr2survl, dim = 1):", size(pd1,1)
      !PRINT *, "---- RIGHT NODE: survival ----"
      !PRINT *, "ROWS pr2survr:", size(pr2survr,1)
      !PRINT *, "COLS pr2survr:", size(pr2survr,2)
      !PRINT *, "Death event: at risk for group 1"
      !PRINT *, "LENGTH (vector) pd2 = sum(pr2survr, dim = 1):", size(pd2,1)
      !PRINT *, "-------------------"
      !PRINT *, "ROWS pr2l:", size(pr2l,1)
      !PRINT *, "COLS pr2l:", size(pr2l,2)
      !PRINT *, "Recurrent event: at risk for group 1"
      !PRINT *, "LENGTH (vector) pd1_m = sum(pr2l, dim = 1):", size(pd1_m,1)
      !PRINT *, "-------------------"
      !PRINT *, "end of testing"

      !TESTINGpd1 = 0.0_dp  ! Initialize to zero
      !PRINT *, "size(leftCases)", size(leftCases)
      !do iiii = 1, 21
      !    TESTINGpd1 = TESTINGpd1 + pr2survl(iiii, :)
      !end do
      !print *, 'Manual sum result (TESTINGpd1) with length ', size(TESTINGpd1)
      !PRINT *, TESTINGpd1
      !PRINT *, "fortran pd1 with length ", size(pd1)
      !PRINT *, pd1
      
      !PRINT *, "leftCases"
      !PRINT *, leftCases
      !PRINT *, "rightCases"
      !PRINT *, rightCases
      !PRINT *, "nt_death", nt_death
      !PRINT *, "nt", nt
      !PRINT *, "--------------------------"
      !PRINT *, "RE: at risk left"
      !PRINT *, atRiskLeft_m
      !PRINT *, "Death in RE: at risk left"
      !PRINT *, atRiskLeft
      !PRINT *, "--------------------------"
      !PRINT *, "splitLeft"
      !PRINT *, splitLeft
      !PRINT *, "splitLeftFinal"
      !PRINT *, splitLeftFinal

    END IF

    ! ============================================
    ! Below is for all Phases
    ! ============================================
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
    
    ! ============================================
    ! Below is for Phase 2 (Endpoint: RE)
    ! ============================================
    IF (isPhase2RE) THEN
      ! FIRST: WE WANT TO OBTAIN RE START/STOP INDICES FOR EACH PERSON
      ! Initialize variables
      uniqueCount = 0
      
      ! Check if the arrays are already allocated and deallocate them if they are
      if (allocated(firstIndex)) then
          deallocate(firstIndex)
      end if

      if (allocated(lastIndex)) then
          deallocate(lastIndex)
      end if
      ALLOCATE(firstIndex(MAXVAL(personID_new)))  ! Allocate maximum size
      ALLOCATE(lastIndex(MAXVAL(personID_new)))   ! Allocate maximum size
      ! Initialize firstIndex and lastIndex to -1 for all potential person IDs
      firstIndex = -1
      lastIndex = -1
      ! Find unique IDs and their indices
      DO doi = 1, nCases
          IF (firstIndex(personID_new(doi)) == -1) THEN
              uniqueCount = uniqueCount + 1
              firstIndex(personID_new(doi)) = doi
          END IF
          lastIndex(personID_new(doi)) = doi
      END DO

      if (allocated(uniqueID)) then
          deallocate(uniqueID)
      end if
      ! Allocate arrays for unique IDs based on the number of unique IDs found
      ALLOCATE(uniqueID(uniqueCount))
      ! Fill the uniqueID array
      doj = 1
      DO doi = 1, nCases
          IF (firstIndex(personID_new(doi)) /= -1) THEN
              ! Check if this personID is already recorded
              IF (doj > 1 .AND. personID_new(doi) == uniqueID(doj - 1)) CYCLE
              uniqueID(doj) = personID_new(doi)
              doj = doj + 1
          END IF
      END DO
      ! Now set the first and last indices for the new arrays
      DO doi = 1, uniqueCount
          firstIndex(doi) = firstIndex(uniqueID(doi))
          lastIndex(doi) = lastIndex(uniqueID(doi))
      END DO

      IF (print_check) THEN
        PRINT *, "Unique person IDs and their corresponding first and last indices:"
        DO doi = 1, uniqueCount
            PRINT *, "Person ", uniqueID(doi), ":::: First Index: ", &
            & firstIndex(doi), ", Last Index: ", lastIndex(doi)
        END DO

        PRINT *, "Summary of splitLeft and splitLeftFinal:"
        PRINT *, "splitLeft:       ", splitLeft
        PRINT *, "splitLeftFinal:  ", splitLeftFinal
        PRINT *
        PRINT *, "PersonID (New):"
        PRINT *, "----------------"
        PRINT *, "personID_new:               ", personID_new
        PRINT *, "personID_new(splitLeft):     ", personID_new(splitLeft)
        PRINT *, "personID_new(splitLeftFinal):", personID_new(splitLeftFinal)
        PRINT *
        PRINT *, "PersonID (Original):"
        PRINT *, "---------------------"
        PRINT *, "personID_og:               ", personID_og
        PRINT *, "personID_og(splitLeft):     ", personID_og(splitLeft)
        PRINT *, "personID_og(splitLeftFinal):", personID_og(splitLeftFinal)
        PRINT *
        PRINT *, "RecordID Details:"
        PRINT *, "------------------"
        PRINT *, "recordID with size:", SIZE(recordID)
        PRINT *, recordID
        PRINT *
      END IF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ============================================
      ! within person do-loop (doi) for Phase 2 RE
      ! ============================================
      DO doi = personID_new(splitLeft), personID_new(splitLeftFinal)
        ! Get start and end indices for the current person
        dostart = firstIndex(doi)
        doend = lastIndex(doi)
        ! Determine the number of cases for the current person
        numCases = doend - dostart + 1

        IF (print_check) THEN

          PRINT *, "------------------------------------------------------------"
          PRINT *, "------------------------------------------------------------"
          PRINT *, "------------------------------------------------------------"
          PRINT *, "Person Summary:"
          PRINT *, "This person (ID: ", doi, ") has ", numCases, " records."
          PRINT *, "Starting at: ", dostart, " and ending at: ", doend
          PRINT *, "Between ID = ", personID_new(splitLeft), " and ID = ", personID_new(splitLeftFinal)
          !PRINT *, " ----------- From personID_new, we have doi: ", doi, "---------------"
          PRINT *
          PRINT *, "Person Summary:"
          PRINT *, "---------------"
          PRINT *, "This person (ID: ", doi, ") has ", numCases, " records."
          PRINT *, "Starting at: ", dostart, " and ending at: ", doend
          PRINT *

          PRINT *, "Case Information:"
          PRINT *, "------------------"
          PRINT *, "Total number of cases: ", nCases
          PRINT *

          PRINT *, "Left Cases Details:"
          PRINT *, "--------------------"
          PRINT *, "Size of leftCases: ", SIZE(leftCases)
          PRINT *, leftCases
          PRINT *

          PRINT *, "Right Cases Details:"
          PRINT *, "---------------------"
          PRINT *, "Size of rightCases: ", SIZE(rightCases)
          PRINT *, rightCases
        END IF
        
        all_delta_same = .true.  ! Initialize to true

        ! Process records for the current person
        DO doj = dostart, doend
            ! Check if current element is different from the first element
            IF (xSorted(doj) /= xSorted(dostart)) THEN
                all_delta_same = .false.  ! Set to false if any difference is found
                PRINT *, "ERROR: Not all records are the same."
                PRINT "(A, I5, A, F6.2)", "doj: xSorted for record", doj, ": ", xSorted(doj)
                PRINT "(A, I5, A, F6.2)", "dstart: xSorted for record", dostart, ": ", xSorted(dostart)
                PRINT "(A, I5, A, F6.2)", "dend: xSorted for record", doend, ": ", xSorted(doend)
                PRINT *, "personID_new", size(personID_new)
                PRINT *, personID_new
                PRINT *, "recordID", size(recordID)
                PRINT *, recordID
                PRINT *, "xSorted", size(xSorted)
                PRINT *, xSorted
                PRINT *, "stopping"
                PRINT *, "ERROR: Not all records are the same."
                STOP
            END IF
        END DO

        ! If all elements are the same, print a success message
        !IF (all_delta_same) THEN
            !PRINT *, "Great! They are all the same."
        !END IF

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!! leftCases_loop and eventsRight_loop !!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! leftCases_loop = leftCases + add each jth case
        ! rightCases_loop = rightCases - rm jth case 
        leftCases_loop = recordID(1:(doend)) !recordID(1:(doend))
        leftPeople_loop = personID_new(1:(doend))
        leftPeople_loop_og = personID_og(1:(doend))
        CALL find_unique(leftPeople_loop, unique_leftPeople_loop, nleftPeople_loop)

        rightCases_loop = recordID(doend+1:nCases) !recordID(doend+1:nCases)
        rightPeople_loop = personID_new(doend+1:nCases)
        rightPeople_loop_og = personID_og(doend+1:nCases)
        CALL find_unique(rightPeople_loop, unique_rightPeople_loop, nrightPeople_loop)

        IF (print_check) THEN

        if (size(recordID) .LT. 242) THEN
          PRINT *, "recordID with size", size(recordID)
          PRINT *, recordID
          PRINT *, "personID_new with size", size(personID_new)
          PRINT *, personID_new
          PRINT *, "doend+1 = ", doend+1, " and nCases = ", nCases
          PRINT *, "testing recordID size:", size(recordID)
        END IF
          PRINT *
          PRINT *, "leftCases_loop with size", size(leftCases_loop)
          PRINT *, leftCases_loop
          PRINT *, "leftPeople_loop with size", size(leftPeople_loop)
          PRINT *, leftPeople_loop
          PRINT *, "nleftPeople_loop:", nleftPeople_loop
          PRINT *, "leftPeople_loop_og with size:", size(leftPeople_loop_og)
          PRINT *, leftPeople_loop_og
          PRINT *
          PRINT *, "rightCases_loop with size", size(rightCases_loop)
          PRINT *, rightCases_loop
          PRINT *, "rightPeople_loop with size", size(rightPeople_loop)
          PRINT *, rightPeople_loop
          PRINT *, "rightPeople_loop_og with size:", size(rightPeople_loop_og)
          PRINT *, rightPeople_loop_og
        END IF
         

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!! eventsLeft_loop and eventsRight_loop !!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !PRINT *, " starting LOOP events and at risk for left and right node"
        ! LATER TO DO: JTH CASE do-loop: eventsRight; eventsLeft
        splitLeft_loop = size(leftCases_loop, dim = 1) ! to give how many records there are

        ! left node 
        ! terminal events
        prl_loop = prsurv(leftCases_loop,:) ! death event times
        eventsLeft_loop = sum(prl_loop * &
        & spread(dSorted(1:splitLeft_loop), 2, nt_death), DIM = 1)
        ! recurrent events
        prl_m_loop = pr(leftCases_loop,:) ! recurrent event times
        eventsLeft_m_loop = sum(prl_m_loop * &
        & spread(dSorted_m(1:splitLeft_loop), 2, nt), DIM = 1)

        ! right node
        ! terminal events
        prr_loop = prsurv(rightCases_loop,:) ! prsurv_sub
        eventsRight_loop = sum(prr_loop * &
        & spread(dSorted(splitLeft_loop+1:nCases), 2, nt_death), DIM = 1)
        ! this is where we are having some issues^
        
        ! recurrent events
        prr_m_loop = pr(rightCases_loop,:) !pr_sub
        eventsRight_m_loop = sum(prr_m_loop * &
        & spread(dSorted_m(splitLeft_loop+1:nCases), 2, nt), DIM = 1)

        IF (print_check) THEN
          ! Print the entire prl_m_loop array in a readable way
          PRINT *, "prr_m_loop array:"
          DO testi = 1, SIZE(prr_m_loop, 1)  ! Loop over rows
              PRINT *, "Person ID: ", rightPeople_loop_og(testi), ": "
              PRINT '(F6.1)', (prr_m_loop(testi, j), j = 1, SIZE(prr_m_loop, 2))  ! Print each element with 1 decimal
          END DO

          PRINT *, "prr_loop array:"
          DO testi = 1, SIZE(prr_loop, 1)  ! Loop over rows
              PRINT *, "Person ID: ", rightPeople_loop_og(testi), ": "
              PRINT '(F6.1)', (prr_loop(testi, j), j = 1, SIZE(prr_loop, 2))  ! Print each element with 1 decimal
          END DO

          PRINT *, "RECURRENT EVENTS: eventsLeft_m_loop array:"
          DO testi = 1, SIZE(eventsLeft_m_loop)  ! Loop over elements in the 1D array
              PRINT '(A, I3, A, F6.1)', "Element(", testi, ") = ", eventsLeft_m_loop(testi)
          END DO

          PRINT *, "prl_loop array:"
          DO testi = 1, SIZE(prr_loop, 1)  ! Loop over rows
              PRINT *, "Person ID: ", rightPeople_loop_og(testi), ": "
              PRINT '(F6.1)', (prr_loop(testi, j), j = 1, SIZE(prr_loop, 2))  ! Print each element with 1 decimal
          END DO

          PRINT *, "---"
          PRINT *, "DEATH EVENTS: eventsRight_loop array:"
          DO testi = 1, SIZE(eventsRight_loop)  ! Loop over elements in the 1D array
              PRINT '(A, I3, A, F6.1)', "Element(", testi, ") = ", eventsRight_loop(testi)
          END DO
          
        PRINT *, "pr2 array:"
        DO testi = 1, SIZE(pr2(:,leftCases_loop), 2)  ! Loop over rows
            PRINT *, "id:", leftCases_loop(testi)
            PRINT '(F6.1)', (pr2(leftCases_loop(testi), j), j = 1, SIZE(pr2, 2))
        END DO

        END IF ! end of print_check

        !PRINT *, shape(prsurv_sub)
        !PRINT *, shape(prsurv)
        !PRINT *, shape(pr2surv_sub)
        !PRINT *, shape(pr2surv)
        !PRINT *, size(dSorted)
        !PRINT *, "prsurv_sub for rightCases_loop across all time points for RE"
        !PRINT *, prsurv_sub(rightCases_loop,:) ! rows are people
        !PRINT *
        !PRINT *, "prsurv_sub for rightCases_loop at 104th timepoint for RE"
        !PRINT *, prsurv_sub(rightCases_loop,104) ! column is timepoint
        !PRINT *
        !PRINT *, "dSorted(splitLeft_loop+1:nCases)"
        !PRINT *, dSorted(splitLeft_loop+1:nCases)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!! atRiskLeft_loop and atRiskRight_loop !!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! LATER TO DO: JTH CASE do-loop at-risk:
        ! pr2
        pr2l_loop = pr2(leftCases_loop,:) !pr2_sub
        pr2r_loop = pr2(rightCases_loop,:) 
        pd1_m_loop = sum(pr2l_loop, DIM = 1)
        pd2_m_loop = sum(pr2r_loop, DIM = 1)
        atRiskLeft_m_loop = pd1_m_loop
        atRiskRight_m_loop = pd2_m_loop

        ! pr2surv
        pr2survl_loop = pr2surv(leftCases_loop,:)  ! pr2surv_sub
        pr2survr_loop = pr2surv(rightCases_loop,:) 
        pd1_loop = sum(pr2survl_loop, DIM = 1)
        pd2_loop = sum(pr2survr_loop, DIM = 1)
        atRiskLeft_loop = pd1_loop
        atRiskRight_loop = pd2_loop

        !PRINT *, "leftCases_Loop with size:", size(leftCases_loop)
        !PRINT *, leftCases_loop
        !PRINT *, "dim of pr2surv is:", SHAPE(pr2surv)
        !PRINT *, "dim of pr2surv_sub is:", SHAPE(pr2surv_sub)
        !PRINT *, "pd1_loop with size:", size(pd1_loop)
        !PRINT *, "pd2_loop with size:", size(pd2_loop)
        !PRINT *, "atRiskLeft_loop with size:", size(atRiskLeft_loop)
        !PRINT *, "------------------------------------"
        !PRINT *, "size(personID_new):", size(personID_new)
        !PRINT *, "size(personID_og):", size(personID_og)
        !PRINT *, "size(recordID):", size(recordID)
        !PRINT *, "size(pr2_sub)", size(pr2_sub,1), " x ", size(pr2_sub,2)
        !PRINT *, "size(pr2surv_sub)", size(pr2surv_sub,1), " x ", size(pr2surv_sub,2)
        !PRINT *, "size(pr_sub)", size(pr_sub,1), " x ", size(pr_sub,2)
        !PRINT *, "size(prsurv_sub)", size(prsurv_sub,1), " x ", size(prsurv_sub,2)
        !PRINT *
        !PRINT *, "size(pr2)", size(pr2,1), " x ", size(pr2,2)
        !PRINT *, "size(pr2surv)", size(pr2surv,1), " x ", size(pr2surv,2)
        !PRINT *, "------------------------------------"
        !PRINT *, "------------------------------------"
        !PRINT *, "~~"

        !!!!!!! =============
        !PRINT "(A, F6.2, A, F6.2)", "xSorted(doend) is: ", xSorted(doend), " and (doend+1) - 1d-8 is: ", xSorted(doend + 1) - 1d-8

        ! if the case is not the last case with this covariate value, cycle
        IF (xSorted(doend) .GE. (xSorted(doend+1) - 1d-8)) CYCLE
 
        !PRINT *, "rule is ", rule
        IF (rule == 5) THEN
          IF (PRINT_CHECK) PRINT *, "~~~~~ SPLITTING TEST: PHASE 2 (RE): Q_LR (extension of gray's for RE) ~~~~~"
          ! atRisk is same for RE and survival within RE test statistic Q_LR
          ! IF (PRINT_CHECK) THEN
          !     PRINT *, "starting RE_INFO"
          !     PRINT *, "leftCases_loop"
          !     PRINT *, leftCases_loop
          !     PRINT *, "group 1: left node"
          !     CALL MeanFreqFunc(nt, nt_death, surv_tp, end_tp, &
          !     atRiskLeft_m_loop, eventsLeft_m_loop, atRiskLeft_loop, eventsLeft_loop, &
          !     survRE_left, dRhat_left, mu_left, dmu_left)

          !     PRINT *, "group 2: right node"
          !     CALL MeanFreqFunc(nt, nt_death, surv_tp, end_tp, &
          !     atRiskRight_m_loop, eventsRight_m_loop, atRiskRight_loop, eventsRight_loop, &
          !     survRE_right, dRhat_right, mu_right, dmu_right)

          !     PRINT *, "======================"

          !     ! Derivative of equation (2.2) in Ghosh and Lin (2000) using product rule
          !     ! Y_bar(t) = sum_i=1^n Y_i(t) = atriskRE vector
          !     ! dMi_hat(i) = eventsRE(i) - atriskRE(i)*dRhat(i)
          !     ! dMi_hat^D(i) = eventsSurv(i) - atriskRE(i) * (eventsSurv(i)/atriskSurv(i))

          !     ! dNi(t)
          !     dNi = prl_m_loop * &
          !     & spread(dSorted_m(1:splitLeft_loop), 2, nt)
          !     ! Y_i(t) * dRhat(t)
          !     YidR = pr2l_loop * &
          !     & spread(dRhat_left, 1, size(leftCases_loop))
          !     ! dMihat
          !     dMi = dNi - YidR

          !     ! dNi^D(t)
          !     dNiD = prl_m_loop * &
          !     & spread(dSorted(1:splitLeft_loop), 2, nt)

          !     PRINT *, "atRiskLeft_m_loop with size ", size(atRiskLeft_m_loop)
          !     PRINT *, atRiskLeft_m_loop

          !     !dLambdahat^D(t)
          !     tmp_events1 = 0.0_dp                          ! Set all elements of tmp_events1 to zero
          !     dlam = 0.0_dp                                  ! Initialize dlam to zero
          !     ! Calculate tmp_events1 as sum over the first dimension (summing across rows)
          !     tmp_events1 = sum(prl_m_loop * &
          !         & spread(dSorted(1:splitLeft_loop), 2, nt), DIM=1)
          !     ! Loop to calculate dlam with a check for zero in atRiskLeft_m_loop
          !     DO tmp_i = 1, nt
          !       PRINT *, "timepoint:", tmp_i
          !       PRINT *, "at risk at this timepoint:"
          !       PRINT *, atRiskLeft_m_loop(tmp_i)
          !         IF (atRiskLeft_m_loop(tmp_i) /= 0.0_dp) THEN
          !             dlam(tmp_i) = tmp_events1(tmp_i) / atRiskLeft_m_loop(tmp_i)
          !         ELSE
          !             dlam(tmp_i) = 0.0_dp 
          !         END IF
          !     END DO
          !     ! Y_i(t) * dLambdahat^D(t)
          !     YidLam = pr2l_loop * &
          !     & spread(dlam, 1, size(leftCases_loop))
          !     ! dMiD
          !     dMiD = dNiD - YidLam

          !     CALL dPsi_indiv(size(leftCases_loop), nleftPeople_loop, nt, nt_death, &
          !             survRE_left, dmu_left, mu_left, dMi, dMiD, &
          !             atRiskLeft_m_loop, dPsi_left)
          !     PRINT *, "end of call dPsi_indiv"
          ! END IF

          !PRINT *, "calling re_info for LEFT!!!"
          !PRINT *, "leftCases_loop with ", size(leftCases_loop), " total records"
          !PRINT *, leftCases_loop
          !PRINT *, "leftPeople_loop_og with size", size(leftPeople_loop_og)
          !PRINT *, leftPeople_loop_og
          !PRINT *, "leftPeople_loop with size", nleftPeople_loop
          !PRINT *, leftPeople_loop
          !CALL PrintSurvival("Survival (Left):", nt_death, surv_tp, &
          !& atRiskLeft_loop, eventsLeft_loop)
          !CALL PrintSurvival("Endpoint (Left):", nt, end_tp, &
          !& atRiskLeft_m_loop, eventsLeft_m_loop)
          
          CALL RE_INFO(nt, nt_death, surv_tp, end_tp, &
            & size(leftCases_loop), leftCases_loop, dSorted_m(1:splitLeft_loop), &
            & dSorted(1:splitLeft_loop), nleftPeople_loop, &
            & atRiskLeft_m_loop, eventsLeft_m_loop, atRiskLeft_loop, eventsLeft_loop, &
            & prl_m_loop, pr2l_loop, &
            & dmu_left, dPsi_left)


          !PRINT *, "calling re_info for RIGHT!!!"
          !PRINT *, "rightCases_loop has ", size(rightCases_loop), " total records."
          !PRINT *, rightCases_loop
          !PRINT *, "rightPeople_loop_og with size", size(rightPeople_loop_og)
          !PRINT *, rightPeople_loop_og
          !PRINT *, "rightPeople_loop with size", nrightPeople_loop
          !PRINT *, rightPeople_loop
          !CALL PrintSurvival("Survival (Right):", nt_death, surv_tp, &
          !& atRiskRight_loop, eventsRight_loop)
          !CALL PrintSurvival("Endpoint (Right):", nt, end_tp, &
          !& atRiskRight_m_loop, eventsRight_m_loop)

          CALL RE_INFO(nt, nt_death, surv_tp, end_tp, &
            & size(rightCases_loop), rightCases_loop, dSorted_m(splitLeft_loop+1:nCases), &
            & dSorted(splitLeft_loop+1:nCases), nrightPeople_loop, &
            & atRiskRight_m_loop, eventsRight_m_loop, atRiskRight_loop, eventsRight_loop, &
            & prr_m_loop, pr2r_loop, &
            & dmu_right, dPsi_right)

          IF (print_check) THEN
              PRINT *, "row    personID_og    10      42      64      71      78       88"
              PRINT '(A, F7.4, A, F7.4, A, F7.4, A, F7.4, A, F7.4, A, F7.4)', "row    personID_og    ", &
                        surv_tp(10), "   ", surv_tp(42), "   ", &
                        surv_tp(64), "   ", surv_tp(71), "   ", &
                        surv_tp(78), "    ", surv_tp(88)
              PRINT *, "-----------------------------------------------------"
              !DO tmp_i = 1, SIZE(rightCases_loop)
              !  caseR = rightCases_loop(tmp_i)  ! Use the integer index directly
              !  ! Adjusted print format for REAL values (F7.2 to display floating-point numbers)
              !  PRINT '(I4, 1X, I10, 3X, F7.2, 1X, F7.2, 1X, F7.2, 1X, F7.2, 1X, F7.2, 1X, F7.2)', &
              !    & tmp_i, rightPeople_loop_og(tmp_i), &
              !    & prsurv_sub(caseR, 10), prsurv_sub(caseR, 42), prsurv_sub(caseR, 64), &
              !    & prsurv_sub(caseR, 71), prsurv_sub(caseR, 78), prsurv_sub(caseR, 88)
              !END DO
              PRINT *, "-----------------------------------------------------"
              PRINT *
            END IF

          !PRINT *, "starting generalized weighted logrank (RE test)"
          CALL GeneralizedWeightedLR_RE(nt, nleftPeople_loop, nrightPeople_loop, &
              & atRiskLeft_m_loop, atRiskRight_m_loop, &
              & leftCases_loop, leftPeople_loop, &
              & rightCases_loop, rightPeople_loop, &
              & dmu_left, dmu_right, dPsi_left, dPsi_right, &
              & valuej_num, valuej_denom, valuej)
          !PRINT *, "end of generalized weighted logrank test: RE"
          !PRINT *, "test stat:", valuej

          IF (ieee_is_nan(valuej)) then
            ! Make sure there are no NaN test statistics
            PRINT *, "The test statistic is NaN. Stopping the execution."
            PRINT *, "Test statistic is: ", valuej
            PRINT *, "Num:", valuej_num
            PRINT *, "Denom:", valuej_denom
            PRINT *, "Variable number #", i, " amongst total ", nv, " variables."
            STOP
          !ELSE
            !PRINT "(A, F10.5)", "Q_LR test statistic valuej = ", valuej
          END IF 
        
        ELSE ! for rule == 5
          PRINT *, "ERROR: RULE SHOULD BE 5 BECAUSE PHASE 2 RE."
          STOP ! do not delete this stop error condition
        END IF ! end of RULE == 5

        if (print_check) PRINT *, "set:", set, " and valuej: ", valuej, " and maxValueXm:", maxValueXm

        IF ((set .EQ. 0) .OR. (valuej .GT. maxValueXm)) THEN
          IF (print_check) PRINT *, "valuej 1"
          ! if first value or value > current max, save
          IF (rUnifSet .EQ. 1) THEN
            cutoff = rUnif
          ELSE
            cutoff = (xSorted(doend) + xSorted(doend+1))/2.d0
          END IF
          maxValueXm = valuej
          tieValue = 1
          set = 1
        ELSE IF (valuej > (maxValueXm - 1d-8)) THEN
          ! if value is a tie, randomly determine if cutoff should be taken
          tieValue = tieValue + 1
          IF (rnd(0.d0, 1.d0) < (1.d0 / REAL(tieValue))) THEN
            cutoff = (xSorted(doend) + xSorted(doend+1))/2.d0
          END IF
        END IF
      ! ============================================
      END DO ! END OF doi do-loop for Phase 2 RE
      ! ============================================
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END IF ! END OF ISPHASE2RE

    ! ============================================
    ! Below is for Phase 1 or Phase2CR
    ! ============================================
    IF (isPhase1 .OR. isPhase2CR) THEN
      cnt = 1
      cnt_m = 1
      splitLeft_doloop = splitLeft
      splitLeftFinal_doloop = splitLeftFinal
      
      !PRINT *, "========== STARTING DO LOOP FOR Jth CASE =========="
      DO j = splitLeft_doloop, splitLeftFinal_doloop
        
        IF (isPhase1 .OR. isPhase2CR) THEN
          ! at risk indicators for jth case
          ! number of events for jth case
          pd1 = prr(cnt,:)
          !write(*,'(F6.2)', advance='no') prr(cnt, :)
          cnt = cnt + 1
          D = pd1*delta(cases(j))
          Rcum(1) = 0.0
          Rcum(2) = pd1(1)
                IF (isPhase2CR) THEN
                  pd1_m = prr_m(cnt_m,:)
                  cnt_m = cnt_m + 1
                  D_m = pd1_m*delta_m(cases(j))
                  Rcum_m(1) = 0.0
                  Rcum_m(2) = pd1_m(1)
                END IF
          DO k = 3, nt
            IF (pd1(k-1) .GT. 1d-8) THEN
              Rcum(k) = Rcum(k-1) + pd1(k-1)
            ELSE 
              Rcum(k) = Rcum(k-1)
            END IF
                IF (isPhase2CR) THEN
                  IF (pd1_m(k-1) .GT. 1d-8) THEN
                    Rcum_m(k) = Rcum_m(k-1) + pd1_m(k-1)
                  ELSE 
                    Rcum_m(k) = Rcum_m(k-1)
                  END IF
                END IF
          END DO
          Rcum = 1.d0 - Rcum 
          ! number at risk
          ! add the jth case to the left node
          atRiskLeft = atRiskLeft + Rcum
          ! remove the jth case from the right node
          atRiskRight = atRiskRight - Rcum
          ! number of events
          !PRINT *, "BEFORE removing jth case: eventsRight"
          !PRINT *, eventsRight
          ! add the jth case to the left node
          eventsLeft = eventsLeft + D
          ! remove the jth case from the right node
          eventsRight = eventsRight - D
          IF (isPhase2CR) THEN
            Rcum_m = 1.d0 - Rcum_m
            ! number at risk
            ! add the jth case to the left node
            atRiskLeft_m = atRiskLeft_m + Rcum_m
            ! remove the jth case from the right node
            atRiskRight_m = atRiskRight_m - Rcum_m
            ! number of events
            ! add the jth case to the left node
            eventsLeft_m = eventsLeft_m + D_m
            ! remove the jth case from the right node
            eventsRight_m = eventsRight_m - D_m
            !PRINT *, "atRiskLeft_m"
            !PRINT *, atRiskLeft_m
            !PRINT *, "atRiskLeft"
            !PRINT *, atRiskLeft
            !PRINT *, "this should be the same because at risk for Survival and same at risk for priority cause."
          END IF
        END IF

        ! if the case is not the last case with this covariate value, cycle
        IF (xSorted(j) .GE. (xSorted(j+1) - 1d-8)) CYCLE

        if (rule .NE. 1 .AND. rule .NE. 2 .AND. rule .NE. 3 .AND. rule .NE. 4) then
            PRINT *, "ERROR: RULE SHOULD BE 1, 2, 3, OR 4 BECAUSE EITHER PHASE 1 SURVIVAL OR PHASE 2 CR."
            STOP
        end if

        ! calculate test statistic
        IF (rule == 1) THEN
          !PRINT *, "~~~~~ SPLITTING TEST: PHASE 1: logrank test ~~~~~"
          CALL logrank(atRiskLeft, atRiskRight, eventsLeft, numJ, &
                    & denJ, valuej)
          !PRINT *, "logrank test statistic valuej = ", valuej
        ELSE IF (rule == 2) THEN
          !PRINT *, "~~~~~ SPLITTING TEST: PHASE 1: truncated mean test ~~~~~"
          ! Call to the meanSplit subroutine
          CALL meanSplit(atRiskLeft, atRiskRight, eventsLeft, eventsRight, valuej)
          ! PRINT *, "mean test statistic valuej = ", valuej
        ELSE IF (rule == 3) THEN
          !PRINT *, "~~~~~ SPLITTING TEST: PHASE 2 (CR): gray's test ~~~~~"
          ! set up for crstm
          IF (isPhase2CR) THEN
            !!!!!!!PRINT *, "CR: Gray's Test Set-Up to split nodes"
            ! below are placeholder/initial values
            rho_set0 = 0;
            ist_set1 = 1
            sorted_cases = cases
            sorted_cases_index = (/(i,i=1,size(sorted_cases))/)
            !!!!!!!PRINT *, "cases"
            !!!!!!!PRINT *, cases
            !!!PRINT *, "old sorted_cases"
            !!!PRINT *, sorted_cases
            !if (isPhase2CR) CALL qsort4(sorted_cases, cases, 1, size(cases))
            if (isPhase2CR) CALL qsort4(sorted_cases, sorted_cases_index, 1, size(sorted_cases))

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
          PRINT *, "THIS SHOULD NEVER BE TRIGGERED, SINCE RULE == 5 ONLY WORKS FOR ISPHASE2RE"
          STOP
        END IF

        ! Make sure there are no NaN test statistics
        IF (ieee_is_nan(valuej)) then
          PRINT *, "The test statistic is NaN. Stopping the execution."
          PRINT *, valuej
          STOP
        END IF

        IF ((set .EQ. 0) .OR. (valuej .GT. maxValueXm)) THEN
          ! if first value or value > current max, save
          IF (rUnifSet .EQ. 1) THEN
            !PRINT *, "rUnifSet: ", rUnifSet
            !PRINT *, "rUnif: ", rUnif
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
    ! ============================================
    END DO ! END OF JTH CASE DO-LOOP
    ! ============================================
    END IF ! end of phase1/phase2cr

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! if not successful, cycle to next covariate
    ! this condition should never be true
    IF (set .EQ. 0) CYCLE

    ! if successful, determine if it yields the maximum value of the
    ! covariates considered
    IF ((splitVar .EQ. -1) .OR. (maxValueXm .GT. maxValueSplit)) THEN
      !PRINT *, "SuCCESSful."
      !PRINT *, "splitVar: ", splitVar, "and kv = ", kv

      ! if first non-zero or largest value, keep cutoff and value and
      ! reset tie counter to 1
      splitVar = kv
      !PRINT *, "splitVar: ", splitVar
      maxValueSplit = maxValueXm
      tieCovariate = 1

      ! count the number of cases in the left node
      lft = count(xSorted .LE. cutoff)
      IF (isPhase1 .OR. isPhase2CR) THEN
        casesOut = cases
      ELSE IF (isPhase2RE) THEN
      !PRINT *, "personID_new with size:", size(personID_new)
      !PRINT *, personID_new
      !PRINT *, "recordID with size:", size(recordID)
      !PRINT *, recordID
      !PRINT *, "personID with size:", size(personID)
      !PRINT *, personID
        casesOut = personID 
        casesOutRE = recordID
        casesOut_people_RE = personID_og
        casesOut_record_RE = recordID_og
        !PRINT *, "casesOut with size:", size(casesOut)
        !PRINT *, casesOut
        !PRINT *, "casesOutRE with size", size(casesOutRE)
        !PRINT *, casesOutRE  
        !PRINT *, "casesOut_people_RE with size", size(casesOut_people_RE)
        !PRINT *, casesOut_people_RE 
        END IF

      
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
        IF (isPhase1 .OR. isPhase2CR) THEN
          casesOut = cases
        ELSE IF (isPhase2RE) THEN
          casesOut = personID
          casesOutRE = recordID
          casesOut_people_RE = personID_og
          casesOut_record_RE = recordID_og
        END IF
        IF (print_check) THEN
          PRINT *, "casesOut at end with size:", size(casesOut)
          PRINT *, casesOut
        END IF

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

  END DO ! end of covariate do-loop nv

  !PRINT *, "splitVar: ", splitVar
  ! if no split was possible return
  if (splitVar .EQ. -1) RETURN

  ! if successful at finding a split set flag and return
  splitFound = 1
  
  IF (print_check) THEN
    PRINT *, "---------- POST -----------"
    PRINT *, "splitFound: ", splitFound
    PRINT *, "tfindsplit casesIn with size ", size(casesIn)
    PRINT *, casesIn
    PRINT *, "tfindsplit casesOut with size ", size(casesOut)
    PRINT *, casesOut
    IF (isPhase2RE) THEN 
      PRINT *, "tfindsplit casesOutRE with size ", size(casesOutRE)
      PRINT *, casesOutRE
    END IF
    PRINT *, "lft:", lft
    PRINT *, "splitVar: ", splitVar
    PRINT *, "cutoffBest: ", cutoffBest
    PRINT "(A, F6.2)", "cutoff:", cutoff
    PRINT *, "nCuts: ", nCuts
    PRINT *, "---------------------"
  END IF

  !IF (isPhase2RE) PRINT *, "============================END OF tfindSplit============================"
  
  RETURN

END SUBROUTINE tfindSplit
! =================================================================================

! Calculate the Kaplan Meier estimator
! ns integer, the number of time points
! nj real(:), at risk
! oj real(:), events
! z real(:), estimator
SUBROUTINE kaplan(ns, nj, oj, z)
  use, intrinsic :: ieee_arithmetic
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ns
  REAL(dp), DIMENSION(1:ns), INTENT(IN) :: nj
  REAL(dp), DIMENSION(1:ns), INTENT(IN) :: oj
  REAL(dp), DIMENSION(1:ns), INTENT(OUT) :: z

  INTEGER :: i

  REAL(dp), DIMENSION(1:ns) :: num

  !PRINT *, "ns:", ns
  !PRINT *, "nj:", nj
  !PRINT *, "oj:", oj
  !PRINT *, "******************** kaplan ********************"
  num = nj - oj

  IF (nj(1) > 1d-8) THEN
    z(1) = num(1) / nj(1)
  ELSE
    !PRINT *, "nj(1) = ", nj(1)
    !PRINT *
    !PRINT *, "nj:", nj
    !PRINT *
    PRINT *, "Error: nj is 0 at time point 1 so stopping (aka time 0)"
    STOP
  END IF

  IF (ns .LT. 2) RETURN

  DO i = 2, ns
    IF (nj(i) > 1d-8) THEN
      z(i) = (num(i)/nj(i)) * z(i-1)
    ELSE
      !IF (isPhase2RE) PRINT *, "z(",i-1,") = ", z(i-1)
      z(i) = z(i-1)
    END IF
    !IF (isPhase2RE) PRINT *, "z(",i,") = ", z(i)
  END DO

  DO i = 1, ns
    IF (z(i) < 0.0 .OR. z(i) > 1.0 .OR. ieee_is_nan(z(i))) THEN
        PRINT *, "CHECK Kaplan Meier estimate!!"
        PRINT *, "z(i) wrong at time point ", i, ": ", z(i)
        PRINT *, "num(i): ", num(i)
        PRINT *, "nj(i): ", nj(i)
        PRINT *, "oj(i): ", oj(i)
        !PRINT *, "delta:", delta
        !PRINT *, "delta_m:", delta_m
        PRINT *, "at risk total:", nj
        PRINT *, "events total:", oj
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
! =================================================================================

SUBROUTINE PrintSurvival(label, count, timePoint, AT_RISK, EVENTS)
  IMPLICIT NONE

  ! Inputs
  CHARACTER(LEN=*), INTENT(IN) :: label
  INTEGER, INTENT(IN) :: count
  REAL(8), INTENT(IN) :: timePoint(:), AT_RISK(:), EVENTS(:)

  ! Local variables
  INTEGER :: tmp_i

  ! Print the header
  PRINT *, label, count
  PRINT *, "-----------------------------------------"
  PRINT *, "index    timePoint    AT_RISK    EVENTS"
  PRINT *, "-----------------------------------------"

  ! Loop through data and print
  DO tmp_i = 1, SIZE(AT_RISK)
    PRINT '(I4, 1X, F20.18, 1X, I4, 1X, I4)', tmp_i, timePoint(tmp_i), INT(AT_RISK(tmp_i)), INT(EVENTS(tmp_i))
  END DO
  PRINT *, "-----------------------------------------"
  PRINT *
END SUBROUTINE PrintSurvival
! =================================================================================

!CALL RE_INFO(nt, nt_death, surv_tp, end_tp, &
!          & size(groupCases_loop), leftCases_loop, &
!          & dSorted_m(1:splitLeft_loop), dSorted(1:splitLeft_loop), nleftPeople_loop, &
!          & atRiskLeft_m_loop, eventsLeft_m_loop, atRiskLeft_loop, eventsLeft_loop, &
!          & prl_m_loop, pr2l_loop, &
!          & dmu_left, dPsi_left)

SUBROUTINE RE_INFO(nt_endpoint, nt_survival, tp_survival, tp_endpoint, &
                    nrecords, groupCases_loop, &
                    dSorted_m_loop, dSorted_loop, ngroupPeople_loop, &
                    Nj_endpoint, Oj_endpoint, Nj_survival, Oj_survival, &
                    prgroup_m_loop, pr2group_loop, &
                    dmu_group, dPsi_group)

  INTEGER, INTENT(IN) :: nt_endpoint
  INTEGER, INTENT(IN) :: nt_survival
  REAL(dp), DIMENSION(1:nt_survival), INTENT(IN) :: tp_survival
  REAL(dp), DIMENSION(1:nt_endpoint), INTENT(IN) :: tp_endpoint
  INTEGER, INTENT(IN) :: nrecords
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: groupCases_loop
  INTEGER, DIMENSION(:), INTENT(IN) :: dSorted_loop
  INTEGER, DIMENSION(:), INTENT(IN) :: dSorted_m_loop
  INTEGER, INTENT(IN) :: ngroupPeople_loop
  REAL(dp), DIMENSION(1:nt_endpoint), INTENT(IN) :: Nj_endpoint ! at risk RE
  REAL(dp), DIMENSION(1:nt_endpoint), INTENT(IN) :: Oj_endpoint ! RE
  REAL(dp), DIMENSION(1:nt_survival), INTENT(IN) :: Nj_survival ! at risk death
  REAL(dp), DIMENSION(1:nt_survival), INTENT(IN) :: Oj_survival ! death events
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: prgroup_m_loop
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: pr2group_loop
  REAL(dp), DIMENSION(1:nt_endpoint), INTENT(OUT) :: dmu_group
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: dPsi_group  ! Output matrix

  REAL(dp), DIMENSION(nt_endpoint) :: tmp_events1, tmp_events2, dlam 
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: dNi, YidR, dMi, dNiD, YidLam, dMiD
  REAL(dp), DIMENSION(1:nt_endpoint) :: survRE_group
  REAL(dp), DIMENSION(1:nt_endpoint) :: dRhat_group
  REAL(dp), DIMENSION(1:nt_endpoint) :: mu_group

  INTEGER :: tmp_i, i

! use below to check, can use checking_atrisk_events.R in R
!!! Check At Risk and Events Inputs
!PRINT *
!PRINT *, nt_survival, " timepoints:"
!PRINT *, "------------------------------------------"
!PRINT '(A)', "TimePoint        AtRisk        TerminalEvents"
!PRINT *, "------------------------------------------"
!!! Loop through each survival time point
!DO i = 1, nt_survival
!    PRINT '(F10.5, 1X, I12, 1X, I12)', tp_survival(i), INT(Nj_survival(i)), INT(Oj_survival(i))
!END DO
!PRINT *, "------------------------------------------"
!PRINT *
!PRINT *, nt_endpoint, " timepoints:"
!PRINT *, "------------------------------------------"
!PRINT '(A)', "TimePoint        AtRisk        RecurrentEvents"
!PRINT *, "------------------------------------------"
!!! Loop through each RE time point
!DO i = 1, nt_endpoint
!    PRINT '(F10.5, 1X, I12, 1X, I12)', tp_endpoint(i), INT(Nj_endpoint(i)), INT(Oj_endpoint(i))
!END DO
!PRINT *, "------------------------------------------"


! FIRST: get GROUP INFORMATION
!PRINT *, "calling meanfreqfunc from RE_info"
CALL MeanFreqFunc(nt_endpoint, nt_survival, tp_survival, tp_endpoint, &
Nj_endpoint, Oj_endpoint, Nj_survival, Oj_survival, &
survRE_group, dRhat_group, mu_group, dmu_group)

! Derivative of equation (2.2) in Ghosh and Lin (2000) using product rule
! dNi(t)
dNi = prgroup_m_loop * &
& spread(dSorted_m_loop, 2, nt_endpoint) ! we want RE dSorted_m
! Y_i(t) * dRhat(t)
YidR = pr2group_loop * &
& spread(dRhat_group, 1, nrecords)
! dMihat
dMi = dNi - YidR
! dNi^D(t)
dNiD = prgroup_m_loop * &
& spread(dSorted_loop, 2, nt_endpoint) ! we want survival dSorted
!dLambdahat^D(t)
tmp_events1 = 0.0_dp                          ! Set all elements of tmp_events1 to zero
dlam = 0.0_dp                                  ! Initialize dlam to zero
! Calculate tmp_events1 as sum over the first dimension (summing across rows)
tmp_events1 = sum(prgroup_m_loop * &
    & spread(dSorted_loop, 2, nt_endpoint), DIM=1) ! we want survival dSorted
! Loop to calculate dlam with a check for zero in Nj_endpoint
DO tmp_i = 1, nt_endpoint
    !PRINT *, "timepoint:", tmp_i
    !PRINT *, "at risk at this timepoint:"
    !PRINT *, Nj_endpoint(tmp_i)
    IF (Nj_endpoint(tmp_i) /= 0.0_dp) THEN
        dlam(tmp_i) = tmp_events1(tmp_i) / Nj_endpoint(tmp_i)
    ELSE
        dlam(tmp_i) = 0.0_dp 
    END IF
END DO
! Y_i(t) * dLambdahat^D(t)
YidLam = pr2group_loop * &
& spread(dlam, 1, nrecords)
! dMiD
dMiD = dNiD - YidLam

!PRINT *, "now obtaining dPsi_group"
CALL dPsi_indiv(nrecords, ngroupPeople_loop, nt_endpoint, nt_survival, &
        survRE_group, dmu_group, mu_group, dMi, dMiD, &
        Nj_endpoint, dPsi_group)
!PRINT *, "end of call dPsi_indiv"
!PRINT *, "dPsi_group with dim:", size(dPsi_group,1), size(dPsi_group,2)
!PRINT *, dPsi_group
!PRINT *, "ngroupPeople_loop:", ngroupPeople_loop, " and nrecords:", nrecords

END SUBROUTINE RE_INFO



! =========================================================

SUBROUTINE MeanFreqFunc(nt_endpoint, nt_survival, tp_survival, tp_endpoint, &
Nj_endpoint, Oj_endpoint, Nj_survival, Oj_survival, &
survRE, dRhat, mu, dmu)
  use, intrinsic :: ieee_arithmetic
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nt_endpoint
  INTEGER, INTENT(IN) :: nt_survival
  REAL(dp), DIMENSION(1:nt_survival), INTENT(IN) :: tp_survival
  REAL(dp), DIMENSION(1:nt_endpoint), INTENT(In) :: tp_endpoint
  REAL(dp), DIMENSION(1:nt_endpoint), INTENT(IN) :: Nj_endpoint ! at risk RE
  REAL(dp), DIMENSION(1:nt_endpoint), INTENT(IN) :: Oj_endpoint ! RE
  REAL(dp), DIMENSION(1:nt_survival), INTENT(IN) :: Nj_survival ! at risk death
  REAL(dp), DIMENSION(1:nt_survival), INTENT(IN) :: Oj_survival ! death events
  REAL(dp), DIMENSION(1:nt_endpoint), INTENT(OUT) :: survRE
  REAL(dp), DIMENSION(1:nt_endpoint), INTENT(OUT) :: dRhat
  REAL(dp), DIMENSION(1:nt_endpoint), INTENT(OUT) :: mu
  REAL(dp), DIMENSION(1:nt_endpoint), INTENT(OUT) :: dmu

  INTEGER :: i, j

  REAL(dp), DIMENSION(1:nt_endpoint) :: muTmp ! dont want to include 0
  REAL(dp), DIMENSION(1:nt_survival) :: survD
  LOGICAL :: print_check

  print_check = .FALSE. 

  if (print_check) then
    PRINT *, "******************** mean frequency function ********************"
    PRINT *, 'nt_endpoint = ', nt_endpoint
      ! Print tp_endpoint array with index, rounded to 1 decimal point
    PRINT *, 'tp_endpoint: '
    DO i = 1, nt_endpoint
        PRINT '(A, I3, A, F6.1)', 'tp_endpoint(', i, ') = ', tp_endpoint(i)
    END DO
    PRINT *, 'nt_survival = ', nt_survival
    ! Print tp_survival array with index, rounded to 1 decimal point
    PRINT *, 'tp_survival: '
    DO i = 1, nt_survival
        PRINT '(A, I3, A, F6.4)', 'tp_survival(', i, ') = ', tp_survival(i)
    END DO
  END IF

if (print_check) then
    PRINT *, "Input Parameters"
    
    PRINT *, 'nt_endpoint = ', nt_endpoint
    ! Print tp_endpoint, Nj_endpoint, and Oj_endpoint arrays together
    PRINT *, 'tp_endpoint, Nj_endpoint, and Oj_endpoint:'
    DO i = 1, nt_endpoint
        PRINT '(A, I3, A, F6.1, A, F6.1, A, F6.1)', &
              'Index(', i, '): tp_endpoint = ', tp_endpoint(i), &
              ', Nj_endpoint = ', Nj_endpoint(i), &
              ', Oj_endpoint = ', Oj_endpoint(i)
    END DO

    PRINT *, 'nt_survival = ', nt_survival
    ! Print tp_survival, Nj_survival, and Oj_survival arrays together
    PRINT *, 'tp_survival, Nj_survival, and Oj_survival:'
    DO i = 1, nt_survival
        PRINT '(A, I3, A, F6.4, A, F6.1, A, F6.1)', &
              'Index(', i, '): tp_survival = ', tp_survival(i), &
              ', Nj_survival = ', Nj_survival(i), &
              ', Oj_survival = ', Oj_survival(i)
    END DO


    PRINT *, "Starting MeanFreqFunc now..."
end if

  mu = 0.d0
  dmu = 0.d0
  dRhat = 0.d0
  ! Initialize survRE to avoid NaN issues
  survRE = 0.0_dp ! make sure this is OK, otherwise delete
  !survD = 0.0_dp ! make sure this is OK, otherwise delete
  
  !PRINT *, "******************** obtain survival KM (deaths) ********************"
if (print_check) then
  IF (Nj_survival(1) .LT. 1) THEN
    PRINT *, 'nt_survival = ', nt_survival
    ! Print tp_survival and Nj_survival arrays together
    PRINT *, 'tp_survival and Nj_survival:'
    DO i = 1, nt_survival
        PRINT '(A, I3, A, F6.4, A, F6.1)', 'Index(', i, '): tp_survival = ', tp_survival(i), ', Nj_survival = ', Nj_survival(i)
    END DO
    ! Print tp_endpoint and Nj_endpoint arrays together
    PRINT *, 'tp_endpoint and Nj_endpoint:'
    DO i = 1, nt_endpoint
        PRINT '(A, I3, A, F6.4, A, F6.1)', 'Index(', i, '): tp_endpoint = ', tp_endpoint(i), ', Nj_endpoint = ', Nj_endpoint(i)
    END DO
  END IF
end if
  
  !PRINT *, "call kaplan within MeanFreqFunc subroutine - LINE 2710"
  call kaplan(nt_survival, Nj_survival, Oj_survival, survD)

  IF (print_check) THEN
    PRINT *, "survD with size ", size(survD)
    PRINT *, survD
    PRINT *, " at risk" 
    PRINT *, Nj_survival
    PRINT *, "events"
    PRINT *, Oj_survival
  DO i = 1, nt_survival
    IF (ieee_is_nan(survD(i))) then
      PRINT *, i
      PRINT *, "timepoint", tp_survival(i)
      PRINT *, "survD(",i,"): ", survD(i)
      PRINT *, "survD(",i-1,"): ", survD(i-1)
      PRINT *, "survD(",i+1,"): ", survD(i+1)
      PRINT *
      PRINT *, "survD with size ", size(survD)
      PRINT *, survD(i)
      PRINT *, Nj_survival(i)
      PRINT *, Oj_survival(i)
      PRINT *, Nj_survival
      STOP
    END IF
  END DO
  END IF
  
  DO i = 1, nt_endpoint ! Loop over all elements in tp_endpoint (or another array of length nt_endpoint)
    !PRINT *, "i: ", i

    !PRINT '(A, I3, A, F6.4, A, F6.1, A, F6.1)', "For i = ", i, ", timepoint is ", tp_endpoint(i), & 
    !  " with ", Nj_endpoint(i), " at risk and ", Oj_endpoint(i), " events."

    ! to get the time points that have individuals at risk
    ! if number at risk is equal to 0, then skip to next time point index (i)
    IF (Nj_endpoint(i) .LT. 1d-8) CYCLE
    !PRINT *, "******************** dRhat ********************"
    dRhat(i) = Oj_endpoint(i) / Nj_endpoint(i) ! dRhat at each timepoint
    IF (ieee_is_nan(dRhat(i))) then
      PRINT *, "NaN: dRhat(",i,") = ", dRhat(i)
    END IF
    !PRINT *, "dRhat(",i,") = ", dRhat(i)

    !PRINT *, "******************** Recurrent KM (using RE times) ********************"
    ! survRE is the survival curve at the recurrent event times
    ! Handle timepoints before tp_survival(1)
    IF (tp_endpoint(i) < tp_survival(1)) THEN
      survRE(i) = survD(1)  ! Assign the survival at the first timepoint
    END IF

    ! Handle timepoints after tp_survival(nt_survival)
    IF (tp_endpoint(i) >= tp_survival(nt_survival)) THEN
      survRE(i) = survD(nt_survival)  ! Assign the survival at the last timepoint
    END IF

    DO j = 1, nt_survival-1  ! Loop over elements in tp_survival, up to tp_survival-1 (need to check interval between nt-1 and nt)
      ! Check if tp_endpoint(i) falls within the interval [tp_survival(j), tp_survival(j+1)]
      IF (tp_survival(j) <= tp_endpoint(i) .AND. tp_endpoint(i) < tp_survival(j+1)) THEN
        survRE(i) = survD(j) ! Assign survD(j) to survRE(i) if condition is met
        EXIT                 ! Exit the inner loop if condition is met
      END IF
    END DO ! j do-loop

  END DO ! i do-loop

  !survRE(nt_endpoint) = survD(nt_survival)

  IF (print_check) THEN
  PRINT *
  PRINT *, "endpoint"
  DO i = 1, nt_endpoint
    PRINT *, i
    PRINT *, tp_endpoint(i)
    PRINT *, survRE(i)
  END DO
  PRINT *
  PRINT *, "survival"
  DO i = 1, nt_survival
    PRINT *, i
    PRINT *, tp_survival(i)
    PRINT *, survD(i)
  END DO
  PRINT *
  END IF

  IF (print_check) then
    PRINT *, "survD has size: ", size(survD) ! 121
    PRINT *, "survRE has size: ", size(survRE) ! 324
  DO i = 1, nt_endpoint
    IF (ieee_is_nan(survRE(i))) then
      PRINT *, "i:",i
      PRINT *, "timepoint:", j, " = ", tp_survival(j)
      PRINT *, tp_survival(j)
      PRINT *, tp_endpoint(i)
      PRINT *, tp_survival(j+1)
      PRINT *, "survRE at timepoint ", tp_endpoint(i), " for ", i, " is: ", survRE(i)
      PRINT *, "survRE at timepoint ", tp_endpoint(i-1), " for ", i-1, "is: ", survRE(i-1)
      PRINT *, "survRE at timepoint ", tp_endpoint(i+1), " for ", i+1, "is: ", survRE(i+1)
      PRINT *, "size survRE", size(survRE)
      PRINT *, j, ": survD(", tp_survival(j),") is:", survD(j)
      PRINT *, "previous survD(", tp_survival(j-1),") is:", survD(j-1)
      PRINT *, "next survD(", tp_survival(j+1),") is:", survD(j+1)
        STOP
    END IF
  END DO
  END IF

  mu(1) = 0
  dmu(1) = 0
  if (nt_endpoint .LT. 2) RETURN
  DO i = 2, nt_endpoint
    dmu(i) = survRE(i) * dRhat(i)
    mu(i) = mu(i-1) + dmu(i)
    !IF (ieee_is_nan(dmu(i))) then
    !  PRINT *, "survRE(",i,"): ", survRE(i)
    !  PRINT *, "dRhat(",i,"): ", dRhat(i)
    !  PRINT *, "dmu(",i,"): ", dmu(i)      
    !END IF
    !if (mu(i) .NE. mu(i-1)) THEN
    !  PRINT *, "mu(",i,") is not equal to mu(",i-1,")."
    !END IF
  END DO

  IF (print_check) THEN
    DO i = 1, nt_endpoint
      IF (ieee_is_nan(dmu(i))) then
        PRINT *, "timepoint: ",i
        PRINT *, "dmu(",i,"): ", dmu(i)
        !PRINT *, "mu(",i-1,"): ", mu(i-1)
        !PRINT *, "mu(",i,"): ", mu(i)
        PRINT *, "survRE(",i,"): ", survRE(i)
        !PRINT *, "dRhat(",i,"): ", dRhat(i)
        STOP
      END IF
    END DO
  END IF


  !PRINT *, "end of MeanFreqFunc"

END SUBROUTINE MeanFreqFunc
! =================================================================================

SUBROUTINE weightKLR(ns, n1, n2, atrisk1, atrisk2, output_weight)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ns ! # time points for unique observed RE times
  INTEGER, INTENT(IN) :: n1 ! sample size in group 1
  INTEGER, INTENT(IN) :: n2 ! sample size in group 2
  REAL(dp), DIMENSION(1:ns), INTENT(IN) :: atrisk1 ! # subj in group 1 at risk at time t
  REAL(dp), DIMENSION(1:ns), INTENT(IN) :: atrisk2 ! # subj in group 2 at risk at time t
  REAL(dp), DIMENSION(1:ns), INTENT(OUT) :: output_weight

  INTEGER :: n 
  INTEGER :: i                       ! loop variable if needed

  n = n1 + n2
!  PRINT *, "atrisk1"
!  PRINT *, atrisk1
!  PRINT *, "atrisk2"
!  PRINT *, atrisk2
!  PRINT *, "n1=",n1, " and n2=", n2, " and n = n1+n2=", n

  DO i = 1, ns
    IF (atrisk1(i) + atrisk2(i) > 0.0_dp) THEN
!      PRINT *, "atrisk1(i) * atrisk2(i)"
!      PRINT *, atrisk1(i) * atrisk2(i)
!      PRINT *, "atrisk1(i) + atrisk2(i)"
!      PRINT *, atrisk1(i) + atrisk2(i)
!      PRINT *, "atrisk1(i) * atrisk2(i) / (atrisk1(i) + atrisk2(i))"
!      PRINT *, atrisk1(i) * atrisk2(i) / (atrisk1(i) + atrisk2(i))
!      PRINT *, "n / (n1 * n2)"
!      PRINT *, REAL(n, dp) / (REAL(n1, dp) * REAL(n2, dp))  ! Casting to REAL to avoid integer division
      ! Update the calculation for output_weight
      output_weight(i) = (atrisk1(i) * atrisk2(i) / (atrisk1(i) + atrisk2(i))) * (REAL(n, dp) / (REAL(n1, dp) * REAL(n2, dp)))
    ELSE
      !PRINT *, "set 0"
      output_weight(i) = 0.0_dp      ! Handle cases with no risk in either group
    END IF
!    PRINT *, "weight for timepoint ", i, " is: ", output_weight(i)
  END DO


END SUBROUTINE weightKLR


! ==================================================
SUBROUTINE CalculateREDenominator(K_LR, dPsi, n_people, n_records, people_loop, outer_sum)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n_people, n_records
    REAL(dp), DIMENSION(:), INTENT(IN) :: K_LR
    REAL(dp), DIMENSION(:,:), INTENT(IN) :: dPsi
    INTEGER, DIMENSION(:), INTENT(IN) :: people_loop
    REAL(dp), INTENT(OUT) :: outer_sum

    REAL(dp), DIMENSION(1:n_records) :: inner_integral
    REAL(dp), DIMENSION(1:n_people) :: inner_sum
    REAL(dp) :: denom_sum
    INTEGER :: person_ind, record_ind, i, j, n_tp, chunk_size
    LOGICAL :: print_check
    INTEGER, DIMENSION(n_records) :: new_people_loop

    print_check = .FALSE.

    !PRINT *, "Starting CalculateREDenominator Subroutine"
    !PRINT *, "n_people:", n_people
    !PRINT *, "n_records:", n_records
    !PRINT *, "shape of K_LR", shape(K_LR)
    !PRINT *, "shape of dPsi", shape(dPsi)

    ! Compute the inner integral for each record
    inner_integral = SUM(SPREAD(K_LR, 1, n_records) * dPsi, 2)
    ! after a certain point, dPsi should be equal to 0 for individuals
    ! but we don't calculate it this way in dPsi_indiv so need to account for that here


    !PRINT *, "shape of dPsi for person 1 across all 496 timepoints with size", shape(dPsi(1,:))
    !PRINT *, dPsi(1,:)
    !PRINT *, "Shape of dPsi for person 1 across all 496 timepoints: ", SHAPE(dPsi(1,:))
 
    n_tp = SIZE(dPsi, 2)       ! Number of timepoints (second dimension size)
    chunk_size = 10            ! Number of timepoints per row for better alignment
    
    !PRINT *, "Timepoints:"
    !DO i = 1, n_tp, chunk_size
    !  PRINT *
      ! Print the timepoint indices (adjust format width for better spacing)
    !  PRINT "(10I8)", (j, j = i, MIN(i+chunk_size-1, n_tp))
      
      ! Print the corresponding dPsi values
    !  PRINT "(10F8.2)", dPsi(1, i:MIN(i+chunk_size-1, n_tp))
      
      ! Print the corresponding K_LR values
    !  PRINT "(10F8.2)", K_LR(i:MIN(i+chunk_size-1, n_tp))

    !END DO


    !DO i = 1, n_tp, chunk_size
    !PRINT *
    !  WRITE(*, "(A5, 10I8)") "TP:", (j, j = i, MIN(i+chunk_size-1, n_tp))
    !  WRITE(*, "(A5, 10F8.5)") "dPSI:", dPsi(1, i:MIN(i+chunk_size-1, n_tp))
    !  WRITE(*, "(A5, 10F8.5)") "K_LR:", K_LR(i:MIN(i+chunk_size-1, n_tp))
    !  WRITE(*, "(A5, 10F8.5)") "K*dP:", (K_LR(i:MIN(i+chunk_size-1, n_tp)) * dPsi(1, i:MIN(i+chunk_size-1, n_tp)))
    !END DO

    !PRINT *, "inner_integral with size:", SHAPE(inner_integral)
    !PRINT *, inner_integral

    ! Initialize inner_sum
    inner_sum = 0.0_dp
    CALL create_new_vector(people_loop, size(people_loop), new_people_loop)
    !PRINT *, "people_loop:", people_loop
    !PRINT *, "new_people_loop:", new_people_loop
    !PRINT *, "people_loop(1):", people_loop(1)

    ! each componenet of sigma^2_LR
    ! Sum over records for each person
    DO person_ind = 1, n_people
        !PRINT *, "********** PERSON ", person_ind, " **********"
        denom_sum = 0.0_dp
        ! Sum contributions for records belonging to the current person
        DO record_ind = 1, n_records
          IF (new_people_loop(record_ind) == person_ind) THEN
            !PRINT *, "       Record #", record_ind
            !PRINT *, "              denom_sum:", denom_sum
            !PRINT *, "              inner_integral(record_ind) = ", inner_integral(record_ind)
            !PRINT *, "              denom_sum + inner_integral(record_ind) = ", denom_sum + inner_integral(record_ind)
            denom_sum = denom_sum + inner_integral(record_ind)
            !print *, "              record-loop denom_sum:", denom_sum
          !ELSE 
            !PRINT *, "*********************************************************************"
            !PRINT *, "********** ", people_loop(record_ind), " IS NOT EQUAL TO ", person_ind, " **********"
            !PRINT *, "*********************************************************************"
          END IF
        END DO
        !print *, "person-loop denom_sum:", denom_sum
        ! Store the square of the sum in inner_sum
        inner_sum(person_ind) = denom_sum**2
    END DO

    ! Calculate the outer sum
    outer_sum = SUM(inner_sum)

    if (print_check) then
    IF (outer_sum .EQ. 0) THEN
      !PRINT *, "dPsi"
      !PRINT *, dPsi
      !PRINT *, "K_LR"
      !PRINT *, K_LR
      PRINT *, "          inner_integral with size:", SIZE(inner_integral)
      PRINT *, inner_integral
      PRINT *, "          inner_sum with size: ", SIZE(inner_sum)
      PRINT *, inner_sum
      PRINT *, "          outer_sum: ", outer_sum
      PRINT *, "          n_people:", n_people
      PRINT *, "          n_records:", n_records
      PRINT *
    END IF
    end if

!PRINT *, "End of CalculateREDenominator"

END SUBROUTINE CalculateREDenominator



SUBROUTINE GeneralizedWeightedLR_RE(ns, n1, n2, atrisk1, atrisk2, &
                                    leftCases_loop, leftPeople_loop, &
                                    rightCases_loop, rightPeople_loop, &
                                    dmu1, dmu2, dPsi1, dPsi2, Q_LR, sigma2_LR, test_statistic)
  use, intrinsic :: ieee_arithmetic
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ns ! # time points for unique observed RE times
  INTEGER, INTENT(IN) :: n1 ! sample size in group 1
  INTEGER, INTENT(IN) :: n2 ! sample size in group 2
  REAL(dp), DIMENSION(1:ns), INTENT(IN) :: atrisk1 ! # subj in group 1 at risk at time t
  REAL(dp), DIMENSION(1:ns), INTENT(IN) :: atrisk2 ! # subj in group 2 at risk at time t
  INTEGER, DIMENSION(:) :: leftCases_loop
  INTEGER, DIMENSION(:) :: leftPeople_loop
  INTEGER, DIMENSION(:) :: rightCases_loop
  INTEGER, DIMENSION(:) :: rightPeople_loop 
  REAL(dp), DIMENSION(1:ns), INTENT(IN) :: dmu1 ! # dmu in group 1 
  REAL(dp), DIMENSION(1:ns), INTENT(IN) :: dmu2 ! # dmu in group 2
  REAL(dp), DIMENSION(1:size(leftCases_loop),1:ns), INTENT(IN) :: dPsi1
  REAL(dp), DIMENSION(1:size(rightCases_loop),1:ns), INTENT(IN) :: dPsi2
  REAL(dp), INTENT(OUT) :: Q_LR
  REAL(dp), INTENT(OUT) :: sigma2_LR
  REAL(dp), INTENT(OUT) :: test_statistic

  INTEGER :: nrecords1 ! sample size in group 1
  INTEGER :: nrecords2 ! sample size in group 2
  REAL(dp), DIMENSION(1:ns) :: K_LR ! K_LR weight
  REAL(dp), DIMENSION(1:ns):: Q_LR_vec

  REAL(dp), DIMENSION(1:size(leftCases_loop)) :: inner_integral1
  REAL(dp), DIMENSION(1:size(rightCases_loop)) :: inner_integral2
  REAL(dp), DIMENSION(1:n1) :: inner_sum1
  REAL(dp), DIMENSION(1:n2) :: inner_sum2
  INTEGER :: n1_ind, n2_ind, ns_ind, i, j, k
  REAL(dp) :: denom_sum, outer_sum1, outer_sum2

  nrecords1 = size(leftCases_loop)
  nrecords2 = size(rightCases_loop)

  !------------------------------------------------------
  !------------------ Calculate Weight ------------------
  !------------------------------------------------------
  CALL weightKLR(ns, n1, n2, atrisk1, atrisk2, K_LR) 
  !PRINT *, "ns: ", ns, "and n1: ", n1, " and n2: ", n2
  !PRINT *, "atrisk1"
  !PRINT *, atrisk1
  !PRINT *, "atrisk2"
  !PRINT *, atrisk2
  !PRINT *, "K_LR with size:", size(K_LR)
  !PRINT *, K_LR
  !PRINT *
  !------------------------------------------------------
  !---------------- Calculate Numerator -----------------
  !------------------------------------------------------
  ! numerator: Q_LR
  Q_LR_vec = K_LR * (dmu1- dmu2)
  Q_LR = SUM(Q_LR_vec)
  !PRINT *, "************** numerator Q_LR is: **************"
  !PRINT *, "SUM(K_LR * dmu1)"
  !PRINT *, SUM(K_LR * dmu1)
  !PRINT *, "SUM(K_LR * dmu2)"
  !PRINT *, SUM(K_LR * dmu2)
  !PRINT *, Q_LR
  !PRINT *, "K_LR"
  !PRINT *, K_LR
  !PRINT *, "dmu1"
  !PRINT *, dmu1
  !PRINT *, "dmu2"
  !PRINT *, dmu2
  !PRINT *, "dmu1-dmu2"
  !PRINT *, dmu1-dmu2
  !PRINT *, "numvec"
  !PRINT *, Q_LR_vec
  !PRINT *
  !------------------------------------------------------
  !---------------- Calculate Denominator ---------------
  !------------------------------------------------------
  ! denominator: variance sigma2_LR
  ! outer_sum1 + outer_sum2 is (3.3) in Ghosh and Lin (without the n fraction parts)
  ! together, they make up the variance sigma2_LR
  ! Call CalculateDenominator for group 1
  CALL CalculateREDenominator(K_LR, dPsi1, n1, nrecords1, leftPeople_loop, &
                              outer_sum1)
  !PRINT *, "denom outer_sum1:", outer_sum1
  !PRINT *, "dPsi1"
  !PRINT *, dPsi1
  ! Call CalculateDenominator for group 2
  CALL CalculateREDenominator(K_LR, dPsi2, n2, nrecords2, rightPeople_loop, &
                              outer_sum2)
  !PRINT *, "denom outer_sum2:", outer_sum2
  !PRINT *, "dPsi2"
  !PRINT *, dPsi2                        
  
  sigma2_LR = REAL(n2, dp) / (REAL(n, dp) * REAL(n1, dp)) * outer_sum1 + REAL(n1, dp) / (REAL(n, dp) * REAL(n2, dp)) * outer_sum2
  !PRINT *, "sigma2_LR denominator: ", sigma2_LR

  ! jan 10, 2025: updated to make test statistic 0 when numerator^2 is very small (essentially 0) and denominator is 0
  IF (Q_LR**2 <= 1.0E-10 .AND. sigma2_LR .EQ. 0.0) THEN
    ! Do something if numerator is very small and denominator is 0
    test_statistic = 0.0
  ELSE 
    test_statistic = (REAL(n1, dp) * REAL(n2, dp) / REAL(n, dp)) * (Q_LR**2)/(sigma2_LR )
  END IF
  !PRINT *, "test statistic z^2 = ", test_statistic
  !test_statistic = (REAL(n1, dp) * REAL(n2, dp) / REAL(n, dp)) * (Q_LR**2)/(sigma2_LR )

  !PRINT *, "test stat:", test_statistic
  !PRINT *, "numerator Q_LR^2:", Q_LR**2
  !PRINT *, "denominator:", sigma2_LR
  !PRINT *, "outer_sum1:", outer_sum1
  !PRINT *, "outer_sum2:", outer_sum2

  IF (sigma2_LR .EQ. 0.0 .AND. Q_LR**2 > 1.0E-10) THEN
    PRINT *, "BIG PROBLEM: DENOMINATOR IS 0 BUT NUMERATOR IS NOT -!!!!!!!!"
    PRINT *, "SUM(K_LR * dmu1)"
    PRINT *, SUM(K_LR * dmu1)
    PRINT *, "SUM(K_LR * dmu2)"
    PRINT *, SUM(K_LR * dmu2)
    PRINT *, "test stat:", test_statistic
    PRINT *, "numerator:", Q_LR**2
    PRINT *, "denominator:", sigma2_LR
    PRINT *, "outer_sum1:", outer_sum1
    PRINT *, "outer_sum2:", outer_sum2
    PRINT *, "stopping to debug"
    STOP
  END IF 

  IF (ieee_is_nan(test_statistic)) then
    !PRINT *, "dPsi1"
    !PRINT *, dPsi1
    !PRINT *, "dPsi2"
    !PRINT *, dPsi2
    PRINT *, "test stat:", test_statistic
    PRINT *, "numerator:", Q_LR**2
    PRINT *, "denominator:", sigma2_LR
    PRINT *, "outer_sum1:", outer_sum1
    PRINT *, "outer_sum2:", outer_sum2
    PRINT *, "stopping to test why test statistic is NaN."
    STOP
  END IF


END SUBROUTINE GeneralizedWeightedLR_RE


! =======================================

SUBROUTINE dPsi_indiv(nrecords, n, ns, ns_death, survRE, dmu, mu, dMi, dMiD, Ybar, dPsi_mat)
  IMPLICIT NONE

  ! Declare the types and intents of the input/output parameters
  INTEGER, INTENT(IN) :: nrecords
  !REAL(dp), INTENT(IN) :: personID ! personID
  !REAL(dp), INTENT(IN) :: recordID ! recordID
  INTEGER, INTENT(IN) :: n        ! Sample size (1:n for individuals)
  INTEGER, INTENT(IN) :: ns       ! Number of recurrent event time points (1:ns)
  INTEGER, INTENT(IN) :: ns_death ! Number of terminal event time points (1:ns_death)
  REAL(dp), DIMENSION(1:ns), INTENT(IN) :: survRE     ! Survival estimates at RE times
  REAL(dp), DIMENSION(1:ns), INTENT(IN) :: dmu
  REAL(dp), DIMENSION(1:ns), INTENT(IN) :: mu
  REAL(dp), DIMENSION(:,:), INTENT(IN) :: dMi
  REAL(dp), DIMENSION(:,:), INTENT(IN) :: dMiD
  REAL(dp), DIMENSION(1:ns), INTENT(IN) :: Ybar
  
  ! Variable Declarations
  INTEGER :: i, j, index_record
  REAL(dp), DIMENSION(1:nrecords, 1:ns) :: termA, termB, int_termB1, termB1

  REAL(dp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: dPsi_mat  ! Output matrix
  LOGICAL :: print_check

  ! Allocate the output matrix
  ALLOCATE(dPsi_mat(1:nrecords, 1:ns))

  ! Initialize the output matrix to zero
  dPsi_mat = 0.0_dp
  print_check = .FALSE.

  ! Compute TERM A for each time point
  DO i = 1, ns
    IF (Ybar(i) > 1d-8) THEN
      termA(:,i) = spread(survRE(i), 1, nrecords) * dMi(:,i) / (spread(Ybar(i), 1, nrecords) / n)
    ELSE
      termA(:,i) = 0.0_dp
    END IF
  END DO
  !below is old code that doesnt work because NaN produced
  !IF (ANY(Ybar /= 0.0_dp)) THEN
  !    termA = spread(survRE, 1, nrecords) * dMi / (spread(Ybar, 1, nrecords) / n)
  !ELSE
  !    termA = 0.0_dp
  !END IF


  ! spread(mu, 1, nrecords)
  DO i = 1, ns ! ns = recurrent event time points
    IF (Ybar(i) > 1d-8) THEN
      termB1(:,i) = dMiD(:,i) / (spread(Ybar(i), 1, nrecords) / n)
    ELSE
      termB1(:,i) = 0.0_dp
    END IF
  END DO
  !BELOW IS OLD CODE THAT DOESNT WORK B/C NAN produced
  !IF (ANY(Ybar /= 0.0_dp)) THEN
  !    PRINT *, "no Ybar is 0"
  !    termB1 = dMiD / (spread(Ybar, 1, nrecords) / n)
  !ELSE
  !    PRINT *, "at least one Ybar is 0"
  !    termB1 = 0.0_dp
  !END IF
  !IF (isPhase2RE .AND. size(termA,1) .EQ. 8) STOP

  ! first timepoint is equal to first column of B1, the first timepoint of B1
  ! initialize first timepoint
  int_termB1(:,1) = termB1(:,1)
  do i = 2, ns
    int_termB1(:,i) = int_termB1(:,i-1) + termB1(:,i)
  end do
  termB = spread(dmu, 1, nrecords) * int_termB1

  dPsi_mat = termA - termB ! changing to negative 1/10/25

if (print_check) then
  ! row is record, column is timepoint
  PRINT *, "size of survRE is: ", size(survRE)
  PRINT *, "size of Ybar is: ", size(Ybar)
  PRINT *, "size of dmu is: ", size(dmu)
  PRINT *, "size of mu is: ", size(mu)
  PRINT *, "row of dMi is: ", size(dMi,1), "and col is: ", size(dMi,2)
  PRINT *
  PRINT *, "size of spreading survRE by ", nrecords, " records is: ", &
          size(spread(survRE, 1, nrecords), 1), "and ", &
          size(spread(survRE, 1, nrecords), 2)
  PRINT *, "size of spreading Ybar by ", nrecords, " records is: ", &
          size(spread(Ybar, 1, nrecords), 1), "and ", &
          size(spread(Ybar, 1, nrecords), 2)
  PRINT *
end if

  IF (print_check) THEN
    PRINT *, "ybar with size:", size(ybar)
    PRINT *, Ybar
    PRINT *, "dmid with size:", size(dMiD)
    PRINT *, dMiD
    PRINT *, "denom"
    PRINT *, spread(Ybar, 1, nrecords) / n
    PRINT *

    DO j = 1, size(termA,1)
      PRINT *, "termA with dim: ", size(termA,1), " x ", size(termA,2)
      PRINT *, "row: ", j
      PRINT *, termA(j,:)
    END DO

    DO j = 1, size(termB1,1)
      PRINT *, "termB1 with dim: ", size(termB1,1), " x ", size(termB1,2)
      PRINT *, "row: ", j
      PRINT *, termB1(j,:)
    END DO

  END IF

if (print_check) then
  PRINT *, "int_termB1"
  PRINT *, int_termB1
  PRINT *, "spread(dmu, 1, nrecords)"
  PRINT *, spread(dmu, 1, nrecords)
  PRINT *
  PRINT *, "termB1 with dim: ", size(termB1,1), " x ", size(termB1,2)
  PRINT *, "int_termB1 with dim: ", size(int_termB1,1), " x ", size(int_termB1,2)
  PRINT *, "spread(dmu, 1, nrecords) with dim: ", size(spread(dmu, 1, nrecords),1), " x ", size(spread(dmu, 1, nrecords),2)
  PRINT *, "termB with dim: ", size(termB,1), " x ", size(termB,2)
  PRINT *, termB
  PRINT *
  PRINT *
  PRINT *
  PRINT *
  PRINT *
  !PRINT *, "termA + termB with size:", size(termA+termB,1), " x ", size(termA+termB,2)  
  !PRINT *, termA + termB
  !PRINT *
  !PRINT *
  !PRINT *
  !PRINT *
  !PRINT *
  PRINT *, "dPsi_mat with size:", size(dPsi_mat,1), " x ", size(dPsi_mat,2)  
  PRINT *, dPsi_mat
  PRINT *
  PRINT *
  PRINT *
  PRINT *
  PRINT *
end if 

END SUBROUTINE dPsi_indiv




!==========================================


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

  !PRINT *, "N1j:", N1j
  !PRINT *, "N2j:", N2j
  !PRINT *, "O1j:", O1j
  !PRINT *, "O2j:", O2j
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
        PRINT *, "STOPPING TO check why jumpcif is <0 or >1"
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

subroutine create_new_vector(personID_og, n, personID_new)
    implicit none
    integer, intent(in) :: n  ! Size of the original array
    integer, dimension(n), intent(in) :: personID_og  ! Original array
    integer, dimension(n), intent(out) :: personID_new  ! New array
    integer :: i, currentID

    currentID = 1  ! Start from 1 for new ID
    personID_new(1) = currentID  ! Assign the first ID
    !PRINT *, "n: ", n

    do i = 2, n
      !PRINT *, "i: ", i
      !PRINT *, personID_og(i)
      !PRINT *, personID_og(i-1)
        if (personID_og(i) /= personID_og(i - 1)) then
          !print *, "not equal"
            currentID = currentID + 1
        end if
        !print *, "personID_new at ", i, "is: ", currentID
        personID_new(i) = currentID
    end do
    !PRINT *, personID_new

end subroutine create_new_vector


!subroutine group_and_sort(vector1, vector2, vector3, vector4, vector5, vector6, n, combinedArray)
!    implicit none
!    INTEGER, INTENT(IN) :: n
!    REAL(dp), DIMENSION(1:n) :: vector1
!    INTEGER, DIMENSION(1:n) :: vector2, vector3, vector4, vector5, vector6
!    REAL(dp), DIMENSION(1:n, 6), INTENT(OUT) :: combinedArray  ! Five columns for combined array
!    INTEGER :: i, j

!    ! Combine vectors into a 2D array
!    combinedArray(:, 1) = vector1
!    combinedArray(:, 2) = REAL(vector2, dp)  ! Cast integers to real(dp)
!    combinedArray(:, 3) = REAL(vector3, dp)
!    combinedArray(:, 4) = REAL(vector4, dp)
!    combinedArray(:, 5) = REAL(vector5, dp)
!    combinedArray(:, 6) = REAL(vector6, dp)

!    print *, "COMBINED ARRAY BEFORE SORTING:"
!    PRINT *, "covar_sorted_RE, personID_og, recordID, personID, delta, delta_m"
!    do i = 1, 5 !n
!        print '(5X, F7.4, 6X, I6, 6X, I6, 6X, I6, 3X, I6, 3X, I6)', &
!              combinedArray(i, 1), INT(combinedArray(i, 2)), INT(combinedArray(i, 3)), &
!              INT(combinedArray(i, 4)), INT(combinedArray(i, 5)), INT(combinedArray(i, 6))
!    end do

!    ! Sort combinedArray ONLY based on vector1, vector2, and vector3
!    do i = 1, n-1
!        do j = i+1, n
!            if (combinedArray(i, 1) > combinedArray(j, 1)) then
!                call swap_rows(combinedArray, i, j, 6)
!            else if (combinedArray(i, 1) == combinedArray(j, 1)) then
!                if (combinedArray(i, 2) > combinedArray(j, 2)) then
!                    call swap_rows(combinedArray, i, j, 6)
!                else if (combinedArray(i, 2) == combinedArray(j, 2)) then
!                    if (combinedArray(i, 3) > combinedArray(j, 3)) then
!                        call swap_rows(combinedArray, i, j, 6)
!                    end if
!                end if
!            end if
!        end do
!    end do

!    ! Output the sorted array
!    print *, "Grouped and Sorted Data:"
!    PRINT *, "covar_sorted_RE, personID_og, recordID, personID, delta, delta_m"
!    do i = 1, 5 !n
!        print '(5X, F7.4, 6X, I6, 6X, I6, 6X, I6, 3X, I6, 3X, I6)', &
!              combinedArray(i, 1), INT(combinedArray(i, 2)), INT(combinedArray(i, 3)), &
!              INT(combinedArray(i, 4)), INT(combinedArray(i, 5)), INT(combinedArray(i, 6))
!    end do

!end subroutine group_and_sort

subroutine group_and_sort(inputVectors, n, combinedArray)
    implicit none
    INTEGER, INTENT(IN) :: n              ! Number of rows
    REAL(dp), DIMENSION(:, :), INTENT(IN) :: inputVectors ! Input vectors (2D array, variable columns)
    REAL(dp), DIMENSION(:, :), INTENT(OUT) :: combinedArray ! Combined array (same size as input)
    INTEGER :: numVectors                 ! Number of input vectors (columns in inputVectors)
    INTEGER :: i, j

    ! Determine the number of vectors (columns in inputVectors)
    numVectors = SIZE(inputVectors, DIM=2)

    ! Copy inputVectors to combinedArray
    combinedArray = inputVectors

    !print *, "COMBINED ARRAY BEFORE SORTING:"
    !do i = 1, MIN(5, n)
    !    write(*, '(1X, 7F10.4)') (combinedArray(i, j), j = 1, numVectors)
    !end do

    ! Sort combinedArray based on the first three columns
    do i = 1, n - 1
        do j = i + 1, n
            if (combinedArray(i, 1) > combinedArray(j, 1)) then
                call swap_rows(combinedArray, i, j, numVectors)
            else if (combinedArray(i, 1) == combinedArray(j, 1)) then
                if (combinedArray(i, 2) > combinedArray(j, 2)) then
                    call swap_rows(combinedArray, i, j, numVectors)
                else if (combinedArray(i, 2) == combinedArray(j, 2)) then
                    if (combinedArray(i, 3) > combinedArray(j, 3)) then
                        call swap_rows(combinedArray, i, j, numVectors)
                    end if
                end if
            end if
        end do
    end do

    !print *, "Grouped and Sorted Data:"
    !do i = 1, MIN(5, n)
    !    write(*, '(1X, 7F10.4)') (combinedArray(i, j), j = 1, numVectors)
    !end do
end subroutine group_and_sort


subroutine swap_rows(array, i, j, numCols)
    implicit none
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: i, j, numCols
    REAL(dp) :: temp
    INTEGER :: k

    ! Swap all columns between row i and row j
    do k = 1, numCols
        temp = array(i, k)
        array(i, k) = array(j, k)
        array(j, k) = temp
    end do
end subroutine swap_rows

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
use, intrinsic :: ieee_arithmetic
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nCases
  INTEGER, DIMENSION(1:nCases), INTENT(IN) :: casesIn
  REAL(dp), DIMENSION(1:nt), INTENT(OUT) :: Func
  REAL(dp), INTENT(OUT) :: mean

  INTEGER :: i, j, t, testj

  REAL(dp), DIMENSION(1:nt_death) :: Nj, Oj, Rb
  REAL(dp), DIMENSION(1:nt) :: Nj_m, Oj_m
  REAL(dp), DIMENSION(1:nt) :: jumpCIF
  REAL(dp), DIMENSION(1:nt) :: survRE
  REAL(dp), DIMENSION(1:nt) :: dRhat
!  REAL(dp), DIMENSION(1:nt) :: mu
  REAL(dp), DIMENSION(1:nt) :: dmu
  ! for Phase1/2CR: nt=nt_death
  ! for Phase2RE: nt != nt_death

  !IF (isPhase2RE) PRINT *, "******************** calcValueSingle ********************"
  !PRINT *, "estimate the survival/cif function and mean survival/cif time"
  Func = 0.d0
  mean = 0.d0

  IF (isPhase1 .OR. isPhase2CR) THEN
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
      !PRINT *, "mean survival time: ", mean
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
  END IF 

  IF (isPhase2RE) THEN

    ! Terminal Event Nj and Oj
    ! We now calculate the number of at risk cases at each time point using pr2 {nt}
    Nj_m = sum(pr2(casesIn,:), DIM = 1) ! at risk
    Nj = sum(pr2surv(casesIn,:), DIM = 1) ! at risk

      !PRINT *, "******************** calcValueSingle ********************"
      !PRINT *, "STEP2RE estimate MFF"
      !PRINT *, "casesIn with size:", size(casesIn)
      !PRINT *, casesIn
      !PRINT *, "nt", nt
      !PRINT *, "nt_death", nt_death
      !PRINT *, "dim of pr2", shape(pr2)
      !PRINT *, "dim of pr2surv", shape(pr2surv)
      !PRINT *, "dim of pr", shape(pr)
      !PRINT *, "dim of prsurv", shape(prsurv)
      !PRINT *, "delta with size", size(delta)
      !PRINT *, "delta_m with size", size(delta_m)
      !PRINT *
    
    IF ((Nj_m(1) < Nj_m(2)) .OR. (Nj_m(1) .LT. 1)) THEN 
      PRINT *, "Looping for pr2surv: at risk for survival for first 10 timepoints"
      do j = 1, size(casesIn)
        testj = casesIn(j)
        do t = 10, 20  ! Loop over time indices 1 to 10
          PRINT '(A, I2, A, I2, A, F8.5)', 'pr2surv(case', testj, ', time', t, ') = ', pr2surv(testj, t)
        end do
        PRINT '(A, I2, A, I2, A, F8.5)', 'pr2surv(case', testj, ', time', nt_death, ') = ', pr2surv(testj, nt_death)
      end do
    END IF


    ! Print Nj_survival array with index, rounded to 1 decimal point
    !PRINT *, 'Nj with size ', size(Nj)
    !DO i = 1, nt_death
    !  PRINT '(A, I2, A, F6.1)', 'Nj(', i, ') = ', Nj(i)
    !END DO
    !PRINT *
    !DO i = 1, size(Nj_m)
    !  PRINT '(A, I2, A, F6.1)', 'Nj_m(', i, ') = ', Nj_m(i)
    !END DO
    !PRINT *


    ! We now calculate the number of events at each time point using pr {nt}
    DO i = 1, nt
      ! here: delta depends on if you are doing survival or CR
      Oj_m(i) = sum(pr(casesIn, i)*delta_m(casesIn))
      !PRINT *, "((((((((((((( CHECKING OJ_m(i) ))))))))))))): ", Oj_m(i)
    END DO

    DO i = 1, nt_death
      Oj(i) = sum(prsurv(casesIn, i)*delta(casesIn)) ! HOW MANY FAILURES AT EACH TIME POINT, * IS ELEMENT-WISE MULTIPLICATION.
      !PRINT *, "((((((((((((( CHECKING OJ(i) ))))))))))))): ", Oj(i)
    END DO

    ! Mean Frequency Function
    ! PRINT *, "calling mean freq func within calcvaluesingle"
    CALL MeanFreqFunc(nt, nt_death, surv_tp, end_tp, &
                    Nj_m, Oj_m, Nj, Oj, &
                    survRE, dRhat, Func, dmu) ! name mu to be Func
    mean = Func(nt)

    IF (mean < 0.0 .OR. IEEE_IS_NAN(mean)) THEN
        PRINT *, "ERROR IN MEAN:", mean
        PRINT *, "dt: "
        PRINT *, dt
        PRINT *, "nt = ", nt
        PRINT *, "mu Estimate MFF: "
        PRINT *, Func
        PRINT *
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

  PRINT *, "******************** getCovariate ********************"
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
    PRINT *, "getCovariate LINE 4170: CALL calcValueSingle"
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

  PRINT *, "tcalculateValue: CALL calcValueSingle"
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
! forestFunc, real(:), function averaged over forest
! forestMean, real, Phase1/2CR: integration from 0 to tau of function averaged over forest, Phase2RE: mff(tau)
! forestProb, real, probability averaged over forest, not applicable for Phase2RE.
SUBROUTINE tsurvTree(forestFunc, forestMean, forestProb)
  IMPLICIT NONE

  REAL(dp), DIMENSION(1:nAll_surv*nt), INTENT(OUT) :: forestFunc
  REAL(dp), DIMENSION(1:nAll_surv), INTENT(OUT) :: forestMean
  REAL(dp), DIMENSION(1:nAll_surv), INTENT(OUT) :: forestProb

  INTEGER :: i, iTree, j, k, lft, m, nc, ncur, splitFound, splitVar
  INTEGER, allocatable :: indices(:), jdex(:), jdexRE(:), xrand(:)
  INTEGER, allocatable, DIMENSION(:) :: indices_surv_RE, jdex_people_RE, jdex_record_RE
  
  INTEGER, DIMENSION(1:np) :: newstat, pindices
  INTEGER, DIMENSION(1:np, 1:nrNodes) :: cstat
  INTEGER, DIMENSION(1:nrNodes, 1:2) :: stm
  INTEGER, DIMENSION(1:nAll_surv) :: allStatus
  INTEGER, DIMENSION(:), ALLOCATABLE :: ind, indOut, indOutRE, indRE, indOut_people_RE, indOut_record_RE
  INTEGER, DIMENSION(:), ALLOCATABLE :: leftCases, rightCases, leftCasesRE, rightCasesRE, pind !leftPeople, rightPeople, 
  INTEGER, DIMENSION(:), ALLOCATABLE :: ind_people_RE, record_surv_RE
  INTEGER, DIMENSION(:), ALLOCATABLE :: unique_ind_people_RE, record_ind, record_ind_trt

  REAL(dp) :: srs
  REAL(dp), DIMENSION(1:nLevs) :: cutoffBest
  REAL(dp), DIMENSION(1:nrNodes) :: mean, Prob
  REAL(dp), DIMENSION(1:nAll) :: xm
  REAL(dp), DIMENSION(1:nt, 1:nrNodes) :: Func
  REAL(dp), DIMENSION(1:nrNodes, 1:(5+nLevs)) :: nMatrix
  REAL(dp), DIMENSION(1:nt, 1:nAll_surv) :: tforestSurvFunc

  LOGICAL, DIMENSION(1:np) :: cand
  LOGICAL, DIMENSION(1:nAll) :: tst

  LOGICAL :: are_equal, print_check

  integer :: n_subj_from_records, size_unique_ind_people_RE    ! Declare the integer that stores the number of unique elements
  INTEGER, DIMENSION(:), ALLOCATABLE ::  unique_id_RE2

  INTEGER :: n_leftCases, n_rightCases
  INTEGER, DIMENSION(:), ALLOCATABLE :: unique_leftCases, unique_rightCases

  INTEGER, ALLOCATABLE :: sampledArray(:)    ! Expanded sampled IDs
  INTEGER, ALLOCATABLE :: sampledArray_index(:) ! Indices of sampled records
  INTEGER, ALLOCATABLE, DIMENSION(:) :: sampledArray_new
  INTEGER, ALLOCATABLE, DIMENSION(:) :: sampledArray_index_new

  print_check = .FALSE.
  are_equal = .TRUE.
  
  !IF (isPhase2RE) PRINT *, "******************** survTree ********************"

  if (print_check) THEN
    PRINT *, "nAll: ", nAll
    PRINT *, "sampleSize", sampleSize
  end if

  tforestSurvFunc = 0.d0
  forestMean = 0.d0
  forestProb = 0.d0


  DO iTree = 1, nTree ! do iTree

    Func = 0.d0
    mean = 0.d0
    Prob = 0.d0
    nMatrix = 0.0
    allStatus = 1

    ! sample data and set local variables x, pr, and delta to the selected subset
    IF (replace .EQ. 1) THEN
      PRINT *, "THIS IS NOT CORRECTLY CODED UP YET FOR PHASE2RE"
      IF (isPhase2RE) STOP
      PRINT *, "sample with replacement"
      xrand = sampleWithReplace(nAll, sampleSize)
      n = sampleSize ! the number of cases/individuals to sample for each tree
      x = xAll(xrand,:)
      pr = prAll(xrand,:)
      IF (isPhase2RE) THEN
        n_surv = sampleSize ! number of people to sample for each tree (only applicable for Phase2RE) ! this differs from n_subj_from_records later, which reflects the equivlant number of people from sampleSize (records) criteria
        id_RE2 = id_RE(xrand)
        pr2 = pr2All(xrand,:)
        prsurv = prsurvAll(xrand,:)
        pr2surv = pr2survAll(xrand,:)
      END IF
      delta = deltaAll(xrand)
      delta_m = deltaAll_m(xrand)

    ELSE IF (nAll .NE. sampleSize) THEN

      IF (isPhase2RE) THEN ! Phase 2 RE sampling w/o replacement
        IF (allocated(sampledArray)) THEN
          deallocate(sampledArray)
        END IF
        IF (allocated(sampledArray_index)) THEN
            deallocate(sampledArray_index)
        END IF
        CALL sampleExpandWithoutReplacement(nAll, nAll_surv, sampleSize, id_RE, size(id_RE), sampledArray, sampledArray_index)
        !PRINT *, "size(sampledArray):", size(sampledArray)
        !PRINT *, "before assingmnet for smapling array"
        !PRINT *, isAlloc_sampledArray
        IF (allocated(sampledArray_new)) DEALLOCATE(sampledArray_new)
        ALLOCATE(sampledArray_new(1:size(sampledArray)))
        isAlloc_sampledArray = .TRUE.
        !PRINT *, "after assignment for smapling array"
        !PRINT *, isAlloc_sampledArray
        CALL create_new_vector(sampledArray, size(sampledArray), sampledArray_new)
        sampledArray_index_new = (/(i, i = 1, size(sampledArray_index))/)
        n = size(sampledArray_new) ! records
        n_surv = sampleSize ! number of people to sample for each tree (only applicable for Phase2RE) ! this differs from n_subj_from_records later, which reflects the equivlant number of people from sampleSize (records) criteria
        id_RE2 = id_RE(sampledArray_index)
        row_RE2 = row_RE(sampledArray_index)
        pr2 = pr2All(sampledArray_index,:)
        prsurv = prsurvAll(sampledArray_index,:)
        pr2surv = pr2survAll(sampledArray_index,:)
        x = xAll(sampledArray_index,:)
        pr = prAll(sampledArray_index,:)
        delta = deltaAll(sampledArray_index)
        delta_m = deltaAll_m(sampledArray_index)
          IF (print_check) THEN
                  IF (isPhase2RE) THEN
        PRINT *, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        PRINT *, "Sampled Array (Expanded IDs) with size:", size(sampledArray)
        PRINT *, sampledArray
        PRINT *, "Sampled Array Indices with size:", size(sampledArray_index)
        PRINT *, sampledArray_index
        PRINT *, "new Sampled Array with size", size(sampledArray_new)
        PRINT *, sampledArray_new
        PRINT *, "new Sampled Array Indices with size", size(sampledArray_index_new)
        PRINT *, sampledArray_index_new
        PRINT *, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        END IF
            PRINT *, "size(id_RE)",size(id_RE)
            PRINT *, "size(id_RE2)",size(id_RE2)
            PRINT *, "size(deltaAll):",size(deltaAll)
            PRINT *, "size(delta):",size(delta)
            PRINT *, "size(pr2All, 1)",size(pr2All, 1)
            PRINT *, "size(pr2,1)",size(pr2,1)
            PRINT *, "size(prsurvAll, 1)",size(prsurvAll, 1)
            PRINT *, "size(prsurv,1)",size(prsurv,1)
            PRINT *, "size(pr2survAll, 1)",size(pr2survAll, 1)
            PRINT *, "size(pr2surv,1)",size(pr2surv,1)
            PRINT *, "nAll_surv:", nAll_surv
            PRINT *, "sampleSize:", sampleSize
          END IF
      ELSE
        xrand = sampleWithoutReplace(nAll, sampleSize)
        n = sampleSize
        x = xAll(xrand,:)
        pr = prAll(xrand,:)
        delta = deltaAll(xrand)
        delta_m = deltaAll_m(xrand)
        !PRINT *, "sample without replacement"
        !PRINT *, "nAll:", nAll
        !PRINT *, "sampleSize:", sampleSize    
        !PRINT *, "xrand has size ", size(xrand) 
        !PRINT *, xrand
      END IF
    ELSE
      IF (isPhase2RE) THEN
        PRINT *, "THIS IS NOT CORRECTLY CODED UP YET FOR PHASE2RE"
        STOP
      ELSE
        !PRINT *, "replace = ", replace, "and nAll = sampleSize:", nAll, " = ", sampleSize
        n = sampleSize
        !PRINT *, "n = ", n
        x = xAll
        pr = prAll
        delta = deltaAll
        delta_m = deltaAll_m
      END IF
      IF (isPhase2RE) THEN
        n_surv = sampleSize ! number of people to sample for each tree (only applicable for Phase2RE) ! this differs from n_subj_from_records later, which reflects the equivlant number of people from sampleSize (records) criteria
        id_RE2 = id_RE
        pr2 = pr2All
        prsurv = prsurvAll
        pr2surv = pr2survAll
      END IF
    END IF

    IF (isPhase2RE) THEN
      !PRINT *, "sampleSize", sampleSize
      if (print_check) PRINT *, "Number of REs:", n
      ! Call the subroutine to find unique elements
      call find_unique(id_RE2, unique_id_RE2, &
      n_subj_from_records)
      !PRINT *, "id_RE2:"
      !PRINT *, id_RE2
      !PRINT *, "unique_id_RE2"
      !PRINT *, unique_id_RE2
      !print *, 'Number of subjects: ', n_subj_from_records
    END IF

    ! cutoff for identifying covariates to be explored
    srs = stratifiedSplit / REAL(np)

    ! indices for all cases in this tree
    IF (isPhase2RE) THEN
      IF (isAlloc_ind) DEALLOCATE(indices, jdex, indices_surv_RE, jdex_people_RE, jdex_record_RE)
      ALLOCATE(indices(1:size(sampledArray_new)))
      ALLOCATE(jdex(1:sampleSize))
      ALLOCATE(indices_surv_RE(1:size(sampledArray_new)))
      ALLOCATE(jdex_people_RE(1:sampleSize))
      ALLOCATE(jdex_record_RE(1:sampleSize))
      isAlloc_ind = .TRUE.

      indices = sampledArray_new
      jdex = indices !(/(i, i = 1, size(indices))/)
      !PRINT *, "indices with size", size(indices)
      !PRINT *, indices
      jdexRE = sampledArray_index_new
    ELSE 
      indices = (/(i,i=1,n)/)
      jdex = indices
      jdexRE = indices
    END IF
    IF (isPhase2RE) THEN
    !! to do: this is wrong and needs to be updated to the correct thing
    indices_surv_RE = id_RE2
    jdex_people_RE = indices_surv_RE
    jdex_record_RE = row_RE2

    !  ! check that indices aka 1:n is the same as the other inputs always
      IF ((SIZE(delta) /= SIZE(delta_m)) .OR. &
      (SIZE(delta) /= SIZE(id_RE2)) .OR. &
      (SIZE(delta) /= SIZE(indices))) THEN
        PRINT *, "Error: Lengths of delta, delta_m, id_RE2, and indices are not equal."
        PRINT *, "Length of delta =", SIZE(delta)
        PRINT *, "Length of delta_m =", SIZE(delta_m)
        PRINT *, "Length of id_RE2 =", SIZE(id_RE2)
        PRINT *, "Length of indices =", SIZE(indices)
        STOP "ERROR: LINE 2604: indices/1:n is not equal to delta,delta_m,id_RE2 lengths!!"
      END IF
    END IF

    ! indices for all covariates
    pindices = (/(i,i=1,np)/)

    !! initialize first node

    ! calculate survival function and mean survival of the first node
    ! IF (isPhase1) PRINT *, "calculate survival function and mean survival of the first node"
    ! IF (isPhase2CR) PRINT *, "calculate CIF function and mean CIF of the first node"
    !PRINT *, "tsurvTree LINE 4388: CALL calcValueSingle"
    !IF (isPhase2RE) THEN
    !  PRINT *, "indices:"
    !  PRINT *, indices
    !  PRINT *, "n:", n
    !END IF

    IF (isPhase2RE) THEN
      !PRINT *, "Initializing calcvaluesingle"
      CALL calcValueSingle(size(sampledArray_index_new), sampledArray_index_new, Func(:,1), mean(1))
    ELSE
      CALL calcValueSingle(n, indices, Func(:,1), mean(1))
    END IF

    IF (isSurvival) THEN
      ! estimate survival/cif probability at SurvivalTime/CIFTime
      Prob(1) = Func(sIndex,1) * (1.d0 - sFraction) + &
                  & Func(sIndex+1,1) * sFraction
      IF (Prob(1) .LT. 1d-8) Prob(1) = 0.d0
    END IF
    ! For endpoint = CR: Phase 1: delta = delta_m
    ! For endpoint = RE: Phase 1: delta_m is just delta
    ! nMatrix: DIMENSION(1:nrNodes, 1:(5+nLevs))

    ! determine if the node can split based on basic minimum requirements
    IF (isPhase1) THEN 
      ! determine if the node can split based on basic minimum requirements
      if (n .LE. nodeSizeSurv .OR. sum(delta(indices)) .LE. 1) THEN
        nMatrix(1,1) = -1
      ELSE
        nMatrix(1,1) = -2
      END IF
    END IF

    IF (isPhase2) THEN
      IF (n .LE. nodeSizeEnd .OR. sum(delta_m(indices)) .LE. 1) THEN
        IF (isPhase2RE) THEN
          IF (n_surv .LE. nodeSizeSurv .OR. sum(delta(indices)) .LE. 1) THEN
            nMatrix(1,1) = -1
          ELSE
            nMatrix(1,1) = -2
          END IF
        ELSE 
          nMatrix(1,1) = -1 ! isPhase2CR:
        END IF
      ELSE
        nMatrix(1,1) = -2
      END IF
    END IF

    cstat(:,1) = 0

    ! start and finish locations of indices in node
    stm(1,1) = 1
    if (isPhase2RE) THEN
      stm(1,2) = size(jdex)
    ELSE 
      stm(1,2) = n
    END IF

    ! location of most recent storage location in matrices/vectors
    ! ncur is incremented when a node successfully splits indicating the
    ! location in the nodes list where base information for the each daughter
    ! is stored
    ncur = 1

    ! --===----= here

    !PRINT *, "nrNodes: ", nrNodes
    !PRINT *, "Starting do loop for each node..."
    DO k = 1, nrNodes
      !if (isPhase2RE) THEN
      !PRINT *, "####################################################"
      !PRINT *, "####################################################"
      !PRINT *, "############### node:", int(k), "##################"
      !PRINT *, "####################################################"
      !PRINT *, "####################################################"
      !end if

      !if (isPhase2RE) PRINT *, "first ncur:", ncur
      !if (isPhase2RE) PRINT *, "nrNodes - 2:", nrNodes-2

      ! if k is beyond current node count or
      ! current node count at limit, break from loop
      IF (k .GT. ncur .OR. ncur .GT. (nrNodes - 2)) EXIT
      !PRINT *, "DID NOT EXIT: k is not beyond current node count"

      ! if node is not to be split, cycle to next node
      IF (nint(nMatrix(k,1)) .EQ. -1) CYCLE
      
      ! indices for cases contained in node
      ind = jdex(stm(k,1):stm(k,2))
      IF (isPhase2RE) THEN
        indRE = jdexRE(stm(k,1):stm(k,2))
      ELSE
        indRE = ind
      END IF
      !PRINT *, "TEST8"

      IF (isPhase2RE) THEN

        ind_people_RE = jdex_people_RE(stm(k,1):stm(k,2))
        record_ind = sampledArray_index_new(stm(k,1):stm(k,2))
        record_ind_trt = sampledArray_index(stm(k,1):stm(k,2))
        record_surv_RE = jdex_record_RE(stm(k,1):stm(k,2)) !row_RE2(stm(k,1):stm(k,2))

        IF (print_check) THEN
        PRINT *, "record_ind with size", size(record_ind)
        PRINT *, record_ind
        PRINT *, "record_ind_trt with size", size(record_ind_trt)
        PRINT *, record_ind_trt
            PRINT *, "stm(k,1); k = ", k
            PRINT *, stm(k,1)
            PRINT *, "stm(k,2); k = ", k
            PRINT *, stm(k,2)
            PRINT *, "jdex with size", size(jdex)
            PRINT *, jdex
            PRINT *, "ind = jdex(stm(k,1):stm(k,2)) with size", size(jdex(stm(k,1):stm(k,2)))
            PRINT *, jdex(stm(k,1):stm(k,2))
            PRINT *, "id_RE2"
            PRINT *, id_RE2
            PRINT *, "id_RE2(stm(k,1):stm(k,2))"
            PRINT *, id_RE2(stm(k,1):stm(k,2))
            PRINT *, "stm(1:k, :)"
            PRINT *, stm(1:k, :)
            PRINT *
            PRINT *, "jdex(stm(k,1):stm(k,2))"
            PRINT *, jdex(stm(k,1):stm(k,2))
            PRINT *
            PRINT *, "jdex_people_RE(stm(k,1):stm(k,2))"
            PRINT *, jdex_people_RE(stm(k,1):stm(k,2))
            PRINT *
            PRINT *, "sampledArray_index"
            PRINT *, sampledArray_index
            PRINT *, "sampledArray_new"
            PRINT *, sampledArray_new
            PRINT *, "ind"
            PRINT *, ind
          END IF
      END IF

      ! if there are deficient variables, use only these variables
      cand = cstat(:,k) .LT. floor(srs * sum(cStat(:,k)))
      pind = pack(pindices, cand)
      IF (size(pind) .EQ. 0) pind = pindices
      !PRINT *, "TEST9"

      ! split cases
      indOut = ind ! elements of casesIn that go left; ind if yes, 0 otherwise

      if (isPhase2RE) THEN
        indOutRE = indRE ! fix this
        indOut_people_RE = ind_people_RE
        indOut_record_RE = record_surv_RE
        ! get the total number of subjects (needed for 2RE. Same as size(ind) for 1/2CR.)
        call find_unique(ind_people_RE, &
        unique_ind_people_RE, size_unique_ind_people_RE)
      ELSE 
        indOutRE = ind
        indOut_people_RE = ind
        indOut_record_RE = ind
        unique_ind_people_RE = ind
        size_unique_ind_people_RE = size(ind)
      END IF      
    
      IF (print_check) THEN
          IF (isPhase2RE) THEN
          PRINT *, "---------- Line 5153: call tfindSplit ----------"
          PRINT *, "ind with size", size(ind)
          !PRINT *, ind
          PRINT *, "indRE with size", size(indRE)
          !PRINT *, indRE
          PRINT *, "ind_people_RE with size", size(ind_people_RE), " and ", size_unique_ind_people_RE
          !PRINT *, ind_people_RE
          !PRINT *, "record_ind with size", size(record_ind)
          !PRINT *, record_ind
          PRINT *, "new Sampled Array Indices with size", size(sampledArray_index_new)
          !PRINT *, sampledArray_index_new
        END IF ! end of isPhase2RE
      END IF ! end of print_check

      CALL tfindSplit(k, size(ind), ind, indRE, &
                    & size_unique_ind_people_RE, ind_people_RE, &
                    & sampledArray_index_new, &
                    & record_ind, record_ind_trt, record_surv_RE, &
                    & size(pind), pind, splitVar, cutoffBest, &
                    & splitFound, indOut, indOutRE, &
                    & indOut_people_RE, indOut_record_RE, nc, lft)
      !IF (isPhase2RE) PRINT *, "---------- Line 5174: end of tfindSplit ----------"


      IF (print_check) THEN
        IF (isPhase2RE) THEN
          IF (splitFound .EQ. 1) THEN
          PRINT *, "^^^^^^^^^^^^^^^^^^^^^^^^^^^"
            PRINT *, "splitFound:", splitFound
            !PRINT *, "size of left group: lft:", lft
            !PRINT *, "### end of tfindSplit for node=", k
            !PRINT *, "indOut with size:", size(indOut)
            !PRINT *, indOut
              !PRINT *, "indOutRE with size:", size(indOutRE)
              !PRINT *, indOutRE
              !PRINT *, "indOut_people_RE with size:", size(indOut_people_RE)
              !PRINT *, indOut_people_RE
            PRINT *, "vvvvvvvvvvvvvvvvvvvvvvvvv"
          END IF
        END IF
      END IF


      IF (splitFound .EQ. 0 ) THEN
        ! if no split available, set node k as terminal node
        !IF (isPhase2RE) PRINT *, "Line 2635: No split found for node ", k, " so nMatrix(k,1) = -1 (terminal node)"
        nMatrix(k,1) = -1
        CYCLE
      END IF
      !PRINT *, "TEST11"

      ! set node k to be interior (i.e. has split)
      !IF (isPhase2RE) PRINT *, "Line 2642: Setting node ", k, "to be interior (has split), so nMatrix(k,1) = -3"
      nMatrix(k,1) = -3
      !PRINT *, "TEST12"

      ! add split information to node
      nMatrix(k,4) = pindices(splitVar)
      nMatrix(k,5) = nc
      nMatrix(k,6:(6+nc-1)) = cutoffBest(1:nc)
      nMatrix(k,2) = ncur + 1
      nMatrix(k,3) = ncur + 2

      IF (k == 2 .AND. isPhase2RE .AND. print_check) THEN
        PRINT *, "nMatrix"
        PRINT *, nMatrix(1:k,:)
      END IF

      !PRINT *, "increment the times the variable was used in a split."
      ! increment the times the variable was used in a split
      newstat = cStat(:,k)
      newstat(nint(nMatrix(k,4))) = newstat(nint(nMatrix(k,4))) + 1

      ! store new case order in jdex
      jdex(stm(k,1):stm(k,2)) = indOut
      IF (isPhase2RE) THEN
        ! store new record order in jdex
        jdexRE(stm(k,1):stm(k,2)) = indOutRE
        ! store original ID labels
        !PRINT *, "shape of jdex_people_RE", shape(jdex_people_RE)
        jdex_people_RE(stm(k,1):stm(k,2)) = indOut_people_RE
        jdex_record_RE(stm(k,1):stm(k,2)) = indOut_record_RE
      END IF

      !! left node
      !PRINT *, "!!!!! left node !!!!!"

      !PRINT *, "left node ncur:", ncur
      ncur = ncur + 1
      !PRINT *, "left node ncur + 1:", ncur

      ! index boundaries for cases in left node
      stm(ncur,1) = stm(k,1)
      stm(ncur,2) = stm(k,1) + lft - 1
      !IF (isPhase2RE) PRINT *, "stm(", ncur, ",2) = ", stm(k,1) + lft - 1

      leftCases = jdex(stm(ncur,1):stm(ncur,2))
      !IF (isPhase2RE) PRINT *, "leftCases with size ", size(leftCases)
      !IF (isPhase2RE) PRINT *, leftCases

      ! get basic node information for left daughter
      ! IF (isPhase1) PRINT *, "get basic node information for left daughter."
      !IF (isPhase2RE) PRINT *, "tsurvTree: Line 4616: CALL calcValueSingle"
      IF (isPhase2RE) THEN
        leftCasesRE = jdexRE(stm(ncur,1):stm(ncur,2))
        !PRINT *, "leftCasesRE with size ", size(leftCasesRE)
        !PRINT *, leftCasesRE
        !PRINT *, "calcValueSingle for leftRE"
        CALL calcValueSingle(size(leftCasesRE), leftCasesRE, Func(:,ncur), &
                          & mean(ncur))
      ELSE
        CALL calcValueSingle(size(leftCases), leftCases, Func(:,ncur), &
                          & mean(ncur))
      END IF

      ! estimate survival probability at SurvivalTime
      IF (isSurvival) THEN
      Prob(ncur) = Func(sIndex,ncur) * (1.d0 - sFraction) + &
             & Func(sIndex+1,ncur) * sFraction
      IF (Prob(ncur) .LT. 1d-8) Prob(1) = 0.d0
      END IF

      IF (isPhase1) THEN
        IF (size(leftCases) .LE. nodeSizeSurv .OR. sum(delta(leftCases)) .LE. 1) THEN
          ! if the number of cases in the node is at or below the minimum required
          ! or the number of uncensored OS/Death event is only 1
          ! status is terminal
          nMatrix(ncur,1) = -1
        ELSE
          nMatrix(ncur,1) = -2
        END IF
      ELSE ! PHASE2
        IF (isPhase2CR) THEN
          IF (size(leftCases) .LE. nodeSizeEnd .OR. sum(delta_m(leftCases)) .LE. 1) THEN
            ! if the number of cases in the node is at or below the minimum required
            ! or the number of uncensored PC/RE event is at most 1
            ! status is terminal
            nMatrix(ncur,1) = -1
          ELSE
            nMatrix(ncur,1) = -2
          END IF
        ELSE ! RE
          CALL find_unique(leftCases, unique_leftCases, n_leftCases)
          IF (n_leftCases .LE. nodeSizeEnd .OR. sum(delta_m(leftCases)) .LE. 1) THEN
            ! if the number of cases in the node is at or below the minimum required
            ! or the number of uncensored PC/RE event is at most 1
            ! status is terminal
            nMatrix(ncur,1) = -1
          ELSE
            nMatrix(ncur,1) = -2
          END IF
        END IF
      END IF

      cstat(:,ncur) = newstat

      !! right node
      !PRINT *, "!!!!! right node !!!!!"
      !PRINT *, "right node ncur:", ncur
      ncur = ncur + 1
      !PRINT *, "right node ncur + 1:", ncur

      ! index boundaries for cases in right node
      stm(ncur,1) = stm(k,1) + lft
      stm(ncur,2) = stm(k,2)
      !IF (isPhase2RE) PRINT *, "lft:", lft
      !IF (isPhase2RE) PRINT *, "stm(",ncur,",1) = ",stm(k,1) + lft

      ! retrieve left and right cases
      rightCases = jdex(stm(ncur,1):stm(ncur,2))
      !IF (isPhase2RE) PRINT *, "rightCases with size", size(rightCases)
      !IF (isPhase2RE) PRINT *, rightCases

      IF (isPhase2RE) THEN
        ! retrieve left and right cases
        rightCasesRE = jdexRE(stm(ncur,1):stm(ncur,2))
        !PRINT *, "rightCasesRE with size", size(rightCasesRE)
        !PRINT *, rightCasesRE
        ! calculate survival function and mean survival time for right node
        CALL calcValueSingle(size(rightCasesRE), rightCasesRE, Func(:,ncur), &
                          & mean(ncur))
        !PRINT *, "----------"   
      ELSE
        ! calculate survival function and mean survival time for right node
        ! IF (isPhase1) PRINT *, "calculate survival function and mean survival time for right node."
        !IF (isPhase2RE) PRINT *, "tsurvTree: Line 5356: CALL calcValueSingle"
        CALL calcValueSingle(size(rightCases), rightCases, Func(:,ncur), &
                          & mean(ncur))
        !PRINT *, "rightCases"
        !PRINT *, rightCases
        !PRINT *, "----------"        
      END IF 

      !PRINT *, "################# Line 4506"
      IF (isSurvival) THEN
        Prob(ncur) = Func(sIndex,ncur) * (1.d0 - sFraction) + &
                       & Func(sIndex+1,ncur) * sFraction
        IF (Prob(ncur) .LT. 1d-8) Prob(1) = 0.d0
      END IF

      IF (isPhase1) THEN
        IF (size(rightCases) .LE. nodeSizeSurv .OR. &
          & sum(delta(rightCases)) .LE. 1) THEN
          ! if the number of cases in the node is at or below the minimum required
          ! or the number of uncensored event (from cause M) is only 1
          ! status is terminal
          nMatrix(ncur,1) = -1
        ELSE
          nMatrix(ncur,1) = -2
        END If
      END IF

      IF (isPhase2) THEN
        IF (isPhase2CR) THEN
          IF (size(rightCases) .LE. nodeSizeEnd .OR. &
            & sum(delta_m(rightCases)) .LE. 1) THEN
            ! if the number of cases in the node is at or below the minimum required
            ! or the number of uncensored event (from cause M) is only 1
            ! status is terminal
            nMatrix(ncur,1) = -1
          ELSE
            nMatrix(ncur,1) = -2
          END IF
        ELSE ! RE
          CALL find_unique(rightCases, unique_rightCases, n_rightCases)
          IF (n_rightCases .LE. nodeSizeEnd .OR. &
            & sum(delta_m(rightCases)) .LE. 1) THEN
            ! if the number of cases in the node is at or below the minimum required
            ! or the number of uncensored event (from cause M) is only 1
            ! status is terminal
            nMatrix(ncur,1) = -1
          ELSE
            nMatrix(ncur,1) = -2
          END IF
        END IF
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

      DO j = 1, nAll_surv
        IF (allStatus(j) .NE. k) CYCLE
        IF (tst(j)) THEN
          allStatus(j) = nint(nMatrix(k,2))
        ELSE
          allStatus(j) = nint(nMatrix(k,3))
        END IF
      END DO

    END DO ! end of k = 1 to nrNodes

    ! ensure that all nodes that are not "interior" are "terminal"
    WHERE (nMatrix(:,1) .EQ. -2) nMatrix(:,1) = -1


    !PRINT *, "end ncur:", ncur
    !PRINT *, "mean"
    !PRINT *, mean
    !PRINT *
    !PRINT *, "allStatus"
    !PRINT *, allStatus

    trees(iTree)%Func = Func(:,1:ncur)
    trees(iTree)%mean = mean(1:ncur)
    trees(iTree)%Prob = Prob(1:ncur)
    trees(iTree)%matrix = nMatrix(1:ncur,:)
    trees(iTree)%nNode = ncur

    tforestSurvFunc = tforestSurvFunc + Func(:,allStatus)
    forestMean = forestMean + mean(allStatus)
    forestProb = forestProb + Prob(allStatus)

    !IF (isPhase2RE) PRINT *, "end of iTree:", iTree

  END DO ! end do-loop for iTree

  ! -- to here
  !IF (isPhase2RE) PRINT *, "################# end of iTree do-loop at line 5452 #################"

  ! This change is to eliminate a strange lto warning from R
  ! forestFunc = reshape(tforestSurvFunc, (/nt*nAll/)) / nTree
  j = 0
  DO i = 1, SIZE(tforestSurvFunc,2)
    forestFunc(j+1:j+SIZE(tforestSurvFunc,1)) = tforestSurvFunc(:,i)
    j = j + SIZE(tforestSurvFunc,1)
  END DO

  forestFunc = forestFunc / nTree
  forestMean = forestMean / nTree
  forestProb = forestProb / nTree

  !PRINT *
  !PRINT *
  !PRINT *
  !PRINT *, "forestFunc: ", forestFunc
  !PRINT *
  !PRINT *
  !PRINT *
  !PRINT *, "forestMean with size", size(forestMean)
  !PRINT *, forestMean
  PRINT *
  PRINT *
  PRINT *
  PRINT *
  PRINT *
  PRINT *
  !PRINT *, "forestProb: ", forestProb

!IF (isPhase2RE) PRINT *, "END OF SUBROUTINE TSURVTREE"
  

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
! t_surv_tp, real(:), the time points for survival (death)
! t_end_tp, real(:), the time points for phase 2 (recurrent events)
! t_nt, integer, the number of observed unique event time points including 0 and tau; death for isPhase1/isPhase2CR and recurrent events for isPhase2RE
! t_nt_death, integer, the number of observed unique death time points including 0 and tau; needed for isPhase2RE
! t_dt, real(:), the time differences between time points
! t_dt_death, real(:), the time differences between death time points
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
SUBROUTINE setUpBasics(t_surv_tp, t_end_tp, t_nt, t_nt_death, t_dt, t_dt_death, t_rs, t_ERT, t_uniformSplit, &
                     & t_nodeSize, t_nodeSizeSurv, &
                     & t_minEvent, t_minEventSurv, &
                     & t_rule, t_sIndex, t_sFraction, &
                     & t_stratifiedSplit, t_replace)

  USE INNERS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: t_nt
  INTEGER, INTENT(IN) :: t_nt_death
  REAL(dp), DIMENSION(1:t_nt), INTENT(IN) :: t_end_tp
  REAL(dp), DIMENSION(1:t_nt_death), INTENT(IN) :: t_surv_tp
  REAL(dp), DIMENSION(1:t_nt), INTENT(IN) :: t_dt
  REAL(dp), DIMENSION(1:t_nt_death), INTENT(IN) :: t_dt_death
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

  !PRINT *, "t_nt:", t_nt
  !PRINT *, "t_nt_death:", t_nt_death
  !PRINT *
  !PRINT *, "t_surv_tp:"
  !PRINT *, t_surv_tp
  !PRINT *
  !PRINT *, "t_end_tp:"
  !PRINT *, t_end_tp
  !PRINT *
  
  nt = t_nt
  nt_death = t_nt_death
  surv_tp = t_surv_tp
  end_tp = t_end_tp

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
      PRINT *, "itrSurv.f90: Line 5127: nt should equal nt_death for Phase1 and Phase2CR so stopping"
      STOP
   END IF
END IF

  !PRINT *, dtAllocated
  IF (dtAllocated) DEALLOCATE(dt)
  ALLOCATE(dt(1:nt))
  dtAllocated = .TRUE.
  !PRINT *, "After assignment"
  !PRINT *, dtAllocated
  ! write(*, *) 'Check2 Fortran'

  dt = t_dt
  dt_death = t_dt_death

  rs = t_rs
  ERT = t_ERT
  uniformSplit = t_uniformSplit
  nodeSizeEnd = t_nodeSize
  nodeSizeSurv = t_nodeSizeSurv
  minEventEnd = t_minEvent
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
! t_rowvec, integer(:), the row label for Phase2 dataset for RE
! t_person_ind, integer(:), indicator for person in RE setting
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
SUBROUTINE setUpInners(t_n, t_n_surv, t_idvec, t_rowvec, t_person_ind, t_np, t_x, &
                      & t_pr, t_pr2, t_pr2surv, t_prsurv, t_ord_causeind, t_ord_response, &
                      & t_delta, t_delta_m, t_mTry, t_nCat, &
                      & t_sampleSize, t_nTree, t_nrNodes)

  USE INNERS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: t_n, t_n_surv
  INTEGER, DIMENSION(1:t_n), INTENT(IN) :: t_idvec
  INTEGER, DIMENSION(1:t_n), INTENT(IN) :: t_rowvec
  INTEGER, DIMENSION(1:t_n), INTENT(IN) :: t_person_ind
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

  INTEGER :: i, testi, j
  LOGICAL :: are_equal


 !PRINT *, "******************** setUpInners ********************"

  isAlloc_sampledArray = .FALSE.
  isAlloc_ind = .FALSE.

  nAll = t_n
  nAll_surv = t_n_surv
  !print *, "nAll is:", nAll
  ord_causeindAll = t_ord_causeind
  id_RE = t_idvec
  row_RE = t_rowvec
  person_ind = t_person_ind
  ord_responseAll = t_ord_response
  np = t_np

  !IF (isPhase2RE) THEN
  !  DO j = 1, MIN(10, t_n)
  !    PRINT *, 't_n =', j, ': '
  !    DO i = 1, nt_death
  !        WRITE(*, '(F5.1, 1X)', ADVANCE='NO') t_prsurv((i - 1) * t_n + j)
  !    END DO
  !    PRINT *  ! New line after each t_n
  !  END DO
  !END IF

  !PRINT *, "Before assignment alloc_list"
  !PRINT *, isAlloc_list
  IF (isAlloc_list) THEN
    DEALLOCATE(xAll, prAll, deltaAll, deltaAll_m, nCat, forest%Func, forest%mean,  &
             & forest%Prob, trees)
  END IF

  ALLOCATE(xAll(1:nAll,1:np))
  ALLOCATE(prAll(1:nAll, 1:nt))
  ALLOCATE(deltaAll(1:nAll))
  ALLOCATE(deltaAll_m(1:nAll))
  ALLOCATE(nCat(1:np))

  isAlloc_list = .TRUE.
  !PRINT *, "After assignment alloc_list"
  !PRINT *, isAlloc_list
  

  xAll = reshape(t_x, (/nAll,np/))
  prAll = reshape(t_pr, (/nAll,nt/))
  pr2All = reshape(t_pr2, (/nAll,nt/))
  prsurvAll = reshape(t_prsurv, (/nAll,nt_death/))
  pr2survAll = reshape(t_pr2surv, (/nAll,nt_death/))
  deltaAll = t_delta
  deltaAll_m = t_delta_m

  nCat = t_nCat
  nLevs = max(maxval(nCat),1)

  mTry = t_mTry
  sampleSize = t_sampleSize !the number of cases to sample for each tree
  !sampleSize_surv = t_sampleSize_surv

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
    !PRINT *, "******************** setUpInners ********************"
    !PRINT *, "Number of cases under consideration: nAll:", nAll
    !PRINT *, "Number of cases survival:", nAll_surv
    !PRINT *, "id_RE"
    !PRINT *, id_RE
    !PRINT *, "Number of Covariates np: ", np
    !PRINT *, "t_delta"
    !PRINT *, t_delta
    !PRINT *, "t_delta_m"
    !PRINT *, t_delta_m
    !PRINT *, "******************************************************"
  END IF

  !PRINT *, "End of SetupInners"

END SUBROUTINE setUpInners
! =================================================================================

! access function for calculating forest
SUBROUTINE survTree(tFunc, mean, Prob)
  USE INNERS
  IMPLICIT NONE

  REAL(dp), DIMENSION(1:nrNodes*nt_death), INTENT(OUT) :: tFunc
  REAL(dp), DIMENSION(1:nrNodes), INTENT(OUT) :: mean
  REAL(dp), DIMENSION(1:nrNodes), INTENT(OUT) :: Prob

  !PRINT *, "******************** survTree ********************"
  CALL tsurvTree(tFunc, mean, Prob)

END SUBROUTINE survTree
! =================================================================================

! access function for calculating forest
SUBROUTINE endpointTree(tFunc, mean, Prob)
  USE INNERS
  IMPLICIT NONE

  REAL(dp), DIMENSION(1:nrNodes*nt), INTENT(OUT) :: tFunc
  REAL(dp), DIMENSION(1:nrNodes), INTENT(OUT) :: mean
  REAL(dp), DIMENSION(1:nrNodes), INTENT(OUT) :: Prob

  !PRINT *, "******************** endpointTree ********************"
  CALL tsurvTree(tFunc, mean, Prob)

END SUBROUTINE endpointTree
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