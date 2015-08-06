PROGRAM doubleDecay
! University of Birmingham
! Ben Palmer
!
!
! Force declaration of all variables
  Implicit None
! Include MPI header
  Include 'mpif.h'
! Kinds
  Integer, Parameter :: SingleReal = Selected_Real_Kind(6,37)         ! single real, 6 decimal precision, exponent range 37
  Integer, Parameter :: DoubleReal = Selected_Real_Kind(15,307)       ! double real, 15 decimal precision, exponent range 307
  Integer, Parameter :: QuadrupoleReal = Selected_Real_Kind(33,4931)  ! quadrupole real
  Integer, Parameter :: TinyInteger = Selected_Int_Kind(1)            ! tiny integer    1 byte
  Integer, Parameter :: SmallInteger = Selected_Int_Kind(4)           ! small integer    4 bytes -2E31 to 2E31-1
  Integer, Parameter :: StandardInteger = Selected_Int_Kind(8)        ! standard integer 8 bytes -2E63 to 2E63-1
  Integer, Parameter :: LongInteger = Selected_Int_Kind(12)           ! long integer
  Integer, Parameter :: VeryLongInteger = Selected_Int_Kind(32)       ! very long integer
! Variables
  Integer(kind=LongInteger) :: mMaths_randomLCG_n=0  
  Integer(kind=LongInteger) :: mMaths_randomLCG_xn   
  Integer(kind=StandardInteger) :: mMaths_processID, mMaths_processCount, mMaths_error  
  Integer(kind=StandardInteger) :: mMaths_ddfRssCount
  Integer(kind=StandardInteger) :: error
  Integer(kind=StandardInteger) :: fitType=2
! Read data
  Real(kind=DoubleReal), Dimension(1:10000,1:2) :: rawData  
  Integer(kind=StandardInteger) :: dataRows
  
! Init MPI
  Call MPI_Init(error)
! Read Data
  Call readData()  
! Run fit
  Call M_expFit(dataRows)  
! Finalise MPI  
  Call MPI_Finalize(error)  
  
  Contains
! ----------------------------------------------------------------------------------------------------------------------------   

  Subroutine readData()
    Implicit None   ! Force declaration of all variables
! Private variables  
    Integer(kind=StandardInteger) :: i
    Integer(kind=StandardInteger) :: ios
    Character(len=255) :: inputFilePath, fileRow
    Character(len=32) :: bufferA, bufferB
! Read in command line arguments
    Call get_command_argument(1,inputFilePath)
    Call get_command_argument(2,bufferA)
    If(bufferA(1:1).ne." ")Then
      Read(bufferA,*) fitType
    End If  
! Read in data
    Open(UNIT=1,FILE=trim(inputFilePath))
    dataRows = 0
    Do i=1,10000
      Read(1,"(A255)",IOSTAT=ios) fileRow
      If(ios /= 0)Then
        EXIT
      End If
      If(fileRow(1:4).ne."    ")Then
        dataRows = dataRows + 1
        Read(fileRow,*) bufferA, bufferB
        Read(bufferA,*) rawData(dataRows,1)
        Read(bufferB,*) rawData(dataRows,2)      
      End If
    End Do
  End Subroutine readData
  
  Subroutine M_expFit(dataRows)
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i, dataRows
    Real(kind=DoubleReal), Dimension(1:dataRows,1:2) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:3) :: outputA
    Real(kind=DoubleReal), Dimension(1:6) :: outputB
    Real(kind=DoubleReal), Dimension(1:8) :: outputC
    Real(kind=DoubleReal) :: programStartTime, programEndTime
! MPI vars    
    Integer(kind=StandardInteger) :: processID, processCount, error
! Start Time
    Call cpu_time(programStartTime)
! Call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)  
! Output    
    If(processID.eq.0)Then
      print "(A30)","------------------------------"
      print "(A30)"," Double Decay Exponential Fit "
      print "(A30)","------------------------------"
      print "(A10,I2)","Fit Type: ",fitType
    End If        
! Fill array
    Do i=1,dataRows
      dataPoints(i,1) = rawData(i,1)
      dataPoints(i,2) = rawData(i,2)
    End Do
    outputA = 0.0D0
    outputB = 0.0D0
    outputC = 0.0D0
! Fit - one exponential term
    If(fitType.eq.1)Then
      outputA = SingleDecayFit(dataPoints)
! End Time    
      Call cpu_time(programEndTime)
! Output    
      If(processID.eq.0)Then
        print "(A21,F14.7)","Time:                ",(programEndTime-programStartTime)
        print "(A21,E14.7)","Final RSS Value:     ",outputA(3)
        print "(A33)","f(x) = a exp(lA x)              "
        print "(A21,E14.7)","a:                   ",outputA(1)
        print "(A21,E14.7)","lA:                  ",outputA(2)    
      End If
    End If  
! Fit - two exponential terms 
    If(fitType.eq.2)Then
      outputB = M_DoubleDecayFit(dataPoints)
! End Time    
      Call cpu_time(programEndTime)
! Output    
      If(processID.eq.0)Then
        print "(A21,F14.7)","Time:                ",(programEndTime-programStartTime)
        print "(A21,I8)","RSS Calculations:    ",Ceiling(1.0D0*outputB(6))
        print "(A21,E14.7)","Final RSS Value:     ",outputB(5)
        print "(A33)","f(x) = a exp(lA x) + b exp(lB x) "
        print "(A21,E14.7)","a:                   ",outputB(1)
        print "(A21,E14.7)","lA:                  ",outputB(2)
        print "(A21,E14.7)","b:                   ",outputB(3)
        print "(A21,E14.7)","lB:                  ",outputB(4)      
      End If
    End If  
! Fit - three exponential terms 
    If(fitType.eq.3)Then
      outputC = M_TripleDecayFit(dataPoints)
! End Time    
      Call cpu_time(programEndTime)
! Output    
      If(processID.eq.0)Then
        print "(A21,F14.7)","Time:                ",(programEndTime-programStartTime)
        print "(A21,I8)","RSS Calculations:    ",Ceiling(1.0D0*outputC(8))
        print "(A21,E14.7)","Final RSS Value:     ",outputC(7)
        print "(A47)","f(x) = a exp(lA x) + b exp(lB x) + c exp(lC x) "
        print "(A21,E14.7)","a:                   ",outputC(1)
        print "(A21,E14.7)","lA:                  ",outputC(2)
        print "(A21,E14.7)","b:                   ",outputC(3)
        print "(A21,E14.7)","lB:                  ",outputC(4)   
        print "(A21,E14.7)","c:                   ",outputC(5)
        print "(A21,E14.7)","lC:                  ",outputC(6)     
      End If
    End If      
    
  End Subroutine M_expFit   
  
! ------------------------------------------------------------------------!
! Functions
! ------------------------------------------------------------------------!

  Function SingleDecayFit(dataPoints) RESULT (output)
! Fits double exponential to data
! f(x) = a exp(lA)
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, iA, m, n, maxLoops
    Real(kind=DoubleReal) :: rss, lastRSS, optRSS, convergence, maxRSSVal
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: R
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1),1:2) :: J
    Real(kind=DoubleReal), Dimension(1:2) :: parameters, parametersOpt, change
    Real(kind=DoubleReal), Dimension(1:3) :: output
    Real(kind=DoubleReal), Dimension(1:100,1:2) :: aRange
    Integer(kind=StandardInteger) :: gridA
    Real(kind=DoubleReal) :: a_T, lA_T
    Real(kind=DoubleReal), Dimension(1:100,1:3) :: topBoxes 
    Integer(kind=StandardInteger) :: topBoxCount, maxRSS 
    Logical :: storeFlag
!--------------------------------------------------
! Find Starting Parameters  
!--------------------------------------------------
! Init
    topBoxCount = 4
    topBoxes = 2.0D20
! Set a ranges
    gridA = 12
    Do i=1,gridA
      aRange(i,1) = -1.0D0*10D0**((gridA-4)-i)
      aRange(i,2) = -1.0D0*10D0**((gridA-5)-i)
    End Do    
    Do i=1,gridA
      aRange(i+gridA,1) = 1.0D0*10D0**(i-6)
      aRange(i+gridA,2) = 1.0D0*10D0**(i-5)
    End Do
    aRange(gridA,2) = 0.0D0
    aRange(gridA+1,1) = 0.0D0 
! Reduce search region into 10 smaller "boxes"
    Do iA=1,(gridA+gridA) 
      Do n=1,1000
        Do m=1,3
          a_T = (0.25D0*m)*(aRange(iA,1)+aRange(iA,2))
          lA_T = 40.0D0*(M_RandomLCG()-0.5D0)
! Calc rss
          rss = SingleDecayFitRSS(dataPoints, a_T, lA_T)  
! Check if better than boxed - update if better than already stored for this box         
          storeFlag = .true.
          Do i=1,topBoxCount
            If(topBoxes(i,2).eq.a_T.and.rss.lt.topBoxes(i,1))Then
              topBoxes(i,1) = rss
              topBoxes(i,2) = a_T
              topBoxes(i,3) = lA_T
              storeFlag = .false.
              Exit
            End If  
            If(i.eq.1)Then
              maxRSS = 1
              maxRSSVal = topBoxes(1,1)
            Else
              If(topBoxes(i,1).gt.maxRSSVal)Then
                maxRSS = i
                maxRSSVal = topBoxes(i,1)
              End If            
            End If
          End Do
! If better than any in box
          If(storeFlag)Then
            If(rss.lt.maxRSSVal)Then
              topBoxes(maxRSS,1) = rss
              topBoxes(maxRSS,2) = a_T
              topBoxes(maxRSS,3) = lA_T
            End If
          End If          
        End Do
      End Do
    End Do
!--------------------------------------------------
! Newton Gauss Elimination
!--------------------------------------------------
    Do n=1,topBoxCount
      parameters(1) = topBoxes(n,2)
      parameters(2) = topBoxes(n,3) 
! NG Opt
      lastRSS = 0
      convergence = 1.0D0
      maxLoops = 0
      Do While(maxLoops.le.100.and.convergence.gt.1.0D-7)
        maxLoops = maxLoops + 1
        rss = 0.0D0
        Do i=1,size(dataPoints,1)
          R(i) = parameters(1)*exp(parameters(2)*dataPoints(i,1))-dataPoints(i,2)   ! f(x)-y
          J(i,1) = exp(parameters(2)*dataPoints(i,1))  ! d/dx1
          J(i,2) = dataPoints(i,1)*parameters(1)*exp(parameters(2)*dataPoints(i,1))  ! d/dx2
          rss = rss + R(i)**2
        End Do
        change = NewtonGaussOpt(J,R)
        Do i=1,size(change)
          parameters(i) = parameters(i) + change(i)
        End Do
        convergence = abs(lastRSS-rss)
        lastRSS = rss
      End Do
      If(n.eq.1)Then
        optRSS = rss
        Do i=1,size(change)
          parametersOpt(i) = parameters(i)
        End Do
      Else
        If(rss.lt.optRSS)Then
          optRSS = rss
          Do i=1,size(change)
            parametersOpt(i) = parameters(i)
          End Do
        End If
      End If  
    End Do  
    output(1) = parametersOpt(1)
    output(2) = parametersOpt(2)
    output(3) = optRSS
  End Function SingleDecayFit    

  Function M_DoubleDecayFit(dataPoints,gridAIn,gridAFactorIn,&
  searchLoopsIn,filterLoopsIn) RESULT (output)
! Fits double exponential to data
! f(x) = a exp(lA) + b exp(lB)
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, n, m, r
    Integer(kind=StandardInteger) :: iA,iB
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal) :: a, b, lA, lB
    Real(kind=DoubleReal) :: a_T, b_T, lA_T, lB_T
    Real(kind=DoubleReal) :: rss, startRSS, optRSS   
    Real(kind=DoubleReal), Dimension(1:6) :: output
    Real(kind=DoubleReal), Dimension(1:100,1:2) :: lRange
    Real(kind=DoubleReal), Dimension(1:100,1:2) :: aRange
    Integer(kind=StandardInteger) :: filterLoops, searchLoops, gridA
    Integer(kind=StandardInteger), Optional :: searchLoopsIn, filterLoopsIn, gridAIn
    Real(kind=DoubleReal) :: gridAFactor
    Real(kind=DoubleReal), Optional :: gridAFactorIn
    Integer(kind=StandardInteger) :: maxRSS
    Real(kind=DoubleReal) :: maxRSSVal, variation
    Integer(kind=StandardInteger) :: topBoxCount
    Real(kind=DoubleReal), Dimension(1:100,1:5) :: topBoxes 
    Real(kind=DoubleReal), Dimension(1:5) :: topBoxesSwap 
    Logical :: sortFlag, storeFlag
! Newton Gauss Opt Vars
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: residuals
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1),1:4) :: jacobian
    Real(kind=DoubleReal), Dimension(1:4) :: parameters, parametersOpt, change
    Real(kind=DoubleReal) :: convergence, lastRSS
    Integer(kind=StandardInteger) :: maxLoops
! Call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,mMaths_processID,mMaths_error)
    Call MPI_Comm_size(MPI_COMM_WORLD,mMaths_processCount,mMaths_error)    
! User options:
! gridA - size of grid -10^(grid-6) to 10^(grid-6) [grid size = gridA^2]
! gridAFactor - span of grid increased/decreased by factor, bt keeps grid size same
! searchLoops - no. random lA and lB points for each (a,b) box
! vary loops - number of refinement loops for final result
! filter loops - number of loops to filter out top results   
!     
! Init
    mMaths_ddfRssCount = 0
! Optional   
    gridA = 10
    If(Present(gridAIn))Then
      gridA = gridAIn
    End If  
    gridAFactor = 1.0D0
    If(Present(gridAFactorIn))Then
      gridAFactor = gridAFactorIn
    End If  
    searchLoops = 500
    If(Present(searchLoopsIn))Then
      searchLoops = searchLoopsIn
    End If  
    filterLoops = 1000
    If(Present(filterLoopsIn))Then
      filterLoops = filterLoopsIn
    End If       
! Init
    topBoxCount = 15
    topBoxes = 2.0D20
! Set first test values
    a = 1.0D0
    b = 1.0D0
    lA = 1.0D0
    lB = 1.0D0
    startRSS = M_DoubleDecayFitRSS(dataPoints, a, b, lA, lB)
    optRSS = startRSS
! Assume:
! -10^5 < a < 10^5 
! -10^5 < b < 10^5 
! -20 < lA < 20 
! -20 < lB < 20 
    Do i=1,5
      lRange(i,1) = -2.0D0*10D0**(2-i)
      lRange(i,2) = -2.0D0*10D0**(1-i)
    End Do    
    Do i=1,5
      lRange(i+5,1) = 2.0D0*10D0**(i-5)
      lRange(i+5,2) = 2.0D0*10D0**(i-4)
    End Do
    lRange(5,2) = 0.0D0
    lRange(6,1) = 0.0D0    
! Set a and b ranges
    Do i=1,gridA
      aRange(i,1) = -1.0D0*gridAFactor*10D0**((gridA-4)-i)
      aRange(i,2) = -1.0D0*gridAFactor*10D0**((gridA-5)-i)
    End Do    
    Do i=1,gridA
      aRange(i+gridA,1) = 1.0D0*gridAFactor*10D0**(i-6)
      aRange(i+gridA,2) = 1.0D0*gridAFactor*10D0**(i-5)
    End Do
    aRange(gridA,2) = 0.0D0
    aRange(gridA+1,1) = 0.0D0 
! Reduce search region into 10 smaller "boxes"
    Do iA=1,(gridA+gridA) 
      Do iB=iA,(gridA+gridA)  
! Set values
        Do m=1,3
          Do n=1,ceiling(1.0D0*searchLoops/mMaths_processCount)
            a_T = (0.25D0*m)*(aRange(iA,1)+aRange(iA,2))
            lA_T = 20.0D0*(M_RandomLCG()-0.5D0)
            b_T = (0.25D0*m)*(aRange(iB,1)+aRange(iB,2))
            lB_T = 20.0D0*(M_RandomLCG()-0.5D0)
! Calculate RSS
            rss = M_DoubleDecayFitRSS(dataPoints, a_T, b_T, lA_T, lB_T)
            storeFlag = .true.
            Do i=1,topBoxCount
              If(topBoxes(i,2).eq.a_T.and.topBoxes(i,4).eq.b_T)Then
                storeFlag = .false.
                If(rss.lt.topBoxes(i,1))Then
                  topBoxes(i,1) = rss
                  topBoxes(i,2) = a_T
                  topBoxes(i,3) = lA_T
                  topBoxes(i,4) = b_T
                  topBoxes(i,5) = lB_T
                End If
                Exit
              End If
            End Do 
! Store if in best 10
            If(storeFlag)Then
              maxRSS = 1
              maxRSSVal = topBoxes(1,1)
              Do i=2,topBoxCount
                If(topBoxes(i,1).gt.maxRSSVal)Then
                  maxRSS = i
                  maxRSSVal = topBoxes(i,1)
                End If
              End Do
              If(rss.lt.maxRSSVal)Then
                topBoxes(maxRSS,1) = rss
                topBoxes(maxRSS,2) = a_T
                topBoxes(maxRSS,3) = lA_T
                topBoxes(maxRSS,4) = b_T
                topBoxes(maxRSS,5) = lB_T
              End If
            End If
          End Do
        End Do
      End Do
    End Do
! Combine results from all mpiProcesses
! Pick best boxes across all processes
    Call M_collDouble2D_Best(topBoxes, 1, topBoxCount)    
! Sort Start-----------
    sortFlag = .true.
    Do While(sortFlag)
      sortFlag = .false.
      Do i=2,topBoxCount
        If(topBoxes(i,1).lt.topBoxes(i-1,1))Then
          sortFlag = .true.
          topBoxesSwap(1) = topBoxes(i,1)
          topBoxesSwap(2) = topBoxes(i,2)
          topBoxesSwap(3) = topBoxes(i,3)
          topBoxesSwap(4) = topBoxes(i,4)
          topBoxesSwap(5) = topBoxes(i,5)
          topBoxes(i,1) = topBoxes(i-1,1)
          topBoxes(i,2) = topBoxes(i-1,2)
          topBoxes(i,3) = topBoxes(i-1,3)
          topBoxes(i,4) = topBoxes(i-1,4)
          topBoxes(i,5) = topBoxes(i-1,5)
          topBoxes(i-1,1) = topBoxesSwap(1)
          topBoxes(i-1,2) = topBoxesSwap(2)
          topBoxes(i-1,3) = topBoxesSwap(3)
          topBoxes(i-1,4) = topBoxesSwap(4)
          topBoxes(i-1,5) = topBoxesSwap(5)
        End If
      End Do
    End Do
! Sort End-----------  
! Start refinement loop
    filterLoops = ceiling(1.0D0*filterLoops/(5.0D0*mMaths_processCount))  ! 5 refinement loops
    Do r=1,5
      variation = 1.0D0/(2.0D0**(r-1))
! Filter out top boxes      
      Do n=1,topBoxCount
        optRSS = topBoxes(n,1)
        a = topBoxes(n,2)
        lA = topBoxes(n,3)
        b = topBoxes(n,4)
        lB = topBoxes(n,5)
        Do i=1,filterLoops
          a_T = topBoxes(n,2)*(1+variation*(M_RandomLCG()-0.5D0))
          lA_T = topBoxes(n,3)*(1+variation*(M_RandomLCG()-0.5D0))
          b_T = topBoxes(n,4)*(1+variation*(M_RandomLCG()-0.5D0))
          lB_T = topBoxes(n,5)*(1+variation*(M_RandomLCG()-0.5D0))
          rss = M_DoubleDecayFitRSS(dataPoints, a_T, b_T, lA_T, lB_T)
          If(rss.lt.optRSS)Then
            optRSS = rss
            a = a_T
            lA = lA_T
            b = b_T
            lB = lB_T
          End If
        End Do 
! Update Values
        topBoxes(n,1) = optRSS
        topBoxes(n,2) = a
        topBoxes(n,3) = lA
        topBoxes(n,4) = b
        topBoxes(n,5) = lB
      End Do 
! Sort Start-----------
      sortFlag = .true.
      Do While(sortFlag)
        sortFlag = .false.
        Do i=2,topBoxCount
          If(topBoxes(i,1).lt.topBoxes(i-1,1))Then
            sortFlag = .true.
            topBoxesSwap(1) = topBoxes(i,1)
            topBoxesSwap(2) = topBoxes(i,2)
            topBoxesSwap(3) = topBoxes(i,3)
            topBoxesSwap(4) = topBoxes(i,4)
            topBoxesSwap(5) = topBoxes(i,5)
            topBoxes(i,1) = topBoxes(i-1,1)
            topBoxes(i,2) = topBoxes(i-1,2)
            topBoxes(i,3) = topBoxes(i-1,3)
            topBoxes(i,4) = topBoxes(i-1,4)
            topBoxes(i,5) = topBoxes(i-1,5)
            topBoxes(i-1,1) = topBoxesSwap(1)
            topBoxes(i-1,2) = topBoxesSwap(2)
            topBoxes(i-1,3) = topBoxesSwap(3)
            topBoxes(i-1,4) = topBoxesSwap(4)
            topBoxes(i-1,5) = topBoxesSwap(5)
          End If
        End Do
      End Do
! Sort End-----------      
      ! Reduce number of top boxes by 25%
      topBoxCount = ceiling(0.75D0*topBoxCount)
      If(topBoxCount.lt.2)Then
        topBoxCount = 2
      End If
    End Do   
! Store results
    optRSS = topBoxes(1,1)
    a = topBoxes(1,2)
    lA = topBoxes(1,3)
    b = topBoxes(1,4)
    lB = topBoxes(1,5)  
! Best result from all processes
    Call M_bestParameters(optRSS,a,lA,b,lB)

! Newton Gauss on root process
    If(mMaths_processID.eq.0)Then    

!--------------------------------------------------
! Newton Gauss Elimination
!--------------------------------------------------  
      parameters = 0.0D0 
      parametersOpt = 0.0D0
! set parameters      
      parameters(1) = a
      parameters(2) = lA
      parameters(3) = b
      parameters(4) = lB
! NG Opt
      lastRSS = 0
      convergence = 1.0D0
      maxLoops = 0
      Do While(maxLoops.le.100.and.convergence.gt.1.0D-7)
        maxLoops = maxLoops + 1
        rss = 0.0D0
        Do i=1,size(dataPoints,1)
          residuals(i) = (parameters(1)*exp(parameters(2)*dataPoints(i,1))&
          +parameters(3)*exp(parameters(4)*dataPoints(i,1)))&
          -dataPoints(i,2)   ! f(x)-y
          jacobian(i,1) = exp(parameters(2)*dataPoints(i,1))  ! d/dx1
          jacobian(i,2) = dataPoints(i,1)*parameters(1)*exp(parameters(2)*dataPoints(i,1))  ! d/dx2
          jacobian(i,3) = exp(parameters(4)*dataPoints(i,1))  ! d/dx1
          jacobian(i,4) = dataPoints(i,1)*parameters(3)*exp(parameters(4)*dataPoints(i,1))  ! d/dx2
          rss = rss + residuals(i)**2          
        End Do
        change = NewtonGaussOpt(jacobian,residuals)
        Do i=1,size(change)
          parameters(i) = parameters(i) + change(i)
        End Do
        convergence = abs(lastRSS-rss)
        lastRSS = rss
      End Do
      If(n.eq.1)Then
        optRSS = rss
        Do i=1,size(change)
          parametersOpt(i) = parameters(i)
        End Do
      Else
        If(rss.lt.optRSS)Then
          optRSS = rss
          Do i=1,size(change)
            parametersOpt(i) = parameters(i)
          End Do
        End If
      End If 
      a = parametersOpt(1)
      lA = parametersOpt(2)
      b = parametersOpt(3)
      lB = parametersOpt(4)
    End If
! Best result from all processes
    Call M_bestParameters(optRSS,a,lA,b,lB)        
! Output
    output(1) = a
    output(2) = lA
    output(3) = b
    output(4) = lB
    output(5) = optRSS
    output(6) = mMaths_ddfRssCount*mMaths_processCount
  End Function M_DoubleDecayFit   
  
  
  
  Function M_TripleDecayFit(dataPoints,gridAIn,gridAFactorIn,&
  searchLoopsIn,filterLoopsIn) RESULT (output)
! Fits double exponential to data
! f(x) = a exp(lA) + b exp(lB)
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, n, m, r
    Integer(kind=StandardInteger) :: iA,iB,iC
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal) :: a, b, c, lA, lB, lC
    Real(kind=DoubleReal) :: a_T, b_T, c_T, lA_T, lB_T, lC_T
    Real(kind=DoubleReal) :: rss, startRSS, optRSS   
    Real(kind=DoubleReal), Dimension(1:8) :: output
    Real(kind=DoubleReal), Dimension(1:100,1:2) :: aRange
    Integer(kind=StandardInteger) :: filterLoops, searchLoops, gridA
    Integer(kind=StandardInteger), Optional :: searchLoopsIn, filterLoopsIn, gridAIn
    Real(kind=DoubleReal) :: gridAFactor
    Real(kind=DoubleReal), Optional :: gridAFactorIn
    Integer(kind=StandardInteger) :: maxRSS
    Real(kind=DoubleReal) :: maxRSSVal, variation
    Integer(kind=StandardInteger) :: topBoxCount
    Real(kind=DoubleReal), Dimension(1:100,1:7) :: topBoxes 
    Real(kind=DoubleReal), Dimension(1:7) :: topBoxesSwap 
    Logical :: sortFlag, storeFlag
! Newton Gauss Opt Vars
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: residuals
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1),1:6) :: jacobian
    Real(kind=DoubleReal), Dimension(1:6) :: parameters, parametersOpt, change
    Real(kind=DoubleReal) :: convergence, lastRSS
    Integer(kind=StandardInteger) :: maxLoops
! Call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,mMaths_processID,mMaths_error)
    Call MPI_Comm_size(MPI_COMM_WORLD,mMaths_processCount,mMaths_error)    
! User options:
! gridA - size of grid -10^(grid-6) to 10^(grid-6) [grid size = gridA^2]
! gridAFactor - span of grid increased/decreased by factor, bt keeps grid size same
! searchLoops - no. random lA and lB points for each (a,b) box
! vary loops - number of refinement loops for final result
! filter loops - number of loops to filter out top results   
!     
! Init
    mMaths_ddfRssCount = 0
! Optional   
    gridA = 10
    If(Present(gridAIn))Then
      gridA = gridAIn
    End If  
    gridAFactor = 1.0D0
    If(Present(gridAFactorIn))Then
      gridAFactor = gridAFactorIn
    End If  
    searchLoops = 5000
    If(Present(searchLoopsIn))Then
      searchLoops = searchLoopsIn
    End If  
    filterLoops = 20000
    If(Present(filterLoopsIn))Then
      filterLoops = filterLoopsIn
    End If       
! Init
    topBoxCount = 30
    topBoxes = 2.0D20
! Set first test values
    a = 1.0D0
    b = 1.0D0
    c = 1.0D0
    lA = 1.0D0
    lB = 1.0D0
    lC = 1.0D0
    startRSS = M_TripleDecayFitRSS(dataPoints, a, b, c, lA, lB, lC)
    optRSS = startRSS
! Assume:
! -10^5 < a < 10^5 
! -10^5 < b < 10^5 
! -10^5 < c < 10^5 
! Set a, b and c ranges
    Do i=1,gridA
      aRange(i,1) = -1.0D0*gridAFactor*10D0**((gridA-4)-i)
      aRange(i,2) = -1.0D0*gridAFactor*10D0**((gridA-5)-i)
    End Do    
    Do i=1,gridA
      aRange(i+gridA,1) = 1.0D0*gridAFactor*10D0**(i-6)
      aRange(i+gridA,2) = 1.0D0*gridAFactor*10D0**(i-5)
    End Do
    aRange(gridA,2) = 0.0D0
    aRange(gridA+1,1) = 0.0D0 
! Reduce search region into 10 smaller "boxes"
    Do iA=1,(gridA+gridA) 
      Do iB=iA,(gridA+gridA)  
        Do iC=iB,(gridA+gridA)  
! Set values
          Do m=1,3
            Do n=1,ceiling(1.0D0*searchLoops/mMaths_processCount)
              a_T = (0.25D0*m)*(aRange(iA,1)+aRange(iA,2))
              lA_T = 20.0D0*(M_RandomLCG()-0.5D0)
              b_T = (0.25D0*m)*(aRange(iB,1)+aRange(iB,2))
              lB_T = 20.0D0*(M_RandomLCG()-0.5D0)
              c_T = (0.25D0*m)*(aRange(iC,1)+aRange(iC,2))
              lC_T = 20.0D0*(M_RandomLCG()-0.5D0)
! Calculate RSS
              rss = M_TripleDecayFitRSS(dataPoints, a_T, b_T, c_T, lA_T, lB_T, lC_T)
              storeFlag = .true.
              Do i=1,topBoxCount
                If(topBoxes(i,2).eq.a_T.and.topBoxes(i,4).eq.b_T.and.topBoxes(i,6).eq.c_T)Then
                  storeFlag = .false.
                  If(rss.lt.topBoxes(i,1))Then
                    topBoxes(i,1) = rss
                    topBoxes(i,2) = a_T
                    topBoxes(i,3) = lA_T
                    topBoxes(i,4) = b_T
                    topBoxes(i,5) = lB_T
                    topBoxes(i,6) = c_T
                    topBoxes(i,7) = lC_T
                  End If
                  Exit
                End If
              End Do 
! Store if in best 10
              If(storeFlag)Then
                maxRSS = 1
                maxRSSVal = topBoxes(1,1)
                Do i=2,topBoxCount
                  If(topBoxes(i,1).gt.maxRSSVal)Then
                    maxRSS = i
                    maxRSSVal = topBoxes(i,1)
                  End If
                End Do
                If(rss.lt.maxRSSVal)Then
                  topBoxes(maxRSS,1) = rss
                  topBoxes(maxRSS,2) = a_T
                  topBoxes(maxRSS,3) = lA_T
                  topBoxes(maxRSS,4) = b_T
                  topBoxes(maxRSS,5) = lB_T
                  topBoxes(maxRSS,6) = c_T
                  topBoxes(maxRSS,7) = lC_T
                End If
              End If
            End Do
          End Do  
        End Do
      End Do
    End Do
    Do i=1,topBoxCount
      print *,i,topBoxes(i,1),topBoxes(i,2),topBoxes(i,3),&
      topBoxes(i,4),topBoxes(i,5),topBoxes(i,6),topBoxes(i,7)
    End Do
! Combine results from all mpiProcesses
! Pick best boxes across all processes
    Call M_collDouble2D_Best(topBoxes, 1, topBoxCount)    
! Sort Start-----------
    sortFlag = .true.
    Do While(sortFlag)
      sortFlag = .false.
      Do i=2,topBoxCount
        If(topBoxes(i,1).lt.topBoxes(i-1,1))Then
          sortFlag = .true.
          topBoxesSwap(1) = topBoxes(i,1)
          topBoxesSwap(2) = topBoxes(i,2)
          topBoxesSwap(3) = topBoxes(i,3)
          topBoxesSwap(4) = topBoxes(i,4)
          topBoxesSwap(5) = topBoxes(i,5)
          topBoxesSwap(6) = topBoxes(i,6)
          topBoxesSwap(7) = topBoxes(i,7)
          topBoxes(i,1) = topBoxes(i-1,1)
          topBoxes(i,2) = topBoxes(i-1,2)
          topBoxes(i,3) = topBoxes(i-1,3)
          topBoxes(i,4) = topBoxes(i-1,4)
          topBoxes(i,5) = topBoxes(i-1,5)
          topBoxes(i,6) = topBoxes(i-1,6)
          topBoxes(i,7) = topBoxes(i-1,7)
          topBoxes(i-1,1) = topBoxesSwap(1)
          topBoxes(i-1,2) = topBoxesSwap(2)
          topBoxes(i-1,3) = topBoxesSwap(3)
          topBoxes(i-1,4) = topBoxesSwap(4)
          topBoxes(i-1,5) = topBoxesSwap(5)
          topBoxes(i-1,6) = topBoxesSwap(6)
          topBoxes(i-1,7) = topBoxesSwap(7)
        End If
      End Do
    End Do
! Sort End-----------  
! Start refinement loop
    filterLoops = ceiling(1.0D0*filterLoops/(5.0D0*mMaths_processCount))  ! 5 refinement loops
    Do r=1,5
      variation = 1.0D0/(2.0D0**(r-1))
! Filter out top boxes      
      Do n=1,topBoxCount
        optRSS = topBoxes(n,1)
        a = topBoxes(n,2)
        lA = topBoxes(n,3)
        b = topBoxes(n,4)
        lB = topBoxes(n,5)
        c = topBoxes(n,6)
        lC = topBoxes(n,7)
        Do i=1,filterLoops
          a_T = topBoxes(n,2)*(1+variation*(M_RandomLCG()-0.5D0))
          lA_T = topBoxes(n,3)*(1+variation*(M_RandomLCG()-0.5D0))
          b_T = topBoxes(n,4)*(1+variation*(M_RandomLCG()-0.5D0))
          lB_T = topBoxes(n,5)*(1+variation*(M_RandomLCG()-0.5D0))
          c_T = topBoxes(n,6)*(1+variation*(M_RandomLCG()-0.5D0))
          lC_T = topBoxes(n,7)*(1+variation*(M_RandomLCG()-0.5D0))
          rss = M_TripleDecayFitRSS(dataPoints, a_T, b_T, c_T, lA_T, lB_T, lC_T)
          If(rss.lt.optRSS)Then
            optRSS = rss
            a = a_T
            lA = lA_T
            b = b_T
            lB = lB_T
            c = c_T
            lC = lC_T
          End If
        End Do 
! Update Values
        topBoxes(n,1) = optRSS
        topBoxes(n,2) = a
        topBoxes(n,3) = lA
        topBoxes(n,4) = b
        topBoxes(n,5) = lB
        topBoxes(n,6) = c
        topBoxes(n,7) = lC
      End Do 
! Sort Start-----------
      sortFlag = .true.
      Do While(sortFlag)
        sortFlag = .false.
        Do i=2,topBoxCount
          If(topBoxes(i,1).lt.topBoxes(i-1,1))Then
            sortFlag = .true.
            topBoxesSwap(1) = topBoxes(i,1)
            topBoxesSwap(2) = topBoxes(i,2)
            topBoxesSwap(3) = topBoxes(i,3)
            topBoxesSwap(4) = topBoxes(i,4)
            topBoxesSwap(5) = topBoxes(i,5)
            topBoxesSwap(6) = topBoxes(i,6)
            topBoxesSwap(7) = topBoxes(i,7)
            topBoxes(i,1) = topBoxes(i-1,1)
            topBoxes(i,2) = topBoxes(i-1,2)
            topBoxes(i,3) = topBoxes(i-1,3)
            topBoxes(i,4) = topBoxes(i-1,4)
            topBoxes(i,5) = topBoxes(i-1,5)
            topBoxes(i,6) = topBoxes(i-1,6)
            topBoxes(i,7) = topBoxes(i-1,7)
            topBoxes(i-1,1) = topBoxesSwap(1)
            topBoxes(i-1,2) = topBoxesSwap(2)
            topBoxes(i-1,3) = topBoxesSwap(3)
            topBoxes(i-1,4) = topBoxesSwap(4)
            topBoxes(i-1,5) = topBoxesSwap(5)
            topBoxes(i-1,6) = topBoxesSwap(6)
            topBoxes(i-1,7) = topBoxesSwap(7)
          End If
        End Do
      End Do
! Sort End-----------      
      ! Reduce number of top boxes by 25%
      topBoxCount = ceiling(0.75D0*topBoxCount)
      If(topBoxCount.lt.2)Then
        topBoxCount = 2
      End If
    End Do   
    Do i=1,1
      print *,i,topBoxes(i,1),topBoxes(i,2),topBoxes(i,3),&
      topBoxes(i,4),topBoxes(i,5),topBoxes(i,6),topBoxes(i,7)
    End Do
! Store results
    optRSS = topBoxes(1,1)
    a = topBoxes(1,2)
    lA = topBoxes(1,3)
    b = topBoxes(1,4)
    lB = topBoxes(1,5)  
    c = topBoxes(1,6)
    lC = topBoxes(1,7) 
! Best result from all processes
    Call M_bestParametersB(optRSS,a,lA,b,lB,c,lC)
    print *,optRSS,a,lA,b,lB,c,lC
! Newton Gauss on root process
    If(mMaths_processID.eq.0)Then   
!--------------------------------------------------
! Newton Gauss Elimination
!--------------------------------------------------  
      parameters = 0.0D0 
      parametersOpt = 0.0D0
! set parameters      
      parameters(1) = a
      parameters(2) = lA
      parameters(3) = b
      parameters(4) = lB
      parameters(5) = c
      parameters(6) = lC
! NG Opt
      lastRSS = 0
      convergence = 1.0D0
      maxLoops = 0
      Do While(maxLoops.le.100.and.convergence.gt.1.0D-7)
        maxLoops = maxLoops + 1
        rss = 0.0D0
        Do i=1,size(dataPoints,1)
          residuals(i) = &
          (parameters(1)*exp(parameters(2)*dataPoints(i,1))&
          +parameters(3)*exp(parameters(4)*dataPoints(i,1))&
          +parameters(5)*exp(parameters(6)*dataPoints(i,1)))&
          -dataPoints(i,2)   ! f(x)-y
          jacobian(i,1) = exp(parameters(2)*dataPoints(i,1))  ! d/dx1
          jacobian(i,2) = dataPoints(i,1)*parameters(1)*exp(parameters(2)*dataPoints(i,1))  ! d/dx2
          jacobian(i,3) = exp(parameters(4)*dataPoints(i,1))  ! d/dx3
          jacobian(i,4) = dataPoints(i,1)*parameters(3)*exp(parameters(4)*dataPoints(i,1))  ! d/dx4     
          jacobian(i,5) = exp(parameters(6)*dataPoints(i,1))  ! d/dx5
          jacobian(i,6) = dataPoints(i,1)*parameters(5)*exp(parameters(6)*dataPoints(i,1))  ! d/dx6
          rss = rss + residuals(i)**2    
        End Do
        change = NewtonGaussOpt(jacobian,residuals)
        print *,""
        Do i=1,size(change)
          parameters(i) = parameters(i) + change(i)
          print *,change(i)
        End Do
        print *,""
        convergence = abs(lastRSS-rss)
        lastRSS = rss
      End Do
      If(n.eq.1)Then
        optRSS = rss
        Do i=1,size(change)
          parametersOpt(i) = parameters(i)
        End Do
      Else
        If(rss.lt.optRSS)Then
          optRSS = rss
          Do i=1,size(change)
            parametersOpt(i) = parameters(i)
          End Do
        End If
      End If 
      a = parametersOpt(1)
      lA = parametersOpt(2)
      b = parametersOpt(3)
      lB = parametersOpt(4)
      c = parametersOpt(5)
      lC = parametersOpt(6)
    End If
! Best result from all processes
    Call M_bestParametersB(optRSS,a,lA,b,lB,c,lC)        
! Output
    output(1) = a
    output(2) = lA
    output(3) = b
    output(4) = lB
    output(5) = c
    output(6) = lC
    output(7) = optRSS
    output(8) = mMaths_ddfRssCount*mMaths_processCount
  End Function M_TripleDecayFit   
  
    
  Function SingleDecayFitRSS(dataPoints, a, lA) RESULT (rss)
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal) :: a, lA, rss, x, y
    rss = 0.0D0
    Do i=1,size(dataPoints,1)
      x = dataPoints(i,1)
      y = a*exp(lA*x)
      rss = rss + (dataPoints(i,2)-y)**2
    End Do
  End Function SingleDecayFitRSS 
  
  Function M_DoubleDecayFitRSS(dataPoints, a, b, lA, lB) RESULT (rss)
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal) :: a, b, lA, lB, rss, x, y
    rss = 0.0D0
    Do i=1,size(dataPoints,1)
      x = dataPoints(i,1)
      y = a*exp(lA*x)+b*exp(lB*x)
      rss = rss + (dataPoints(i,2)-y)**2
    End Do
    mMaths_ddfRssCount = mMaths_ddfRssCount + 1
  End Function M_DoubleDecayFitRSS   
  
  Function M_TripleDecayFitRSS(dataPoints, a, b, c, lA, lB, lC) RESULT (rss)
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal) :: a, b, c, lA, lB, lC, rss, x, y
    rss = 0.0D0
    Do i=1,size(dataPoints,1)
      x = dataPoints(i,1)
      y = a*exp(lA*x)+b*exp(lB*x)+c*exp(lC*x)
      rss = rss + (dataPoints(i,2)-y)**2
    End Do
    mMaths_ddfRssCount = mMaths_ddfRssCount + 1
  End Function M_TripleDecayFitRSS   
  
! ------------------------------------------------------------------------!
! Random Number Related Functions
! ------------------------------------------------------------------------!  
    
  Function M_RandomLCG(seedIn) RESULT (output) 
! Random number - linear congruential generator
! Different number on each mpi process
    Implicit None ! Force declaration of all variables
! Declare variables  
    Integer(kind=LongInteger) :: m, a, c
    Integer(kind=StandardInteger) :: seed    
    Integer(kind=StandardInteger), Optional :: seedIn    
    Real(kind=DoubleReal) :: output
! init
    seed = 0
    output = 0.0D0
    m = 4294967296_LongInteger
    a = 1103515245_LongInteger
    c = 12345_LongInteger
! Read input, reset counter    
    If(Present(seedIn))Then
      seed = seedIn
      mMaths_randomLCG_n = 0
    End If
! If first iteration
    If(mMaths_randomLCG_n.eq.0)Then    
      If(seed.eq.0)Then
        seed = 12791244 ! Use default seed
      End If  
! Call mpi subroutines
      Call MPI_Comm_rank(MPI_COMM_WORLD,mMaths_processID,mMaths_error)
      Call MPI_Comm_size(MPI_COMM_WORLD,mMaths_processCount,mMaths_error)  
! Use default seed
      seed = mod(seed+(mMaths_processID*12791244),733515237)
! Set different seed for each mpi process      
      mMaths_randomLCG_n = 0
      mMaths_randomLCG_xn = seed
! Apply a couple of times to stop 1st random number being too similar on all processes     
      mMaths_randomLCG_xn = mod((a*mMaths_randomLCG_xn+c),m)
      mMaths_randomLCG_xn = mod((a*mMaths_randomLCG_xn+c),m)
      mMaths_randomLCG_xn = mod((a*mMaths_randomLCG_xn+c),m)
    End If
! Increment counter    
    mMaths_randomLCG_n = mMaths_randomLCG_n + 1
! calculate
    mMaths_randomLCG_xn = mod((a*mMaths_randomLCG_xn+c),m)
    output = (1.0D0*mMaths_randomLCG_xn)/(1.0D0*m)
  End Function M_RandomLCG

! ------------------------------------------------------------------------!
! Newton Gauss
! ------------------------------------------------------------------------!  

  Function NewtonGaussOpt(J,R) RESULT (P)!   
    Implicit None  !Force declaration of all variables
! Declare variables 
    Real(kind=DoubleReal), Dimension(:,:) :: J   ! Jacobian
    Real(kind=DoubleReal), Dimension(:) :: R      ! Residuals
    Real(kind=DoubleReal), Dimension(1:size(J,2),1:size(J,1)) :: JT     ! Transpose Jacobian
    Real(kind=DoubleReal), Dimension(1:size(J,2),1:size(J,2)) :: JTJ    ! (Jacobian Transpose * Jacobian)
    Real(kind=DoubleReal), Dimension(1:size(J,2)) :: JTR                ! (Jacobian Transpose * Residuals)
    Real(kind=DoubleReal), Dimension(1:size(J,2)) :: P      ! Change
!***********     
! P = (JTJ)^(-1)(-1*JTR)   
!***********      
! Transpose Jacobian
    JT = TransposeMatrix(J)
    JTJ = matmul(JT,J)
    JTJ = InvertMatrix(JTJ) ! store inverse (recycle JTJ var)
    JTR = matmul(JT,R)
    JTR = -1.0D0*JTR ! Recycle JTR var
    P = matmul(JTJ,JTR)  
  End Function NewtonGaussOpt
  
  
! ------------------------------------------------------------------------!
! Matrix Functions
! ------------------------------------------------------------------------!   
  
  Function InvertMatrix(xMatrix) RESULT (xMatrixInverse)
! Invert square matrix
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: xMatrix
    Integer(kind=StandardInteger) :: row,col,rowb
    Integer(kind=StandardInteger) :: matrixSize
    Real(kind=DoubleReal), Dimension(1:size(xMatrix,1),1:2*size(xMatrix,1)) :: xMatrixWorking
    Real(kind=DoubleReal), Dimension(1:size(xMatrix,1),1:size(xMatrix,1)) :: xMatrixInverse
    Real(kind=DoubleReal), Dimension(1:2*size(xMatrix,1)) :: xMatrixRow
! matrix(row,column)
! Initialise variables
    row = 0
    rowb = 0
    col = 0
    matrixSize = size(xMatrix,1)
    xMatrixWorking = 0.0D0
    xMatrixInverse = 0.0D0
    xMatrixRow = 0.0D0
! if a square matrix
    If(size(xMatrix,1).eq.size(xMatrix,2))Then
! Fill working array
      Do row=1,matrixSize
        Do col=1,matrixSize
          xMatrixWorking(row,col) = 1.0D0*xMatrix(row,col)
        End Do
      End Do
      Do row=1,matrixSize
        Do col=1,matrixSize
          If(row.eq.col)Then
            xMatrixWorking(row,col+matrixSize) = 1.0D0
          End If
        End Do
      End Do
! make lower triangle of zeros
      Do row=1,matrixSize-1
        Do rowb=row+1,matrixSize
          If(xMatrixWorking(rowb,row).ne.0.0D0)Then !Only do if necessary
            Do col=1,(2*matrixSize) !loop over all columns
              xMatrixRow(col) = 1.0D0*&
              ((1.0D0*xMatrixWorking(row,row))/(1.0D0*xMatrixWorking(rowb,row)))*&
              xMatrixWorking(rowb,col)-1.0D0*xMatrixWorking(row,col)
            End Do
! replace row values
            Do col=1,(2*matrixSize) !loop over all columns
              xMatrixWorking(rowb,col) = 1.0D0 * xMatrixRow(col)
            End Do
          End If
        End Do
! force zeros in the lower triangle
        Do rowb=row+1,matrixSize
          xMatrixWorking(rowb,row) = 0.0D0
        End Do
      End Do
! re-force zeros in the lower triangle
      Do row=1,matrixSize
        Do col=1,matrixSize
          If(row.gt.col)Then
            xMatrixWorking(row,col) = 0.0D0
          End If
        End Do
      End Do
! make upper triangle of zeros
      Do row=matrixSize,2,-1
        Do rowb=row-1,1,-1
          If(xMatrixWorking(rowb,row).ne.0.0D0)Then !Only do if necessary
            Do col=1,(2*matrixSize) !loop over all columns
              xMatrixRow(col) = 1.0D0*&
              ((1.0D0*xMatrixWorking(row,row))/(1.0D0*xMatrixWorking(rowb,row)))*&
              xMatrixWorking(rowb,col)-1.0D0*xMatrixWorking(row,col)
            End Do
! replace row values
            Do col=1,(2*matrixSize) !loop over all columns
              xMatrixWorking(rowb,col) = 1.0D0 * xMatrixRow(col)
            End Do
          End If
        End Do
! force zeros in the upper triangle
        Do rowb=row-1,1,-1
          xMatrixWorking(rowb,row) = 0.0D0
        End Do
      End Do
! Divide rhs by diagonal on lhs and store in inverse
      Do row=1,matrixSize
        Do col=1,matrixSize
          xMatrixInverse(row,col) = 1.0D0*&
          xMatrixWorking(row,col+matrixSize)/xMatrixWorking(row,row)
        End Do
      End Do
    End If
  End Function InvertMatrix
  
  Function TransposeMatrix(xMatrix) RESULT (xMatrixTranspose)
! Invert square matrix
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: xMatrix
    Real(kind=DoubleReal), Dimension(1:size(xMatrix,2),1:size(xMatrix,1)) :: xMatrixTranspose
    Integer(kind=StandardInteger) :: row,col
! Transpose
    Do row=1,size(xMatrix,1)
      Do col=1,size(xMatrix,2)
        xMatrixTranspose(col,row) = xMatrix(row,col)
      End Do
    End Do    
  End Function TransposeMatrix
  

! ------------------------------------------------------------------------!
! Subroutines
! ------------------------------------------------------------------------!   
  
  Subroutine M_collDouble2D_Best(sendIn, startKeyIn, endKeyIn)
! Specific to the exp decay fit
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger), optional :: startKeyIn, endKeyIn
    Integer(kind=StandardInteger) :: startKey, endKey
    Integer(kind=StandardInteger) :: arrayX, arrayY
    Integer(kind=StandardInteger) :: maxRow
    Real(kind=DoubleReal), Dimension(:,:) :: sendIn
    Real(kind=DoubleReal), Dimension(1:size(sendIn,2)) :: send, receive
    Real(kind=DoubleReal) :: maxRowVal
    Integer(kind=StandardInteger) :: n, i, j, k, processTo, processFrom, tag
    Integer(kind=StandardInteger) :: processID,processCount,error
    Integer, Dimension(MPI_STATUS_SIZE) :: status
! Init variables
    maxRow = -1
    maxRowVal = 1.0D20
    n = 0
    send = 0.0D0
    receive = 0.0D0
    arrayX = size(sendIn,1)
    arrayY = size(sendIn,2)    
! Optional variables
    startKey = 1
    endKey = arrayX
    If(present(startKeyIn))Then
      startKey = startKeyIn
    End If
    If(present(endKeyIn))Then
      endKey = endKeyIn
    End If
! Call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
    If(processCount.gt.1)Then
! All processes back to root
      Do i=startKey,endKey ! loop through selected rows and put in a 1D array
! Send from workers
        If(processID.gt.0)Then  ! Worker processes
! Prep send array
          Do j=1,arrayY  ! Loop through columns (j) in array
            send(j) = sendIn(i,j) 
          End Do
          processTo = 0 ! Root process
          tag = 1000 + processID  ! process id - different for each MPI process
          Call MPI_SEND(send,arrayY,MPI_DOUBLE_PRECISION,processTo,tag,&
          MPI_COMM_WORLD,status,error)  ! Send to root
        End If
! Collect from workers by root
        If(processID.eq.0)Then   ! If root process 
          Do n=1,processCount-1  ! Loop through incoming worker processes
            processFrom = n
            tag = 1000 + n
            Call MPI_RECV(receive,arrayY,MPI_DOUBLE_PRECISION,processFrom,tag,&
            MPI_COMM_WORLD,status,error)  ! receive array
! Replace lowest, if better
            Do k=startKey,endKey
              If(k.eq.startKey)Then
                maxRow = 1
                maxRowVal = sendIn(k,1)
              Else
                If(sendIn(k,1).gt.maxRowVal)Then
                  maxRow = k
                  maxRowVal = sendIn(k,1)
                End If
              End If
            End Do
            If(receive(1).lt.maxRowVal)Then
              Do j=1,arrayY
                sendIn(maxRow,j) = receive(j)
              End Do
            End If  
          End Do
        End If
      End Do
      Call MPI_Barrier(MPI_COMM_WORLD,error)
! Send out to all processes
      Do i=startKey,endKey ! loop through selected rows and put in a 1D array
! Send from root
        If(processID.eq.0)Then  ! Root processes
! Prep send array
          Do j=1,arrayY  ! Loop through columns (j) in array
            send(j) = sendIn(i,j) 
          End Do    
          Do n=1,processCount-1  ! Loop through incoming worker processes
            processTo = n
            tag = 3000 + n
            Call MPI_SEND(send,arrayY,MPI_DOUBLE_PRECISION,processTo,tag,&
            MPI_COMM_WORLD,error)
          End Do
        End If       
! Collect from root by workers   
        If(processID.gt.0)Then  ! Worker processes 
          processFrom = 0
          tag = 3000 + processID
          Call MPI_RECV(receive,arrayY,MPI_DOUBLE_PRECISION,processFrom,tag,&
          MPI_COMM_WORLD,status,error)
          Do j=1,arrayY
            sendIn(i,j) = receive(j)
          End Do
        End If
      End Do
      Call MPI_Barrier(MPI_COMM_WORLD,error)
    End If
  End Subroutine M_collDouble2D_Best
  
  
  Subroutine M_bestParameters(rss, a, lA, b, lB)
! Specific to the exp decay fit
    Implicit None   ! Force declaration of all variables
! Private variables
    Real(kind=DoubleReal) :: rss, a, lA, b, lB
    Real(kind=DoubleReal), Dimension(1:5) :: sendArr, recvArr
    Integer(kind=StandardInteger) :: n, i, processTo, processFrom, tag
    Integer(kind=StandardInteger) :: processID,processCount,error
    Integer, Dimension(MPI_STATUS_SIZE) :: status
! Prepare send array    
    sendArr(1) = rss
    sendArr(2) = a
    sendArr(3) = lA
    sendArr(4) = b
    sendArr(5) = lB
! Call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
    If(processCount.gt.1)Then
! Collect to root process and pick best
      If(processID.gt.0)Then
        processTo = 0
        tag = 7000 + processID
        Call MPI_SEND(sendArr,5,MPI_DOUBLE_PRECISION,processTo,tag,&
        MPI_COMM_WORLD,status,error)
      End If
! Collect from workers by root
      If(processID.eq.0)Then
        Do n=1,processCount-1
          processFrom = n
          tag = 7000 + n
          Call MPI_RECV(recvArr,5,MPI_DOUBLE_PRECISION,processFrom,tag,&
          MPI_COMM_WORLD,status,error)
! If better settings, use them
          If(recvArr(1).lt.sendArr(1))Then
            Do i=1,5
              sendArr(i) = recvArr(i)
            End Do
          End If
        End Do
      End If
      Call MPI_Barrier(MPI_COMM_WORLD,error)
! Distribute to workers
      If(processID.eq.0)Then
        Do n=1,processCount-1
          processTo = n
          tag = 7000 + n
          Call MPI_SEND(sendArr,5,MPI_DOUBLE_PRECISION,processTo,tag,&
          MPI_COMM_WORLD,status,error)        
        End Do
      End If
      If(processID.gt.0)Then
        processFrom = 0
        tag = 7000 + processID
        Call MPI_RECV(recvArr,5,MPI_DOUBLE_PRECISION,processFrom,tag,&
        MPI_COMM_WORLD,status,error)
        sendArr = recvArr
      End If
      Call MPI_Barrier(MPI_COMM_WORLD,error)
    End If    
! Reassign values
    rss = sendArr(1)
    a = sendArr(2)
    lA = sendArr(3)
    b = sendArr(4)
    lB = sendArr(5)
  End Subroutine M_bestParameters
  
    Subroutine M_bestParametersB(rss, a, lA, b, lB, c, lC)
! Specific to the exp decay fit
    Implicit None   ! Force declaration of all variables
! Private variables
    Real(kind=DoubleReal) :: rss, a, lA, b, lB, c, lC
    Real(kind=DoubleReal), Dimension(1:7) :: sendArr, recvArr
    Integer(kind=StandardInteger) :: n, i, processTo, processFrom, tag
    Integer(kind=StandardInteger) :: processID,processCount,error
    Integer, Dimension(MPI_STATUS_SIZE) :: status
! Prepare send array    
    sendArr(1) = rss
    sendArr(2) = a
    sendArr(3) = lA
    sendArr(4) = b
    sendArr(5) = lB
    sendArr(6) = c
    sendArr(7) = lC
! Call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
    If(processCount.gt.1)Then
! Collect to root process and pick best
      If(processID.gt.0)Then
        processTo = 0
        tag = 7000 + processID
        Call MPI_SEND(sendArr,7,MPI_DOUBLE_PRECISION,processTo,tag,&
        MPI_COMM_WORLD,status,error)
      End If
! Collect from workers by root
      If(processID.eq.0)Then
        Do n=1,processCount-1
          processFrom = n
          tag = 7000 + n
          Call MPI_RECV(recvArr,7,MPI_DOUBLE_PRECISION,processFrom,tag,&
          MPI_COMM_WORLD,status,error)
! If better settings, use them
          If(recvArr(1).lt.sendArr(1))Then
            Do i=1,7
              sendArr(i) = recvArr(i)
            End Do
          End If
        End Do
      End If
      Call MPI_Barrier(MPI_COMM_WORLD,error)
! Distribute to workers
      If(processID.eq.0)Then
        Do n=1,processCount-1
          processTo = n
          tag = 7000 + n
          Call MPI_SEND(sendArr,7,MPI_DOUBLE_PRECISION,processTo,tag,&
          MPI_COMM_WORLD,status,error)        
        End Do
      End If
      If(processID.gt.0)Then
        processFrom = 0
        tag = 7000 + processID
        Call MPI_RECV(recvArr,7,MPI_DOUBLE_PRECISION,processFrom,tag,&
        MPI_COMM_WORLD,status,error)
        sendArr = recvArr
      End If
      Call MPI_Barrier(MPI_COMM_WORLD,error)
    End If    
! Reassign values
    rss = sendArr(1)
    a = sendArr(2)
    lA = sendArr(3)
    b = sendArr(4)
    lB = sendArr(5)
    c = sendArr(6)
    lC = sendArr(7)
  End Subroutine M_bestParametersB
  
End Program doubleDecay
