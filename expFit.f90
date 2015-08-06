PROGRAM expFit
! University of Birmingham
! Ben Palmer
!
! Advice from Claude Leibovici
! Book by Jean Jacquelin:  https://www.scribd.com/doc/14674814/Regressions-et-equations-integrales
!
! 
! Example run:
! ./expFit.x data1.in 2      ! reads data1.in, fits double exp term function
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
      fileRow = trim(adjustl(fileRow))
      If(fileRow(1:1).ne."!")Then   ! Skip comment rows
        If(fileRow(1:4).ne."    ")Then
          dataRows = dataRows + 1
          Read(fileRow,*) bufferA, bufferB
          Read(bufferA,*) rawData(dataRows,1)
          Read(bufferB,*) rawData(dataRows,2)        
        End If  
      End If
    End Do
  End Subroutine readData
  
  Subroutine M_expFit(dataRows)
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i, dataRows
    Real(kind=DoubleReal), Dimension(1:dataRows,1:2) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:3) :: outputA
    Real(kind=DoubleReal), Dimension(1:5) :: outputB
    Real(kind=DoubleReal), Dimension(1:7) :: outputC
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
    If(fitType.eq.1)Then  ! Not working yet
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
      outputB = DoubleDecayFit(dataPoints)
! End Time    
      Call cpu_time(programEndTime)
! Output    
      If(processID.eq.0)Then
        print "(A21,F14.7)","Time:                ",(programEndTime-programStartTime)
        print "(A21,E14.7)","Final RSS Value:     ",outputB(5)
        print "(A33)","f(x) = a exp(lA x) + b exp(lB x) "
        print "(A21,E14.7)","a:                   ",outputB(1)
        print "(A21,E14.7)","lA:                  ",outputB(2)
        print "(A21,E14.7)","b:                   ",outputB(3)
        print "(A21,E14.7)","lB:                  ",outputB(4)      
      End If
    End If  
! Fit - three exponential terms 
    If(fitType.eq.3)Then  ! Not working yet
      outputC = TripleDecayFit(dataPoints)
! End Time    
      Call cpu_time(programEndTime)
! Output    
      If(processID.eq.0)Then
        print "(A21,F14.7)","Time:                ",(programEndTime-programStartTime)
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

 
  Function SingleDecayFit(dataPoints, convergenceThresholdIn) RESULT (output)
! Fits double exponential to data
! f(x) = a exp(lA)
    Implicit None  !Force declaration of all variables
! In
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal), Optional :: convergenceThresholdIn
    Real(kind=DoubleReal) :: convergenceThreshold
! Out    
    Real(kind=DoubleReal), Dimension(1:3) :: output
! Private    
    Integer(kind=StandardInteger) :: n
    Real(kind=DoubleReal) :: lnY
    Real(kind=DoubleReal) :: sumY, sumX_Y, sumX_X_Y
    Real(kind=DoubleReal) :: sumY_LnY, sumX_Y_LnY
! Linear regression
    Real(kind=DoubleReal), Dimension(1:2,1:2) :: xMatrix
    Real(kind=DoubleReal), Dimension(1:2) :: yMatrix, cMatrix
! LMA    
    Integer(kind=StandardInteger) :: i, k
    Real(kind=DoubleReal) :: convergence, rss, bestRSS, testRSS, lambda
    Real(kind=DoubleReal), Dimension(1:2) :: parameters, parameters_Last
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: R
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1),1:2) :: J
    Real(kind=DoubleReal), Dimension(1:2,1:size(dataPoints,1)) :: JT     ! Transpose Jacobian
    Real(kind=DoubleReal), Dimension(1:2,1:2) :: JTJ    ! (Jacobian Transpose * Jacobian)
    Real(kind=DoubleReal), Dimension(1:2,1:2) :: JTJ_Diag
    Real(kind=DoubleReal), Dimension(1:2) :: JTR                ! (Jacobian Transpose * Residuals)
    Real(kind=DoubleReal), Dimension(1:2) :: P      ! Change  
! Optional argument
    convergenceThreshold = 1.0D-8
    If(Present(convergenceThresholdIn))Then
      convergenceThreshold = convergenceThresholdIn
    End If    
! --------------------
! Linear regression
! --------------------
    sumY = 0.0D0
    sumX_Y = 0.0D0
    sumX_X_Y = 0.0D0
    sumY_LnY = 0.0D0
    sumX_Y_LnY = 0.0D0
    Do n=1,size(dataPoints,1)
      lnY = log(dataPoints(n,2))
      sumY = sumY + dataPoints(n,2)    
      sumX_Y = sumX_Y + dataPoints(n,1)*dataPoints(n,2)    
      sumX_X_Y = sumX_X_Y + dataPoints(n,1)*dataPoints(n,1)*dataPoints(n,2)    
      sumY_LnY = sumY_LnY + dataPoints(n,2)*lnY
      sumX_Y_LnY = sumX_Y_LnY + dataPoints(n,1)*dataPoints(n,2)*lnY
    End Do    
! x matrix
    xMatrix(1,1) = sumY
    xMatrix(1,2) = sumX_Y
    xMatrix(2,1) = sumX_Y
    xMatrix(2,2) = sumX_X_Y
! y matrix    
    yMatrix(1) = sumY_LnY
    yMatrix(2) = sumX_Y_LnY
! solve    
    xMatrix = InvertMatrix(xMatrix)
    cMatrix = matmul(xMatrix,yMatrix)
! --------------------
! LMA
! --------------------    
! set parameters    
    parameters(1) = exp(cMatrix(1))
    parameters(2) = cMatrix(2) 
    Do i=1,2   
      If(isnan(parameters(i)))Then  ! If fitting fails
        parameters(i) = 1.0D0
      End If  
    End Do
! LMA Opt
    convergence = 1.0D0
    lambda = 1.0D0
!----------------
! Start LMA loop
    Do n=1,1000  ! maximum 1000 loops
      rss = 0.0D0
! Make Jacobian and Residuals matrix
      Do i=1,size(dataPoints,1)
        R(i) = (parameters(1)*exp(parameters(2)*dataPoints(i,1)))-dataPoints(i,2)   ! f(x)-y
        J(i,1) = exp(parameters(2)*dataPoints(i,1))  ! d/dx1
        J(i,2) = dataPoints(i,1)*parameters(1)*exp(parameters(2)*dataPoints(i,1))  ! d/dx2
        rss = rss + R(i)**2
      End Do
      Do k=1,50 ! max 50
! calculate change matrix
        !***********     
        ! P = (JTJ+L*diag(JTJ))^(-1)(-1*JTR)   
        !***********      
! Transpose Jacobian
        JT = TransposeMatrix(J)
        JTJ = matmul(JT,J)
        JTJ_Diag = lambda*DiagMatrix(JTJ) ! Dampening Matrix
        JTJ = MatAdd(JTJ,JTJ_Diag) ! Recycle JTJ
        JTJ = InvertMatrix(JTJ) ! store inverse (recycle JTJ var)      
        JTR = matmul(JT,R)
        JTR = -1.0D0*JTR ! Recycle JTR var
        P = matmul(JTJ,JTR)  
! Store last loop values
        parameters_Last = parameters
! Update parameters      
        Do i=1,size(P)
          parameters(i) = parameters(i) + P(i)
        End Do       
! Calc RSS
        testRSS = 0.0D0
        Do i=1,size(dataPoints,1)
          testRSS = testRSS + &
          ((parameters(1)*exp(parameters(2)*dataPoints(i,1)))-dataPoints(i,2))**2
        End Do
! Delayed gratification scheme - 1.5*lambda or 0.2*lambda    
        If(testRSS.gt.rss)Then  ! If worse
          lambda = lambda * 1.5D0   
          parameters = parameters_Last      
          bestRSS = rss
        Else  ! If better
          lambda = lambda * 0.2D0  
          bestRSS = testRSS
          Exit
        End If
      End Do
      convergence = abs(testRSS-rss) 
  ! Breakout if convergence threshold met      
      If(convergence.lt.convergenceThreshold)Then
        Exit
      End If
! End LMA loop
!----------------
    End Do  
! Output   f(x) = a exp(b x)
    output(1) = parameters(1)  ! a
    output(2) = parameters(2)  ! b
    output(3) = bestRSS
  End Function SingleDecayFit    
    
  Function DoubleDecayFit(dataPoints, convergenceThresholdIn) RESULT (output)
! Fits double exponential to data
! f(x) = a exp(lA) + b exp(lB)
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Optional :: convergenceThresholdIn
    Real(kind=DoubleReal) :: convergenceThreshold
    Integer(kind=StandardInteger) :: i, n, k
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:5) :: output
! Lin Reg approx vars    
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: S
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: SS
    Real(kind=DoubleReal), Dimension(1:4) :: cMatrix, yMatrix
    Real(kind=DoubleReal), Dimension(1:4,1:4) :: xMatrix
    Real(kind=DoubleReal) :: sumX, sumY, sumX_X, sumX_Y, sumS_X, sumS_Y, sumSS_X, sumSS_Y
    Real(kind=DoubleReal) :: sumS, sumSS, sumS_S, sumS_SS, sumSS_SS
    Real(kind=DoubleReal) :: p1, q1
    Real(kind=DoubleReal), Dimension(1:2) :: cbMatrix, ybMatrix
    Real(kind=DoubleReal), Dimension(1:2,1:2) :: xbMatrix
    Real(kind=DoubleReal) :: sumB_B, sumB_N, sumN_N, sumB_Y, sumN_Y
    Real(kind=DoubleReal) :: beta, eta  
! LMA Vars
    Real(kind=DoubleReal) :: rss, testRSS, bestRSS, convergence, lambda  
    Real(kind=DoubleReal), Dimension(1:4) :: parameters, parameters_Last
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: R
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1),1:4) :: J
    Real(kind=DoubleReal), Dimension(1:4,1:size(dataPoints,1)) :: JT     ! Transpose Jacobian
    Real(kind=DoubleReal), Dimension(1:4,1:4) :: JTJ    ! (Jacobian Transpose * Jacobian)
    Real(kind=DoubleReal), Dimension(1:4,1:4) :: JTJ_Diag
    Real(kind=DoubleReal), Dimension(1:4) :: JTR                ! (Jacobian Transpose * Residuals)
    Real(kind=DoubleReal), Dimension(1:4) :: P      ! Change  
! Optional argument
    convergenceThreshold = 1.0D-8
    If(Present(convergenceThresholdIn))Then
      convergenceThreshold = convergenceThresholdIn
    End If
!   
!-----------------------------------------
! Approximate linear regression 
!-----------------------------------------
! Advice from Claude Leibovici
! Book by Jean Jacquelin:  https://www.scribd.com/doc/14674814/Regressions-et-equations-integrales
! 
! Init vars
   sumX = 0.0D0
   sumY = 0.0D0
   sumX_X = 0.0D0
   sumX_Y = 0.0D0
   sumS_X = 0.0D0
   sumS_Y = 0.0D0
   sumSS_X = 0.0D0
   sumSS_Y = 0.0D0
   sumS = 0.0D0
   sumSS = 0.0D0
   sumS_S = 0.0D0
   sumS_SS = 0.0D0
   sumSS_SS = 0.0D0
! Numeric integration to calc S and SS array
    n = size(dataPoints,1)      ! number of data points
    Do i=1,n
      If(i.eq.1)Then
        S(i) = 0.0D0
        SS(i) = 0.0D0
      Else  
        S(i) = S(i-1) + 0.5D0*(dataPoints(i,2)+dataPoints(i-1,2))*&
        (dataPoints(i,1)-dataPoints(i-1,1)) ! Numeric integration
        SS(i) = SS(i-1) + 0.5D0*(S(i)+S(i-1))*&
        (dataPoints(i,1)-dataPoints(i-1,1)) ! Numeric integration
      End If 
    End Do
! Sum    
    Do i=1,n
      sumX = sumX + dataPoints(i,1)
      sumY = sumY + dataPoints(i,2)
      sumX_X = sumX_X + (dataPoints(i,1)*dataPoints(i,1))
      sumX_Y = sumX_Y + (dataPoints(i,1)*dataPoints(i,2))
      sumS_X = sumS_X + S(i)*dataPoints(i,1)
      sumS_Y = sumS_Y + S(i)*dataPoints(i,2)
      sumSS_X = sumSS_X + SS(i)*dataPoints(i,1)
      sumSS_Y = sumSS_Y + SS(i)*dataPoints(i,2)
      sumS = sumS + S(i)
      sumSS = sumSS + SS(i)
      sumS_S = sumS_S + S(i)*S(i)
      sumS_SS = sumS_SS + S(i)*SS(i)
      sumSS_SS = sumSS_SS + SS(i)*SS(i)
    End Do
! Make y matrix    
    yMatrix(1) = sumSS_Y
    yMatrix(2) = sumS_Y
    yMatrix(3) = sumX_Y
    yMatrix(4) = sumY
! Make xMatrix
    xMatrix(1,1) = sumSS_SS
    xMatrix(1,2) = sumS_SS
    xMatrix(1,3) = sumSS_X
    xMatrix(1,4) = sumSS
    xMatrix(2,1) = sumS_SS
    xMatrix(2,2) = sumS_S
    xMatrix(2,3) = sumS_X
    xMatrix(2,4) = sumS
    xMatrix(3,1) = sumSS_X
    xMatrix(3,2) = sumS_X
    xMatrix(3,3) = sumX_X
    xMatrix(3,4) = sumX
    xMatrix(4,1) = sumSS
    xMatrix(4,2) = sumS
    xMatrix(4,3) = sumX
    xMatrix(4,4) = n 
! Solve set of equations
    xMatrix = InvertMatrix(xMatrix)
    cMatrix = matmul(xMatrix,yMatrix)
! calculate P and Q for next regression    
    p1 = 0.5D0*(cMatrix(2)+sqrt(cMatrix(2)*cMatrix(2)+4*cMatrix(1)))
    q1 = 0.5D0*(cMatrix(2)-sqrt(cMatrix(2)*cMatrix(2)+4*cMatrix(1)))
! Sum
    sumB_B = 0.0D0 
    sumB_N = 0.0D0
    sumN_N = 0.0D0
    sumB_Y = 0.0D0
    sumN_Y = 0.0D0
    Do i=1,n
      beta = exp(p1*dataPoints(i,1))
      eta = exp(q1*dataPoints(i,1))
      sumB_B = sumB_B + beta*beta
      sumB_N = sumB_N + beta*eta
      sumN_N = sumN_N + eta*eta
      sumB_Y = sumB_Y + beta*dataPoints(i,2)
      sumN_Y = sumN_Y + eta*dataPoints(i,2)
    End Do
! Make next x matrix
    xbMatrix(1,1) = sumB_B
    xbMatrix(1,2) = sumB_N
    xbMatrix(2,1) = sumB_N
    xbMatrix(2,2) = sumN_N
! Make next y matrix    
    ybMatrix(1) = sumB_Y
    ybMatrix(2) = sumN_Y
! Calc cb
    xbMatrix = InvertMatrix(xbMatrix)
    cbMatrix = matmul(xbMatrix,ybMatrix) 
!--------------------------------------------------
! LMA
!--------------------------------------------------    
    parameters(1) = cbMatrix(1)
    parameters(2) = p1
    parameters(3) = cbMatrix(2)
    parameters(4) = q1
    Do i=1,4    
      If(isnan(parameters(i)))Then  ! If fitting fails
        parameters(i) = 1.0D0
      End If  
    End Do
! LMA Opt
    convergence = 1.0D0
    lambda = 1.0D0
!----------------
! Start LMA loop
    Do n=1,1000  ! maximum 1000 loops
      rss = 0.0D0
! Make Jacobian and Residuals matrix
      Do i=1,size(dataPoints,1)
        R(i) = (parameters(1)*exp(parameters(2)*dataPoints(i,1))+&
        parameters(3)*exp(parameters(4)*dataPoints(i,1)))-dataPoints(i,2)   ! f(x)-y
        J(i,1) = exp(parameters(2)*dataPoints(i,1))  ! d/dx1
        J(i,2) = dataPoints(i,1)*parameters(1)*exp(parameters(2)*dataPoints(i,1))  ! d/dx2
        J(i,3) = exp(parameters(4)*dataPoints(i,1))  ! d/dx3
        J(i,4) = dataPoints(i,1)*parameters(3)*exp(parameters(4)*dataPoints(i,1))  ! d/dx4
        rss = rss + R(i)**2
      End Do
      Do k=1,50 ! max 50
! calculate change matrix
        !***********     
        ! P = (JTJ+L*diag(JTJ))^(-1)(-1*JTR)   
        !***********      
! Transpose Jacobian
        JT = TransposeMatrix(J)
        JTJ = matmul(JT,J)
        JTJ_Diag = lambda*DiagMatrix(JTJ) ! Dampening Matrix
        JTJ = MatAdd(JTJ,JTJ_Diag) ! Recycle JTJ
        JTJ = InvertMatrix(JTJ) ! store inverse (recycle JTJ var)      
        JTR = matmul(JT,R)
        JTR = -1.0D0*JTR ! Recycle JTR var
        P = matmul(JTJ,JTR)  
! Store last loop values
        parameters_Last = parameters
! Update parameters      
        Do i=1,size(P)
          parameters(i) = parameters(i) + P(i)
        End Do       
! Calc RSS
        testRSS = 0.0D0
        Do i=1,size(dataPoints,1)
          testRSS = testRSS + &
          ((parameters(1)*exp(parameters(2)*dataPoints(i,1))+&
          parameters(3)*exp(parameters(4)*dataPoints(i,1)))-dataPoints(i,2))**2
        End Do
! Delayed gratification scheme - 1.5*lambda or 0.2*lambda   
        If(testRSS.gt.rss)Then  ! If worse
          lambda = lambda * 1.5D0   
          parameters = parameters_Last  
          bestRSS = rss        
        Else  ! If better
          lambda = lambda * 0.2D0
          bestRSS = testRSS  
          Exit
        End If
      End Do
      convergence = abs(testRSS-rss) 
  ! Breakout if convergence threshold met      
      If(convergence.lt.convergenceThreshold)Then
        Exit
      End If
! End LMA loop
!----------------
    End Do  
! Output   f(x) = a exp(b x) + c exp(d x)
    output(1) = parameters(1)  ! a
    output(2) = parameters(2)  ! b
    output(3) = parameters(3)  ! c
    output(4) = parameters(4)  ! d
    output(5) = bestRSS
  End Function DoubleDecayFit 
   
  Function TripleDecayFit(dataPoints, convergenceThresholdIn) RESULT (output)
! Fits double exponential to data
! f(x) = a exp(lA) + b exp(lB)
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Optional :: convergenceThresholdIn
    Real(kind=DoubleReal) :: convergenceThreshold
    Integer(kind=StandardInteger) :: i, n, k
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:7) :: output
! Lin Reg approx vars    
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: S
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: SS
    Real(kind=DoubleReal), Dimension(1:4) :: cMatrix, yMatrix
    Real(kind=DoubleReal), Dimension(1:4,1:4) :: xMatrix
    Real(kind=DoubleReal) :: sumX, sumY, sumX_X, sumX_Y, sumS_X, sumS_Y, sumSS_X, sumSS_Y
    Real(kind=DoubleReal) :: sumS, sumSS, sumS_S, sumS_SS, sumSS_SS
    Real(kind=DoubleReal) :: p1, q1
    Real(kind=DoubleReal), Dimension(1:2) :: cbMatrix, ybMatrix
    Real(kind=DoubleReal), Dimension(1:2,1:2) :: xbMatrix
    Real(kind=DoubleReal) :: sumB_B, sumB_N, sumN_N, sumB_Y, sumN_Y
    Real(kind=DoubleReal) :: beta, eta  
! LMA Vars
    Real(kind=DoubleReal) :: rss, testRSS, bestRSS, convergence, lambda  
    Real(kind=DoubleReal), Dimension(1:6) :: parameters, parameters_Last
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: R
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1),1:6) :: J
    Real(kind=DoubleReal), Dimension(1:6,1:size(dataPoints,1)) :: JT     ! Transpose Jacobian
    Real(kind=DoubleReal), Dimension(1:6,1:6) :: JTJ    ! (Jacobian Transpose * Jacobian)
    Real(kind=DoubleReal), Dimension(1:6,1:6) :: JTJ_Diag
    Real(kind=DoubleReal), Dimension(1:6) :: JTR                ! (Jacobian Transpose * Residuals)
    Real(kind=DoubleReal), Dimension(1:6) :: P      ! Change  
! Optional argument
    convergenceThreshold = 1.0D-8
    If(Present(convergenceThresholdIn))Then
      convergenceThreshold = convergenceThresholdIn
    End If
!   
!-----------------------------------------
! Approximate linear regression 
!-----------------------------------------
! Advice from Claude Leibovici
! Book by Jean Jacquelin:  https://www.scribd.com/doc/14674814/Regressions-et-equations-integrales
! 
! Init vars
   sumX = 0.0D0
   sumY = 0.0D0
   sumX_X = 0.0D0
   sumX_Y = 0.0D0
   sumS_X = 0.0D0
   sumS_Y = 0.0D0
   sumSS_X = 0.0D0
   sumSS_Y = 0.0D0
   sumS = 0.0D0
   sumSS = 0.0D0
   sumS_S = 0.0D0
   sumS_SS = 0.0D0
   sumSS_SS = 0.0D0
! Numeric integration to calc S and SS array
    n = size(dataPoints,1)      ! number of data points
    Do i=1,n
      If(i.eq.1)Then
        S(i) = 0.0D0
        SS(i) = 0.0D0
      Else  
        S(i) = S(i-1) + 0.5D0*(dataPoints(i,2)+dataPoints(i-1,2))*&
        (dataPoints(i,1)-dataPoints(i-1,1)) ! Numeric integration
        SS(i) = SS(i-1) + 0.5D0*(S(i)+S(i-1))*&
        (dataPoints(i,1)-dataPoints(i-1,1)) ! Numeric integration
      End If 
    End Do
! Sum    
    Do i=1,n
      sumX = sumX + dataPoints(i,1)
      sumY = sumY + dataPoints(i,2)
      sumX_X = sumX_X + (dataPoints(i,1)*dataPoints(i,1))
      sumX_Y = sumX_Y + (dataPoints(i,1)*dataPoints(i,2))
      sumS_X = sumS_X + S(i)*dataPoints(i,1)
      sumS_Y = sumS_Y + S(i)*dataPoints(i,2)
      sumSS_X = sumSS_X + SS(i)*dataPoints(i,1)
      sumSS_Y = sumSS_Y + SS(i)*dataPoints(i,2)
      sumS = sumS + S(i)
      sumSS = sumSS + SS(i)
      sumS_S = sumS_S + S(i)*S(i)
      sumS_SS = sumS_SS + S(i)*SS(i)
      sumSS_SS = sumSS_SS + SS(i)*SS(i)
    End Do
! Make y matrix    
    yMatrix(1) = sumSS_Y
    yMatrix(2) = sumS_Y
    yMatrix(3) = sumX_Y
    yMatrix(4) = sumY
! Make xMatrix
    xMatrix(1,1) = sumSS_SS
    xMatrix(1,2) = sumS_SS
    xMatrix(1,3) = sumSS_X
    xMatrix(1,4) = sumSS
    xMatrix(2,1) = sumS_SS
    xMatrix(2,2) = sumS_S
    xMatrix(2,3) = sumS_X
    xMatrix(2,4) = sumS
    xMatrix(3,1) = sumSS_X
    xMatrix(3,2) = sumS_X
    xMatrix(3,3) = sumX_X
    xMatrix(3,4) = sumX
    xMatrix(4,1) = sumSS
    xMatrix(4,2) = sumS
    xMatrix(4,3) = sumX
    xMatrix(4,4) = n 
! Solve set of equations
    xMatrix = InvertMatrix(xMatrix)
    cMatrix = matmul(xMatrix,yMatrix)
! calculate P and Q for next regression    
    p1 = 0.5D0*(cMatrix(2)+sqrt(cMatrix(2)*cMatrix(2)+4*cMatrix(1)))
    q1 = 0.5D0*(cMatrix(2)-sqrt(cMatrix(2)*cMatrix(2)+4*cMatrix(1)))
! Sum
    sumB_B = 0.0D0 
    sumB_N = 0.0D0
    sumN_N = 0.0D0
    sumB_Y = 0.0D0
    sumN_Y = 0.0D0
    Do i=1,n
      beta = exp(p1*dataPoints(i,1))
      eta = exp(q1*dataPoints(i,1))
      sumB_B = sumB_B + beta*beta
      sumB_N = sumB_N + beta*eta
      sumN_N = sumN_N + eta*eta
      sumB_Y = sumB_Y + beta*dataPoints(i,2)
      sumN_Y = sumN_Y + eta*dataPoints(i,2)
    End Do
! Make next x matrix
    xbMatrix(1,1) = sumB_B
    xbMatrix(1,2) = sumB_N
    xbMatrix(2,1) = sumB_N
    xbMatrix(2,2) = sumN_N
! Make next y matrix    
    ybMatrix(1) = sumB_Y
    ybMatrix(2) = sumN_Y
! Calc cb
    xbMatrix = InvertMatrix(xbMatrix)
    cbMatrix = matmul(xbMatrix,ybMatrix) 
!--------------------------------------------------
! LMA
!--------------------------------------------------    
    parameters(1) = cbMatrix(1)
    parameters(2) = p1
    parameters(3) = cbMatrix(2)
    parameters(4) = q1
    parameters(5) = 1.0D0 
    parameters(6) = 1.0D0
    Do i=1,4    
      If(isnan(parameters(i)))Then  ! If fitting fails
        parameters(i) = 1.0D0
      End If  
    End Do
! LMA Opt
    convergence = 1.0D0
    lambda = 1.0D0
!----------------
! Start LMA loop
    Do n=1,1000  ! maximum 1000 loops
      rss = 0.0D0
! Make Jacobian and Residuals matrix
      Do i=1,size(dataPoints,1)
        R(i) = (parameters(1)*exp(parameters(2)*dataPoints(i,1))+&
        parameters(3)*exp(parameters(4)*dataPoints(i,1))+&
        parameters(5)*exp(parameters(6)*dataPoints(i,1))&
        )-dataPoints(i,2)   ! f(x)-y
        J(i,1) = exp(parameters(2)*dataPoints(i,1))  ! d/dx1
        J(i,2) = dataPoints(i,1)*parameters(1)*exp(parameters(2)*dataPoints(i,1))  ! d/dx2
        J(i,3) = exp(parameters(4)*dataPoints(i,1))  ! d/dx3
        J(i,4) = dataPoints(i,1)*parameters(3)*exp(parameters(4)*dataPoints(i,1))  ! d/dx4
        J(i,5) = exp(parameters(6)*dataPoints(i,1))  ! d/dx5
        J(i,6) = dataPoints(i,1)*parameters(5)*exp(parameters(6)*dataPoints(i,1))  ! d/dx6
        rss = rss + R(i)**2
      End Do
      Do k=1,50 ! max 50
! calculate change matrix
        !***********     
        ! P = (JTJ+L*diag(JTJ))^(-1)(-1*JTR)   
        !***********      
! Transpose Jacobian
        JT = TransposeMatrix(J)
        JTJ = matmul(JT,J)
        JTJ_Diag = lambda*DiagMatrix(JTJ) ! Dampening Matrix
        JTJ = MatAdd(JTJ,JTJ_Diag) ! Recycle JTJ
        JTJ = InvertMatrix(JTJ) ! store inverse (recycle JTJ var)      
        JTR = matmul(JT,R)
        JTR = -1.0D0*JTR ! Recycle JTR var
        P = matmul(JTJ,JTR)  
! Store last loop values
        parameters_Last = parameters
! Update parameters      
        Do i=1,size(P)
          parameters(i) = parameters(i) + P(i)
        End Do       
! Calc RSS
        testRSS = 0.0D0
        Do i=1,size(dataPoints,1)
          testRSS = testRSS + &
          ((parameters(1)*exp(parameters(2)*dataPoints(i,1))+&
          parameters(3)*exp(parameters(4)*dataPoints(i,1))+&
          parameters(5)*exp(parameters(6)*dataPoints(i,1))&
          )-dataPoints(i,2))**2
        End Do
! Delayed gratification scheme - 1.5*lambda or 0.2*lambda  
        If(testRSS.gt.rss)Then  ! If worse
          lambda = lambda * 1.5D0   
          parameters = parameters_Last  
          bestRSS = rss        
        Else  ! If better
          lambda = lambda * 0.2D0
          bestRSS = testRSS  
          Exit
        End If
      End Do
      convergence = abs(testRSS-rss) 
  ! Breakout if convergence threshold met      
      If(convergence.lt.convergenceThreshold)Then
        Exit
      End If
! End LMA loop
!----------------
    End Do  
! Output   f(x) = a exp(b x) + c exp(d x)
    output(1) = parameters(1)  ! a
    output(2) = parameters(2)  ! b
    output(3) = parameters(3)  ! c
    output(4) = parameters(4)  ! d
    output(5) = parameters(5)  ! c
    output(6) = parameters(6)  ! d
    output(7) = bestRSS
  End Function TripleDecayFit   
    
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
  
  Function DoubleDecayFitRSS(dataPoints, a, b, lA, lB) RESULT (rss)
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
  End Function DoubleDecayFitRSS 

  
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
  
    Function DiagMatrix(iMatrix) RESULT (oMatrix)
! Takes diagonal of a matrix
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: row,col
    Real(kind=DoubleReal), Dimension(:,:) :: iMatrix
    Real(kind=DoubleReal), Dimension(1:size(iMatrix,1),1:size(iMatrix,2)) :: oMatrix
! Transpose
    If(size(iMatrix,1).ne.size(iMatrix,2))Then
      oMatrix = iMatrix
    Else
      Do row=1,size(iMatrix,1)
        Do col=1,size(iMatrix,2)
          If(row.eq.col)Then
            oMatrix(row,col) = 1.0D0*iMatrix(row,col)
          Else
            oMatrix(row,col) = 0.0D0
          End If
        End Do
      End Do  
    End If      
  End Function DiagMatrix
  
  Function MatAdd(xMatrix,yMatrix) RESULT (oMatrix)
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: i,j
    Real(kind=DoubleReal), Dimension(:,:) :: xMatrix, yMatrix
    Real(kind=DoubleReal), Dimension(1:size(xMatrix,1),1:size(xMatrix,2)) :: oMatrix
! Initialise variables
    oMatrix = 0.0D0
! Add matrices
    If(size(xMatrix,1).eq.size(yMatrix,1).and.size(xMatrix,2).eq.size(yMatrix,2))Then
      Do i=1,size(xMatrix,1)
        Do j=1,size(xMatrix,2)
          oMatrix(i,j) = xMatrix(i,j) + yMatrix(i,j)
        End Do
      End Do
    End If
  End Function MatAdd
  
End Program expFit
