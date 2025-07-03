
    
program jump 
    
	use phys_data                                                               !import libraries
	use utilities
	use marche

	implicit none                                                               !set program to not imply variable type 

    ! ####################################################################
    ! #################### 1. Define variables ###########################    
    ! ####################################################################
    integer             :: switch,duplicate,unique_count,max_species, max_duplicates,z,p
    integer, allocatable:: unique_species(:),event_species(:)
  	integer			    :: ii, iii,iiij,jj, jjj
  	integer			    :: kk, check
    integer             :: i_ind, j_ind, k_ind, print_accretion, file_index     ! Indexing for file generation
    integer             :: print_grid_data,print_evap_data,print_empty_data,print_thickness_data, print_nn_grid_data
    integer             :: print_event_count, print_event_interval, event_count
    double precision    :: TPD_step, time_per_TPD
  	integer			    :: nempty, na
  	integer			    :: i, j, k 
    integer			    :: nmol, surftot  
  	integer			    :: ngot    			
    INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(18)
    INTEGER(KIND=i8) :: print_interval, tmax
  	integer			    :: ev, evna, numev1       
    integer             :: mlycount 
  	integer			    :: npt(nspecies), neigh(4)
  	integer			    :: tri, stored_na,q
  	integer			    :: print
    integer             :: iatnew, jatnew, katnew, try, delta_i, delta_j, delta_k, par
    integer             :: max_TPD_temp                                         ! Max temperature to perform TPD
    integer             :: start_count, end_count, count_rate, elapsed_ticks    ! Variables for runtime counter
    integer             :: move_i, move_j, move_k                               ! Variables for moving species
    integer             :: sumval, iatl,jatl,k_froml, di,dj,dk
    integer             :: shuffle(4), locations(4), i_loc(4),j_loc(4), k_loc(4),moved_dimer(4),temp
    
  	double precision	:: tps, thick, Tgnext, Tg_max, temp_evan
  	double precision	:: dens,surface,volume
    double PRECISION    :: mly
    
    character(len=50)   :: fname                
    character(len=256)  :: path
    
    logical             :: path_exists
    
    real                :: elapsed_time, ang                                        ! Elapsed runtime
    
    external            :: system
    
    ! ####################################################################
    ! ##################### 2. Input variables ###########################    
    ! ####################################################################
    
    ! #################### 2.1 Print conditions ##########################
    
    file_index = 1
    
    print_event_count = 0                                                       ! Set to 1 to print at every xxx events                                             
    print_event_interval = 10000                                                 ! Print every xxxx events
    
    print_accretion = 1                                                         ! Set to 1 to turn on printing the accretion .dat files
    print_grid_data = 1                                                         ! Set to 1 to turn on printing the grid .dat files
    print_evap_data = 1                                                         ! Set to 1 to turn on printing the evap data
    print_empty_data = 1                                                        ! Set to 1 to turn on printing the empty grid data
    print_thickness_data = 1                                                    ! Set to 1 to turn on printing the thickness data
    print_interval = INT(60,i8) * INT(60,i8) * INT(24,i8) * INT(365,i8) * INT(1000,i8) 

    ! ################### 2.2 Gas/Dust variables #########################  
    
    Th=293.15d0                                                                    ! Gas temperature
    Tg=10d0                                                                     ! Dust temperature
    Tgnext= Tg                                                                 ! temperature during TPD (is increased by ramp)
    Tg_max = 200d0
    nH2O=1d4                                                                 ! density cm-3 of water
    nCO2=0.13*nH2O!
    nAr=0.01*nH2O                                                              ! Density Ar
    nKr=0.01*nH2O                                                                    ! Density Kr
    nXe=0.01*nH2O                                                                    ! Density Xe
     
    ! ##################### 2.3 Grid Variables ########################### 
    
    nptg=40                                                                     ! Dimension of the grid
    nsite=(0.5*nptg)*(0.5*nptg)                                                         ! Number of sites in x,y plan
    mly=50 !15                                                                     ! number of monolayers — Brown: 0.55, Krucz: 2, Fraser: 3.7
    top=int(20*mly)                                                                     ! dimension of the grid in z [Have 100 per 3 ml]
    na_max=floor(nsite*mly)                                                     ! number of atoms/molecules that will be send on the grid from the gas phase
    
    ! ##################### 2.4 TPD variables ############################
    
    max_TPD_temp = 200                                                          ! Max temperature for TPD
    time_per_TPD = INT(60,i8) * INT(60,i8) * INT(24,i8) * INT(365,i8) * INT(1,i8)  !1 year                                                           ! Set time per TPD
    TPD_step = 1                                                                ! Temperature step made each time
    tmax = INT(60,i8) * INT(60,i8) * INT(24,i8) * INT(365,i8) * INT(10000,i8) 
    ! ####################################################################
    ! #################### 3. Allocate arrays ############################    
    ! ####################################################################
    
    allocate (event_species(10000),unique_species(10000))
    allocate (iat(0:na_max),jat(0:na_max),kat(0:na_max))                        ! i,j,k coordinates of each species na
	allocate (dina(0:na_max),djna(0:na_max),dkna(0:na_max))                     !
	allocate (iat0(0:na_max),jat0(0:na_max))                                    ! initial i,j coordinates of each species na
	allocate (k_from(0:na_max))                                  !                                                  !
	allocate (nn(0:na_max))                                                     !
    allocate (nnlist(0:na_max))
	allocate (nb(0:na_max,0:27))                                                !
	allocate (evt(na_max))                                                      ! ev=1 move, ev=5 evap, ev=7 accretion 
	allocate (tlist1(na_max), num(na_max))                                      !
    allocate (Ei(0:na_max))                                                     !
    allocate (grid(0:nptg,0:nptg,-1:top-1))                                     ! Define grid size
    allocate (empty(0:nptg,0:nptg,-1:top-1))                                    ! Equal to grid size, empty grid for data
    allocate (attempt_rate(na_max,5))                                             ! Track species evap attempt rate
    allocate (hist_x(na_max,10),hist_y(na_max,10),hist_z(na_max,10))             ! Track species history
    allocate (evan_list(na_max,4))
    allocate (dimer(na_max,4),dimer_weak(na_max))
    allocate (theta(na_max), phi(na_max))
    allocate (tlast(na_max),t_evan(na_max,4))
    allocate (Eia(0:na_max))
    allocate (last_move(na_max,3))
    ! ####################################################################
    ! ####################### 4. Open Files ##############################    
    ! ####################################################################

     
    path = 'grid_data'                                                          ! Set file path
 
    INQUIRE(FILE=path, EXIST=path_exists)                                       ! Check if the folder exists
            
    if (.not. path_exists) then                                                 ! If the folder does not exist, create it
        call system('mkdir ' // trim(path))
    endif
           
	open (45, file="grid_data/taccret.out", status="unknown")                   ! File that lists species, their type and time to accrete
	open (46, file="grid_data/taccreti.out", status="unknown")                  ! File that orders species and assigns timestamp of accretion
    open (47, file="grid_data/thickness3.out", status="unknown",access="append")! File to define Ice thickness
    open (53, file="grid_data/ice1evap.out", status="unknown",access="append")  ! File to define evaporating molecules
        
    close (47)                                                                  ! Close files
    close (53)

    ! ####################################################################
    ! ###################### 5. Start Program ############################    
    ! ####################################################################
        
	print *, "------------------------------------------------- "
	print *, "JUMP> jump!"
        
	print*,'JUMP> Number of atoms sent on the surface', na_max

127 format(A6)   

    ! ####################################################################
    ! ############## 6. Initialisation and Definition ####################    
    ! ####################################################################
    
    ! ########################## 6.1 initalise ###########################

	idum=1                                                                     ! this is the seed for random numbers
    mlycount=0  
    nnlist = 0! Count number of monolayers
	grid=0                                                                      ! initialise grid
	empty=0                                                                     ! intialise empty grid
	evap=0.0_dp                                                                 ! intialise evaporation grid
	npt=0                                                                       ! 
    nspec=0                                                                     ! Set intial number of species to 0
    nmolstay=0                                                                  !
    nmolevap=0                                                                  !
    iat=-1                                                                      ! when the atoms are in the gas (before arriving on the grid) they are declared out of grid with iat=jat=-1
    jat=-1                                                                      ! out of grid
    numev1=0                                                                    !
    iiij=0                                                                      !
    event_count = 0                                                             ! Initialise counter for events occuring
    attempt_rate = 0
    hist_x = 0
    hist_y = 0
    hist_z = 0
    check = 0
    dimer(:,:) = 0.0d0
    dimer_weak(:) = 0.0d0
    tlast(:) = 0
    t_evan(:,:) = 0
    switch = 0
    duplicate = 0
    unique_count = 0
    max_species = 0
    max_duplicates = 0
    z = 0
    p = 0

    ii = 0
    iii = 0
    jj = 0
    jjj = 0
    kk = 0
    check = 0

    i_ind = 0
    j_ind = 0
    k_ind = 0

    print_nn_grid_data = 0

    tri = 0
    stored_na = 0
    q = 0
    print = 0

    iatnew = 0
    jatnew = 0
    katnew = 0
    try = 0

    delta_i = 0
    delta_j = 0
    delta_k = 0
    par = 0

    start_count = 0
    end_count = 0
    count_rate = 0
    elapsed_ticks = 0

    move_i = 0
    move_j = 0
    move_k = 0

    sumval = 0
    iatl = 0
    jatl = 0
    k_froml = 0
    di = 0
    dj = 0
    dk = 0
    
    npt(:) = 0
    neigh(:) = 0
    shuffle(:) = 0
    locations(:) = 0
    i_loc(:) = 0
    j_loc(:) = 0
    k_loc(:) = 0
    moved_dimer(:) = 0
    last_move(:,:) = 0

    call SYSTEM_CLOCK(COUNT_RATE=count_rate)    
    call SYSTEM_CLOCK(COUNT=start_count)                                        ! Initialise timer
    
    elapsed_time = 0

    
    ! ################### 6.3 Accretion rate parameters ##################
    
	sigma =(nptg)*(nptg)*(100*a_pp)**2		                                    ! cross section of grain in cm2
	sigma_at=(100*a_pp)**2                                                      !
	stick = 1d0               			                                        ! sticking coefficient
	v(1) = sqrt(8*k_bz*Th/(pi*amu))*100                                         ! velocity of H atoms cm/s
	v(10) = v(1)/sqrt(18.0)                                                     ! velocity of H2O  cm/s
    v(36) = v(1)/sqrt(44.0)             	                                    ! velocity of CH4  cm/s
    v(37) = v(1)/sqrt(16.0)             	                                    ! velocity of CH4  cm/s
    v(38) = v(1)/sqrt(40.0)             	                                    ! velocity of Ar  cm/s
    v(39) = v(1)/sqrt(84.0)              	                                    ! velocity of Kr  cm/s
    v(40) = v(1)/sqrt(131.0)               	                                    ! velocity of Xe  cm/s
 
	k_accP=1d-99                                                                ! Set all others close to zero
	k_accP(10) = nH2O * sigma * stick * v(10)                                   ! H2O
    k_accP(36) = nCO2 * sigma * stick * v(36)                                   ! CH4
    k_accP(38) = nAr * sigma * stick * v(38)                                    ! Ar
    k_accP(39) = nKr * sigma * stick * v(39)                                    ! Kr
    k_accP(40) = nXe * sigma * stick * v(40)                                    ! Xe
        
    ! ####################################################################
    ! ######################### 7. Simulation ############################    
    ! ####################################################################

    ! ############### 7.1 Generate accretion data ########################
 
    call accretion 
    
	do kk=1,na_max                                                              ! Store all accretion data to grids and files
		write(46,*)kk, taccreti(kk)                                             ! Open file taccreti and write species # and time of accretion
        tlist1(kk)=taccreti(kk)                                                 ! Set list1 values equal to taccreti
		num(kk)=kk                                                              ! Store species numbers
		evt(kk)=8-mod(indice(kk),2)                                             ! Set and store event numbers (7 = accretion)
		kat(kk)=(indice(kk)+mod(indice(kk),2))/2                                ! Set and store species ID (e.g. 10 = H2O)
        k_from(kk)=-2                                                           ! Set initial vertical position out of the grid
    enddo
    
    print*, ' --- '
    print*, 'accretion data generated'                                          ! Confirm data generation
    print*, ' Enter phase 0 - Accretion '
    print*, ' --- '
  
    ! ######### 7.2 Initialise for accretion of first atom ###############

	na=1                                                                        ! Set species ID to 1
	ngot=0                                                                      !
	ii=200                                                                      !
	jj=na_max                                                                   !
    tri=2                                                                       !
	evna=evt(1)                                                                 !
    bind=0                                                                      !
    tlist0=0                                                                    !
    iphase=0                                                                    ! Set phase to accretion (iphase=0) (iphase=1 is TPD and iphase=2 is waiting time)
    jjj=0

    ! ########## 7.3 Start the loop for accretion and TPD ################

    do while (Tg<max_TPD_temp .AND. tlist0<tmax)                             ! Run untill gas temperature Tg hits max TPD temperature
12      continue    
    ! ################ 7.3.1 Order molecules in time #####################
        call piksr3(na_max,tlist1,num,evt)                       ! Subroutine to order molecules in time    
        na=num(1)                                                               ! Species that is acting first
        ev=evt(1)                                                               ! first event
        evna=ev
    ! ############### 7.3.2 ramp temperature in TPD ######################
        
        if (na==na_max .and. iphase==0) then   
            print=1
        endif
        
        event_count = event_count + 1
        
        if (mod(event_count,print_event_interval)==0) then
            print*, event_count
            print*, 'simulated time: ', tlist1(1)/(60*60*24*365.25), 'years | n_sp: ', nmolstay(38:40)
        endif

        if (iphase==0 .and. event_count==print_event_interval .and. print_event_count==1) then  ! Print every xxx events
            print=1
        endif

        if (iphase==1 .and.   tlist1(1)-tlist0>time_per_TPD) then               ! TPD phase - Set next temperature and print
            print = 1
            attempt_rate(:,:) = 0
            tlist0=tlist0+time_per_TPD
            if (Tg .ge. Tg_max) then
                iphase = 2
            else
                 Tg=Tg+TPD_step     
            endif
            do iii=1,na_max                                                     ! for all na molecules on the surface , the new time is calculated
                na=num(iii)              
                stored_na = na
        
                if (k_from(na)>-1) then                                         ! if their ifrom is -1, then the species are already evaporated
                    call add_one(na, tps, ev)         ! Call subroutine to calculate new event time for na                                           
                    evna=ev
                    evt(iii)=ev                                                 ! Store event type
                    tlist1(iii)=tlist0+tps                                      ! Store new event time for reordeing
                    num(iii)=na                                                 ! Store species number
                    if (stored_na==na) then
                    else
                        print*, stored_na,na, 'mismatch'
                    endif
                endif
            enddo
            neigh(:) = 0
            do i=1,na_max
                if (nnlist(i)>0) then
                    neigh(nnlist(i)) = neigh(nnlist(i)) + 1
                endif
            enddo   
            print*, 'neighbour distribution = ', neigh

            do i = 1, SUM(nmolstay) - 1
                do j = i + 1, SUM(nmolstay)
                    if (evan_list(i,1) > evan_list(j,1)) then
                        temp_evan = evan_list(i,1)
                        evan_list(i,1) = evan_list(j,1)
                        evan_list(j,1) = temp_evan
                    end if
                end do
            end do
            call piksr3(na_max,tlist1,evt,num)                                  ! reset the order of the event list
            
            na=num(1)
            ev=evt(1)
            evna=ev 
             
        elseif (iphase==2 .and. tlist1(1)-tlist0>print_interval) then           ! waiting phase - print every xxs
            print=1
            tlist0=tlist0+print_interval
        endif

    ! #################### 7.3.3 Print data ##############################
        if (print==1) then    
            print = 0
            call SYSTEM_CLOCK(COUNT=end_count)                                  ! Stop count

            elapsed_ticks = end_count - start_count                             ! calculate elapsed time
            elapsed_time = elapsed_time + real(elapsed_ticks) / real(count_rate)
            
            if (iphase==0) then
                print*, 'nspecies: ', SUM(nmolstay)
                print*, ' Temperature: ', Tg, 'Simulated time: ',tlist1(1)/(60*60*24*365.25), 'years'
            elseif (iphase .ge. 1) then
                print*, 'Temperature: ', Tg, 'Half-life: ', (evan_list(max(1,INT(SUM(nmolstay)/2)),1)/(60*60*24*365.25))
                print*, 'Number of accreted species: ', SUM(nmolstay), 'of which H2O & CO2: ', nmolstay(10), nmolstay(36), 'of which Vol (Ar,Kr,Xe) : ', nmolstay(38:40), 'Simulated time: ',tlist1(1)/(60*60*24*365.25), 'years'
            endif
            print *, 'Events occured: ', event_count
            print *, "Elapsed time (seconds):", elapsed_time
            print*, " --- "
            event_count = 0                                                     ! Reset count
            call SYSTEM_CLOCK(COUNT=start_count)                                ! restart counter
            
            nn(0)=0
            
            call counting(nempty,surftot,thick,nmol)                            ! 

           
    ! ################## 7.3.3.4 Print thickness data ####################
            if (print_thickness_data==1) then
                open(47,file="grid_data/thickness3.out", status="old",access="append")       
                write (47,"(1E14.4,5I6,1E14.6,10I8,3E14.6,I8)")Tg,&
                    thick/nsite,nempty,surface,volume,dens,nmol
                close(47)
            endif  
     ! ################## 7.3.3.4 Print thickness data ####################
            if (print_grid_data==1) then
                WRITE(fname, '(a,a,i4.4,a,a)') trim(path) // '/', 'grid0', file_index, '.dat'
                OPEN(56, file=fname, form='formatted')
                do k_ind = -1, top -1  ! Loop over the 3rd dimension
                    do j_ind = 0, nptg  ! Loop over the 2nd dimension
                        do i_ind = 0, nptg  ! Loop over the 1st dimension
                            if (grid(i_ind, j_ind, k_ind) .ne. 0) then
                                write(56, '(3I8, I8, F40.6, F8.2)') i_ind, j_ind, k_ind, kat(grid(i_ind, j_ind, k_ind)), tlist0, Tg  ! Write each element
                            endif
                        end do
                    end do
                end do
                close(56)
            endif
       
            file_index=file_index+1                                             ! this is the file counter

            dens=nmol*18*1.66d-24/(thick*(1.58d-8)**3)                          ! density of the ice (see paper MC for description)
            surface=surftot*8.647                                               ! See Cazaux 2015
            volume=nempty*10.74                                                 ! See Cazaux 2015
            
            if (tlist1(1) >1d98) then
                print*, 'all species desorbed'
                stop
            endif
        endif
        
    ! ####################################################################
    ! ################# 8. EVAPORATION (Ev=5) ############################
    ! ####################################################################
    ! If atom evaporates, remove from list.. and look for new event

        if (evna == 5) then                                                     ! the molecule will evaporate
            
            tlast(na)=tlist1(1)
37          continue
            if (tlist1(1) >1d98) then                                           ! all species desorbed, end simulation
                goto 12
            endif
            move_i = 0
            move_j = 0
            move_k = 0
            
            if (theta(na) == 0.0d0 .AND. phi(na)==0.0d0) then
                ang = 2*pi
            
                if (k_from(na) == 0) then
                    ang = pi
                endif
            
                call ran2(idum, xrand)                                              ! random angle of i descent between -80,80 degrees
		        theta(na) = (2*pi - 2*pi*xrand)
                call ran2(idum, xrand)                                              ! random angle of j descent between -80,80 degrees
                phi(na) = (ang/2 - ang*xrand)
            endif
            
            check = 0
            ii = 0

            do while (check==0)                                                       ! the grid is scanned to the top to see if na is on the surface layer
88              continue
                ii = ii + 1
                delta_i = nint((ii)*SIN(phi(na))*COS(theta(na)))
                delta_j = nint((ii)*SIN(phi(na))*SIN(theta(na)))
                delta_k = nint((ii)*COS(phi(na)))
                
                if (delta_k == 0 .AND. ii>1) then
                    call ran2(idum,xrand)
                    if (xrand < 0.5) then
                        delta_k = nint(ii*0.01)
                    else
                        delta_k = nint(-ii*0.01)
                    endif
                endif
                
                do while (iat(na) + delta_i < 0) 
                    delta_i = delta_i + nptg
                enddo
                do while (jat(na) + delta_j < 0) 
                    delta_j = delta_j + nptg
                enddo
                
                katnew = k_from(na) + delta_k
                
                if (katnew<0) then
                    katnew = 0
                endif
                
                if (katnew .ge. top - 4) then
                    check = 2
                    exit
                endif
                
                if (mod(katnew,4)==0) then                                          ! Layer 0
                    iatnew=mod(nptg+iat(na)+delta_i+mod(iat(na)+delta_i,2),nptg)              
                    if (mod(iatnew,4)==0) then                                ! if I=0,4,8,12 then J=0,4,8,12
                        jatnew=mod(nptg+jat(na)+delta_j-mod(jat(na)+delta_j,4),nptg)
                    else                                                        ! if I=2,6,10,14 then J=2,6,10 (k=0)
                        jatnew=mod(nptg+jat(na)+delta_j+2-mod(jat(na)+delta_j,4),nptg)
                    endif

                elseif(mod(katnew,4)==1) then                                       ! Layer 1
                    iatnew=mod(nptg+iat(na)+delta_i+1-mod(iat(na)+delta_i,2),nptg)
                    if (mod(iatnew,4)==1) then                                ! center atom I=1,5,9 J=1,5,9          
                        jatnew=mod(nptg+(jat(na)+delta_j+1-mod(jat(na)+delta_j,4)),nptg)
                    else                                                        ! center atom I=3,7,11 J=3,7,11
                        jatnew=mod(nptg+(jat(na)+delta_j+3-mod(jat(na)+delta_j,4)),nptg)
                    endif

                elseif(mod(katnew,4)==2) then                                       ! Layer 2
                    iatnew=mod(nptg+iat(na)+delta_i+mod(iat(na)+delta_i,2),nptg)
                    if (mod(iatnew,4)==0) then                                ! if I=0,4,8,12 then J=2,6,10 (k=2) 
                        jatnew=mod(nptg+jat(na)+delta_j+2-mod(jat(na)+delta_j,4),nptg)
                    else                                                        ! if I=2,6,10,14 then J=0,4,8 (k=2) 
                        jatnew=mod(nptg+jat(na)+delta_j-mod(jat(na)+delta_j,4),nptg)
                    endif
                    
                elseif(mod(katnew,4)==3) then                                       ! layer 3
                    iatnew=mod(nptg+iat(na)+delta_i+1-mod(iat(na)+delta_i,2),nptg)
                    if (mod(iatnew,4)==1) then                                ! center atom I=1,5,9 J=3,7,11
                        jatnew=mod(nptg+(jat(na)+delta_j+3-mod(jat(na)+delta_j,4)),nptg)
                    else                                                        ! center atom I=3,7,11 J=1,5,9
                        jatnew=mod(nptg+(jat(na)+delta_j+1-mod(jat(na)+delta_j,4)),nptg) 
                    endif
                endif  
                
                if (abs(iatnew-iat(na))+abs(jatnew-jat(na))+ABS(katnew-k_from(na))==0) then
                    goto 88
                endif
                
                if ((grid(iatnew,jatnew,katnew) .ne. 0)) then  
                    check = 1
                    
                    bind=-1
                    if (mod(katnew,2)==0) then                                      ! even numbers +1,+1, +1 and -1, -1, +1 and -1, +1, -1 and +1, -1, -1
                        par=+1
                    else
                        par=-1                                                  ! odd numbers -1,-1, -1 and +1, +1, -1 and -1, +1, +1 and +1, -1, +1
                    endif 
                    
                    locations = (/ grid(mod(nptg+iatnew+par,nptg),mod(nptg+jatnew+par,nptg),katnew+par),grid(mod(nptg+iatnew-par,nptg),mod(nptg+jatnew-par,nptg),katnew+par),grid(mod(nptg+iatnew-par,nptg),mod(nptg+jatnew+par,nptg),katnew-par),grid(mod(nptg+iatnew+par,nptg),mod(nptg+jatnew-par,nptg),katnew-par) /)
                    i_loc = (/ mod(nptg+iatnew+par,nptg), mod(nptg+iatnew-par,nptg), mod(nptg+iatnew-par,nptg), mod(nptg+iatnew+par,nptg) /)
                    j_loc = (/ mod(nptg+jatnew+par,nptg), mod(nptg+jatnew-par,nptg), mod(nptg+jatnew+par,nptg), mod(nptg+jatnew-par,nptg)  /)
                    k_loc = (/ katnew+par, katnew+par, katnew-par, katnew-par /)
                    shuffle = (/ 1, 2, 3, 4 /)
                    
                    do i = 4, 2, -1
                        call ran2(idum,xrand)
                        j = INT(xrand * i) + 1
                        temp = shuffle(i)          ! Swap elements
                        shuffle(i) = shuffle(j)
                        shuffle(j) = temp
                    end do
                    do i=1,4
                        if (locations(shuffle(i)) == 0 .OR. locations(shuffle(i))== na) then
                            if ((shuffle(i)<3 .and. katnew+par .ge. 0) .OR. (shuffle(i)>2 .AND. katnew-par .ge. 0)) then
                                move_i=i_loc(shuffle(i))-iat(na)
                                move_j=j_loc(shuffle(i))-jat(na)
                                move_k=k_loc(shuffle(i))-k_from(na)
                                bind=1    
                            endif
                        endif
                    enddo 
                    if (bind==-1) then
                        theta(na) = 0.0d0
                        phi(na) = 0.0d0
                        goto 37
                    endif
                    exit 
                    
                    
                elseif (grid(iatnew,jatnew,katnew) == 0) then                    
                    if (mod(katnew,2)==0) then 
                        par=+1                                                  ! even numbers +1,+1, +1 and -1, -1, +1 and -1, +1, -1 and +1, -1, -1
                    else
                        par=-1                                                  ! odd numbers -1,-1, -1 and +1, +1, -1 and -1, +1, +1 and +1, -1, +1
                    endif
! +1,+1, +1
                    if (grid(mod(nptg+iatnew+par,nptg),mod(nptg+jatnew+par,nptg),katnew+par)>0 .and. katnew+par .ge. 0) then
                        if (grid(mod(nptg+iatnew+par,nptg),mod(nptg+jatnew+par,nptg),katnew+par) .ne. na) then
                            move_i = iatnew - iat(na)
                            move_j = jatnew - jat(na)
                            move_k = katnew - k_from(na)
                            check = 1
                            exit
                        endif
!-1, -1, +1
                    elseif (grid(mod(nptg+iatnew-par,nptg),mod(nptg+jatnew-par,nptg),katnew+par)>0 .and. katnew+par .ge. 0) then
                        if (grid(mod(nptg+iatnew-par,nptg),mod(nptg+jatnew-par,nptg),katnew+par) .ne. na) then
                            move_i = iatnew - iat(na)
                            move_j = jatnew - jat(na)
                            move_k = katnew - k_from(na)
                            check = 1
                            exit
                        endif       
!-1, +1, -1
                    elseif (grid(mod(nptg+iatnew-par,nptg),mod(nptg+jatnew+par,nptg),katnew-par)>0  .and. katnew-par .ge. 0) then
                        if (grid(mod(nptg+iatnew-par,nptg),mod(nptg+jatnew+par,nptg),katnew-par) .ne. na) then
                            move_i = iatnew - iat(na)
                            move_j = jatnew - jat(na)
                            move_k = katnew - k_from(na)
                            check = 1
                            exit
                        endif      
! +1, -1, -1
                    elseif (grid(mod(nptg+iatnew+par,nptg),mod(nptg+jatnew-par,nptg),katnew-par)>0  .and. katnew-par .ge. 0) then
                        if (grid(mod(nptg+iatnew+par,nptg),mod(nptg+jatnew-par,nptg),katnew-par) .ne. na) then
                            move_i = iatnew - iat(na)
                            move_j = jatnew - jat(na)
                            move_k = katnew - k_from(na)
                            check = 1
                            exit
                        endif 
                        
                    elseif (katnew == 0) then
                        move_i = iatnew - iat(na)
                        move_j = jatnew - jat(na)
                        move_k = katnew - k_from(na)
                        check = 1
                        exit
                    endif
                endif
            enddo
            attempt_rate(na,2) = attempt_rate(na,2)+1
            if (check==1) then                                          !na is not on the surface
                evna = 1
                evt(1) = 1
                dina(na) = move_i                           !Assign new position
                djna(na) = move_j  
                dkna(na) = move_k 
                
            elseif (check==2) then                                      ! na is on the surface, it can evaporate
                nmolstay(kat(na))=nmolstay(kat(na))-1                   ! one less stay
                nmolevap(kat(na))=nmolevap(kat(na))+1                   ! one more evaporate  
                
                iatl = iat(na)
                jatl = jat(na)
                k_froml = k_from(na)
            
                if (dimer_weak(na) > 0) then
                    do i=1,4
                        if (dimer(dimer_weak(na),i) == na) then
                            dimer(dimer_weak(na),i) = 0.0d0
                            exit
                        endif
                    enddo
                    dimer_weak(na) = 0.0d0     
                endif
                                
                if (sum(dimer(na,:))>0) then
                    do i=1,4
                        if (dimer(na,i)>0) then
                            theta(dimer(na,i)) = theta(na)
                            phi(dimer(na,i)) = phi(na)
                            do q=1,na_max
                                if (num(q)==dimer(na,i)) then
                                    exit
                                endif
                            enddo
                            evt(q) = 5
                            tlist1(q) = minval(tlist1)-1d-8
                            tps = 0.0
                            
                            moved_dimer(i) = dimer(na,i)

                            dimer_weak(dimer(na,i)) = 0.0d0
                            dimer(na,i) = 0.0d0
                        endif
                    enddo
                endif 
                call remove(na,tlist1)                   ! na is removed from the list 
                
                nnlist(na)=0
                do di=-1,1                                                          !Update neighbours of old position
                    do dj=-1,1
                        do dk = -1,1
                            if ((k_froml+dk >=0) .AND. (k_froml+dk<=top-1)) then
                                if (ANY(moved_dimer==grid(mod(nptg+(iatl+di),nptg),mod(nptg+(jatl+dj),nptg),(k_froml+dk)))) then

                                elseif (grid(mod(nptg+(iatl+di),nptg),mod(nptg+(jatl+dj),nptg),(k_froml+dk))>0) then
                                    call add_one(grid(mod(nptg+(iatl+di),nptg),mod(nptg+(jatl+dj),nptg),(k_froml+dk)) ,tps,ev)
                                endif   
                            endif
                        enddo
                    enddo
                enddo
                
                call piksr3(na_max,tlist1,num,evt)
                na=num(1)
                ev=evt(1)
                evna=ev
                if (evna==5) then                                       ! if what is next is an evapoartion, go back to the beginning of the loop
                    event_count = event_count + 1
                    goto 37
                endif   
            endif
            
            theta(na) = 0.0d0
            phi(na) = 0.0d0
        endif
        
    ! ####################################################################
    ! ################### 9. Accretion (Ev=7) ############################
    ! ####################################################################
        call bouge(na,evna, tps)                      ! Call subroutine to move the species
       
    ! ####################################################################
    ! ####################### 9. Next event  #############################
    ! ####################################################################   

        call add_one(na, tps, ev)                     ! Add new event 

        if (na==(na_max) .and. iphase==0) then                                 ! Print end of phase 0
           print*,'fin accretion of',na, ' species in ',taccreti(na), '[s]'
           iphase=1                                                             ! update phase
           tlist0=taccreti(na)
           print*,'Entering TPD phase'
           Tg=Tgnext                                                            ! Set temperature to TPD temp
           print*,' --- '
        endif
    enddo
	print *, "JUMP> Bye !"
end program jump