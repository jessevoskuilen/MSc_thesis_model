
    
    ! ####################################################################
    ! ################# 1. Module: Physical data #########################
    ! ####################################################################
    ! ###################################F#################################
    
module phys_data

    implicit none                                                               !set program to not imply variable type 

    ! ####################################################################
    ! ########################## Constants ###############################
    ! ####################################################################

    integer, parameter        :: dp = SELECTED_REAL_KIND(P=15)
    
    real (kind=dp), parameter :: k_bz = 1.380658d-23                            ! Boltzmann constant
    real (kind=dp), parameter :: amu = 1.6605402d-27
    real (kind=dp), parameter :: a_pp =1.58d-10  !distance between 2 sites 
    real (kind=dp), parameter :: pi = 3.1415926535897932384_dp                  ! PI
    
    integer, parameter :: nspecies = 40
    
    character*8, public :: nameesp(0:40)                                        ! name of the species

    ! ############# mobility, evaporation and photon: rates ##############

    real (kind=dp)            :: evap(nspecies,0:nspecies) 	                    ! evaporation rate s-1
    real (kind=dp), public    :: kpp(nspecies,0:nspecies)                       ! mobility s-1
    real (kind=dp), public    :: k_accP(nspecies)                               ! accretion rate in phys
    real (kind=dp), public    :: k_accC(nspecies)                               ! accretion rates in chem

    ! ################### Grid and walk parameters #######################
  
    integer, public         :: na_max, iphase                                             ! number of atoms sent
    integer, public         :: nptg, nsite, top                                   ! size grid
    integer, public         :: bind  

    integer, public, dimension(:), allocatable      :: iat, jat, kat            ! position and nature of atoms
    integer, public, dimension(:), allocatable      :: iat0, jat0               ! previous position of atoms
    integer, public, dimension(:), allocatable      :: dina, djna, dkna         ! position and nature of atoms
    integer, public, dimension(:,:), allocatable    :: nb                       ! neighbors
    integer, public, dimension(:), allocatable      :: nn, nnlist               ! how many neigbours (1 to 4)
    integer, public, dimension(:,:,:), allocatable  :: grid                     ! grid
    integer, public, dimension(:,:,:), allocatable  :: empty                    ! empty grid
    integer, public, dimension(:), allocatable      :: evt, num, indice, ind
    integer, public, dimension(:), allocatable      :: k_from                   ! vertical position in the grid (from, to, original)
    double precision, dimension(:), allocatable     :: Eia 
    integer, allocatable                            :: attempt_rate(:,:), dimer(:,:),dimer_weak(:), last_move(:,:),hist_x(:,:),hist_y(:,:),hist_z(:,:) ! Track number of attempts to diffuse from self
    real (kind=dp), allocatable  :: theta(:), phi(:)

    real (kind=dp), public, dimension(:), allocatable   :: taccreti, tlist1
    real (kind=dp), public, dimension(:), allocatable   :: Ei,evan_list(:,:), tlast(:), t_evan(:,:)  		            ! matrice that contains the barriers to calculate the probabilities of reactions
 
    integer, public       :: ipos(27), jpos(27), kpos(27)


    ! ######## Physical characteristics of the gas and the grains ########

    real (kind=dp) :: Tg, Th, T0		                                        ! Grain and Gas temperatures
    real (kind=dp) :: nH2O, nAr, nKr, nXe, nCO2
    real (kind=dp) :: v(nspecies)		                                        ! velocity of atoms    
    real (kind=dp) :: stick            	                                        ! sticking coefficient
    real (kind=dp) :: sigma            	                                        ! cross section of the dust grains
    real (kind=dp) :: sigma_at                                                  ! cross section of a site


    ! ################## counting what is on the grain ###################

    integer        :: nmolstay(nspecies)                                        ! species that stay 
    integer        :: nmolevap(nspecies)                                        ! species that evaporate
    integer        :: nspec(nspecies)                                           ! sigma

    ! ################### Random numbers generator #######################

    real (kind=dp)  :: xrand                                                    ! random number 
    integer, public :: idum                                                     ! Random Number seed

    ! ########### Surface characteristics (defined in input.f) ###########

    real (kind=dp)  :: nu_p(nspecies)
    real (kind=dp)  :: tlist0


end module phys_data

    ! ####################################################################
    ! ##################### 2. Module: Utilities #########################
    ! ####################################################################
    
module utilities                                                                ! Provides the random number generator (From Numerical Recipes)

    use phys_data

    implicit none
    
    public ran2

    private

    contains

    ! ####################################################################
    ! #################### 2.1. Subroutine: Ran2 #########################
    ! ####################################################################

    subroutine ran2(idum, ran2_res)	                                            ! Creates a random number - From recipe, do not alter

        integer, intent (inout)         :: idum
        real (kind = dp), intent (out)  :: ran2_res

        integer, parameter :: im1 = 2147483563
        integer, parameter :: im2 = 2147483399
        integer, parameter :: imm1 = im1-1
        integer, parameter :: ia1 = 40014
        integer, parameter :: ia2 = 40692
        integer, parameter :: iq1 = 53668
        integer, parameter :: iq2 = 52774
        integer, parameter :: ir1 = 12211
        integer, parameter :: ir2 = 3791
        integer, parameter :: ntab = 32
        integer, parameter :: ndiv = 1+im1/ntab
   
        real (kind = dp), parameter :: am = 1.0_dp/im1
        real (kind = dp), parameter :: eps = 1.2e-12_dp
        real (kind = dp), parameter :: rnmx = 1.0e0_dp-eps

        integer, save :: idum2 = 123456789
        integer, save :: iy = 0
   
        integer, dimension (ntab), save :: iv = 0
   
        integer :: j, k

   
        if (idum <= 0) then
            idum = max(-idum,1)
            idum2 = idum
            do j=ntab+8,1,-1
                k = idum / iq1
                idum = ia1 * (idum-k*iq1) - k * ir1
                if (idum < 0) then
                    idum = idum + im1
                endif
                if (j <= ntab) then
                    iv(j) = idum
                endif
            end do
            iy = iv(1)
        endif

        k = idum / iq1
        idum = ia1 * (idum - k * iq1) - k * ir1

        if (idum < 0) then
           idum = idum + im1
        endif
        
        k = idum2 / iq2
        idum2 = ia2 * (idum2 - k * iq2) - k * ir2
        
        if (idum2 < 0) then
           idum2 = idum2 + im2
        endif
        
        j = 1 + iy / ndiv
        iy = iv(j) - idum2
        iv(j) = idum
        
        if (iy < 1) then
          iy = iy + imm1
        endif

        ran2_res = min(am*iy,rnmx)

    end subroutine ran2

end module utilities
    

    ! ####################################################################
    ! ######################## 3. Module: Marche #########################
    ! ####################################################################

module marche                                                                   ! Provides all elemental processes

    use phys_data
    use utilities
    
    implicit none
  
    public bouge, add_one, accretion, piksr3, remove, neighbors, counting, sample_gaussian

    private

    contains

    ! ####################################################################
    ! ##################### 3.1. Subroutine: Bouge #######################
    ! ####################################################################

	subroutine bouge (na, ev, tps)                    ! Moves an atom na from one site to the next, or accrete atom - Return info on new state
     
	    use phys_data

	    implicit none

	    integer         	    :: ev                                           ! Type of event, max height to accrete from
        double PRECISION        :: tps                                          ! Current elapsed time 
	    integer                 :: na, in_iat, in_jat                           ! Identity of moving atom
        real (kind = dp)        :: fix_angle, spread_angle, iangle, jangle, delta_i, delta_j
	    integer       		    :: di, dj, dk, ii,jj, par                       ! step
        integer                 :: dii,djj,dkk                                  ! Delta steps for neighbour updates
        integer                 :: floating_check
        integer                 :: temp, shuffle(4), locations(4), i,j, i_loc(4),j_loc(4), k_loc(4), moved_dimer(4),moved_dimer_loc(4,3)
        REAL (kind=dp)          :: r
        integer                 :: q,w

58      continue
    ! ####################################################################
    ! ######################## 3.1.1 Diffusion ###########################
    ! ####################################################################
        
        if (ev ==1) then                                                        ! Event is diffusion - Move in the ice
            moved_dimer = 9999999
                            
            grid(iat(na),jat(na),k_from(na))=0
            di=dina(na)
            dj=djna(na)
            dk=dkna(na)

            floating_check = 0
            do dii=-1,1                                                
                do djj=-1,1
                    do dkk = -1,1 
                        if (k_from(na)+dk+dkk .ge. 0 .AND. k_from(na)+dk+dkk<top-4) then
                            if (grid(mod(nptg+(iat(na)+di+dii),nptg),mod(nptg+(jat(na)+dj+djj),nptg),k_from(na)+dk+dkk) .ne. 0) then
                                if ((grid(mod(nptg+(iat(na)+di+dii),nptg),mod(nptg+(jat(na)+dj+djj),nptg),k_from(na)+dk+dkk) .ne. na)) then
                                    floating_check = 1
                                endif
                            endif
                        endif
                    enddo
                enddo
            enddo
            if (k_from(na)+dk==0) then
                floating_check=1
            endif


            if ((ABS(di)+ABS(dj)+ABS(dk)) .ne. 0) then
                if ((grid(mod(nptg+(iat(na)+di),nptg),mod(nptg+(jat(na)+dj),nptg),k_from(na)+dk) .ne. 0) .OR. (floating_check == 0)) then       ! If species already there, then do not go there (can happen as the add one time to move here was xx s ago and in the mean time other species could have moved
                    di=0
                    dj=0
                    dk=0
                endif
                if ((ABS(di)+ABS(dj)+ABS(dk)) .ne. 0) then   
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
                            if (dimer(na,i) > 0) then
                                if ((ABS(mod(nptg+(iat(na)+di),nptg)-iat(dimer(na,i))) > 1 .AND. ABS(mod(nptg+(iat(na)+di),nptg)-iat(dimer(na,i))) <39) .or. (ABS(mod(nptg+(jat(na)+dj),nptg)-jat(dimer(na,i))) > 1 .AND. ABS(mod(nptg+(jat(na)+dj),nptg)-jat(dimer(na,i))) <39) .OR.ABS(k_from(na)+dk - k_from(dimer(na,i))) > 1 ) then
                        
                                    if (mod(k_from(na)+dk,2)==0) then                                      ! even numbers +1,+1, +1 and -1, -1, +1 and -1, +1, -1 and +1, -1, -1
                                        par=+1
                                    else
                                        par=-1                                                  ! odd numbers -1,-1, -1 and +1, +1, -1 and -1, +1, +1 and +1, -1, +1
                                    endif 
                    
                                    locations = (/ grid(mod(nptg+iat(na)+di+par,nptg),mod(nptg+jat(na)+dj+par,nptg),k_from(na)+dk+par),grid(mod(nptg+iat(na)+di-par,nptg),mod(nptg+jat(na)+dj-par,nptg),k_from(na)+dk+par),grid(mod(nptg+iat(na)+di-par,nptg),mod(nptg+jat(na)+dj+par,nptg),k_from(na)+dk-par),grid(mod(nptg+iat(na)+di+par,nptg),mod(nptg+jat(na)+dj-par,nptg),k_from(na)+dk-par) /)
                                    i_loc = (/ mod(nptg+iat(na)+di+par,nptg), mod(nptg+iat(na)+di-par,nptg), mod(nptg+iat(na)+di-par,nptg), mod(nptg+iat(na)+di+par,nptg) /)
                                    j_loc = (/ mod(nptg+jat(na)+dj+par,nptg), mod(nptg+jat(na)+dj-par,nptg), mod(nptg+jat(na)+dj+par,nptg), mod(nptg+jat(na)+dj-par,nptg)  /)
                                    k_loc = (/ k_from(na)+dk+par, k_from(na)+dk+par, k_from(na)+dk-par, k_from(na)+dk-par /)
                                    shuffle = (/ 1, 2, 3, 4 /)
                    
                                    do q = 4, 2, -1
                                        call ran2(idum,xrand)
                                        j = INT(xrand * q) + 1
                                        temp = shuffle(q)          ! Swap elements
                                        shuffle(q) = shuffle(j)
                                        shuffle(j) = temp
                                    end do
                                    do q=1,4
                                        bind=-1
                                        if (i>1) then
                                            do w=1,4
                                                if (i_loc(shuffle(q))== moved_dimer_loc(w,1) .AND. j_loc(shuffle(q))== moved_dimer_loc(w,2) .AND. k_loc(shuffle(q))== moved_dimer_loc(w,3)) then
                                                    bind=-2
                                                endif
                                            enddo
                                        endif
                                        if (locations(shuffle(q)) == 0 .AND. bind==-1) then    
                                            if ((shuffle(q)<3 .and. k_from(na)+dk+par .ge. 0) .OR. (shuffle(q)>2 .AND. k_from(na)+dk-par .ge. 0)) then
                                                dina(dimer(na,i))=i_loc(shuffle(q))-iat(dimer(na,i))
                                                djna(dimer(na,i))=j_loc(shuffle(q))-jat(dimer(na,i))
                                                dkna(dimer(na,i))=k_loc(shuffle(q))-k_from(dimer(na,i))
                                                bind=1
                                                exit
                                            endif
                                        endif
                                    enddo        
                                  
                                    
                                    do q=1,na_max
                                        if (num(q)==dimer(na,i)) then
                                            exit
                                        endif
                                    enddo
                                    evt(q) = 1
                                    tlist1(q) = tlist1(1)-1d-8
                                    tps = 0.0
                                    
                                     if (bind<0) then                                              ! if place is occupied and no place around it is found, restart loop                       
                                        r = SQRT(real(di)**2 + real(dj)**2 + real(dk)**2)
                                        theta(dimer(na,i)) = ACOS(real(dk) / r)
                                        phi(dimer(na,i)) = ATAN2(real(dj), real(di))
                                        evt(q) = 5
                                     endif
                                    
                                    moved_dimer(i) = dimer(na,i)
                                    moved_dimer_loc(i,:) = (/ &
                                        iat(dimer(na,i)) + dina(dimer(na,i)), &
                                        jat(dimer(na,i)) + djna(dimer(na,i)), &
                                        k_from(dimer(na,i)) + dkna(dimer(na,i)) /)                                    
                                    dimer_weak(dimer(na,i)) = 0.0d0
                                    dimer(na,i) = 0.0d0
                                    
                                endif
                            endif
                        enddo
                    endif 
                endif
            endif

            if (k_from(na)<0) then
                print*, 'Error: k_from(na) < 0'
                stop
            endif

            
            iat(na) = mod(nptg+(iat(na)+di),nptg)                               ! update iat If the atome is out of the grain,
            jat(na) = mod(nptg+(jat(na)+dj),nptg)                               ! update jat put it inside (mirror)      
            k_from(na)=k_from(na)+dk                                            ! update ifrom
            grid(iat(na),jat(na),k_from(na))=na                                   ! Update grid data
            
            attempt_rate(na,1) = attempt_rate(na,1)+1                               ! Count number of attempts to move
            attempt_rate(na,5) = attempt_rate(na,5)+1
            
            if ((ABS(di)+ABS(dj)+ABS(dk)) .ne. 0) then 
                
                if (attempt_rate(na,4) == 1) then
                    attempt_rate(na,3:4) = 0
                endif

                last_move(na,1) = di
                last_move(na,2) = dj
                last_move(na,3) = dk
            
                hist_x(na,1:9) = hist_x(na,2:10)
                hist_y(na,1:9) = hist_y(na,2:10)
                hist_z(na,1:9) = hist_z(na,2:10) 
                
                hist_x(na,10) = iat(na)
                hist_y(na,10) = jat(na)
                hist_z(na,10) = k_from(na)
            endif

            do dii=-1,1                                                          !Update neighbours of old position
                do djj=-1,1
                    do dkk = -1,1
                        if ((k_from(na)-dk+dkk)>=0) then
                            if (grid(mod(nptg+(iat(na)-di+dii),nptg),mod(nptg+(jat(na)-dj+djj),nptg),(k_from(na)-dk+dkk))==na) then
                            elseif (ANY(moved_dimer==grid(mod(nptg+(iat(na)-di+dii),nptg),mod(nptg+(jat(na)-dj+djj),nptg),(k_from(na)-dk+dkk)))) then

                            elseif (grid(mod(nptg+(iat(na)-di+dii),nptg),mod(nptg+(jat(na)-dj+djj),nptg),(k_from(na)-dk+dkk))>0) then
                                call add_one(grid(mod(nptg+(iat(na)-di+dii),nptg),mod(nptg+(jat(na)-dj+djj),nptg),(k_from(na)-dk+dkk)) ,tps,ev)
                            endif
                        endif
                    enddo
                enddo
            enddo
            if ((ABS(di)+ABS(dj)+ABS(dk)) .ne. 0) then
                do dii=-1,1                                                          !Update neighbours of new position
                    do djj=-1,1
                        do dkk = -1,1
                            if ((k_from(na)+dkk)>=0) then
                                if (grid(mod(nptg+(iat(na)+dii),nptg),mod(nptg+(jat(na)+djj),nptg),(k_from(na)+dkk))==na) then
                                elseif (ANY(moved_dimer==grid(mod(nptg+(iat(na)+dii),nptg),mod(nptg+(jat(na)+djj),nptg),(k_from(na)+dkk)))) then

                                elseif (grid(mod(nptg+(iat(na)+dii),nptg),mod(nptg+(jat(na)+djj),nptg),(k_from(na)+dkk))>0) then
                                    call add_one(grid(mod(nptg+(iat(na)+dii),nptg),mod(nptg+(jat(na)+djj),nptg),(k_from(na)+dkk)) ,tps,ev)
                                endif
                            endif
                        enddo
                    enddo
                enddo 
            endif
        endif

    ! ####################################################################
    ! ######################## 3.1.2 Accretion ###########################
    ! ####################################################################
    
        if (ev == 7) then                                                       ! Accretion of na
39          continue

            tlast(na) = tlist1(1)
            
		    call ran2(idum, xrand)                                              ! random location of i is chosen
		    iat(na) = (nptg-1) * xrand
            
		    call ran2(idum, xrand)                                              ! random location of j is chosen
		    jat(na) = (nptg-1) * xrand 
            
		    k_from(na) = -2                                                     ! na comes from the gas phase

		    call ran2(idum, xrand)                                              ! here random di is picked up in case the species does not land on the right location
            if (xrand<0.33) then
                di=-1
            elseif (xrand<0.66) then
                di=0
            else
                di=1
            endif
		    
            call ran2(idum, xrand)                                              ! here random dj is picked up in case the species does not land on the right location
            if (xrand<0.33) then
                dj=-1
            elseif (xrand<0.66) then
                dj=0
            else
                dj=1
            endif
            
99          continue
            
            fix_angle = 0
            spread_angle = 85 - fix_angle
                
            call ran2(idum, xrand)                                              ! random angle of i descent between -80,80 degrees
		    iangle = (TAN((fix_angle + spread_angle - 2*spread_angle*xrand)*(pi/180)))
            call ran2(idum, xrand)                                              ! random angle of j descent between -80,80 degrees
            jangle = (TAN((fix_angle + spread_angle - 2*spread_angle*xrand)*(pi/180)))
           
            in_iat = iat(na)
            in_jat = jat(na)
            delta_i = 0
            delta_j = 0

            do ii=0,top-4                                                      ! species are coming from the top of the grid 
                
                jj=top-4-ii
                delta_i = delta_i + iangle
                delta_j = delta_j + jangle
                do while (in_jat + delta_j < 0) 
                    delta_j = delta_j + nptg
                enddo
                do while (in_iat + delta_i < 0) 
                    delta_i = delta_i + nptg
                enddo

                iat(na) = (mod(nptg+in_iat+NINT(delta_i),nptg))                       ! update iat
                jat(na) = (MOD(nptg+in_jat+NINT(delta_j),nptg))                       ! update jat
                
                if (mod(jj,4)==0) then                                          ! Layer 0
                    iat0(na)=mod(nptg+iat(na)+mod(iat(na),2),nptg)              
                    if (mod(iat0(na),4)==0) then                                ! if I=0,4,8,12 then J=0,4,8,12
                        jat0(na)=mod(nptg+jat(na)-mod(jat(na),4),nptg)
                    else                                                        ! if I=2,6,10,14 then J=2,6,10 (k=0)
                        jat0(na)=mod(nptg+jat(na)+2-mod(jat(na),4),nptg)
                    endif
                    
                elseif(mod(jj,4)==1) then                                       ! Layer 1
                    iat0(na)=mod(nptg+(iat(na)+1-mod(iat(na),2)),nptg)
                    if (mod(iat0(na),4)==1) then                                ! center atom I=1,5,9 J=1,5,9          
                        jat0(na)=mod(nptg+(jat(na)+1-mod(jat(na),4)),nptg)
                    else                                                        ! center atom I=3,7,11 J=3,7,11
                        jat0(na)=mod(nptg+(jat(na)+3-mod(jat(na),4)),nptg)
                    endif

                elseif(mod(jj,4)==2) then                                       ! Layer 2
                    iat0(na)=mod(nptg+iat(na)+mod(iat(na),2),nptg)
                    if (mod(iat0(na),4)==0) then                                ! if I=0,4,8,12 then J=2,6,10 (k=2) 
                        jat0(na)=mod(nptg+jat(na)+2-mod(jat(na),4),nptg)
                    else                                                        ! if I=2,6,10,14 then J=0,4,8 (k=2) 
                        jat0(na)=mod(nptg+jat(na)-mod(jat(na),4),nptg)
                    endif
                    
                elseif(mod(jj,4)==3) then                                       ! layer 3
                    iat0(na)=mod(nptg+(iat(na)+1-mod(iat(na),2)),nptg)
                    if (mod(iat0(na),4)==1) then                                ! center atom I=1,5,9 J=3,7,11
                        jat0(na)=mod(nptg+(jat(na)+3-mod(jat(na),4)),nptg)
                    else                                                        ! center atom I=3,7,11 J=1,5,9
                        jat0(na)=mod(nptg+(jat(na)+1-mod(jat(na),4)),nptg) 
                    endif
                endif
                
                bind=0      
                                                                                ! the location of na is found in iat and jat.
                if (grid(iat0(na),jat0(na),jj)==0) then                         ! now where will it go. If there is a spot available grid(iat,jat,jj)=0
                    if (mod(jj,2)==0) then 
                        par=+1                                                  ! even numbers +1,+1, +1 and -1, -1, +1 and -1, +1, -1 and +1, -1, -1
                    else
                        par=-1                                                  ! odd numbers -1,-1, -1 and +1, +1, -1 and -1, +1, +1 and +1, -1, +1
                    endif
                                                                                ! look around the empty spot to see if na is attaced to another species (it cannot float in space)
                    
! +1,+1, +1
                    if (grid(mod(nptg+iat0(na)+par,nptg),mod(nptg+jat0(na)+par,nptg),jj+par)>0 .and. jj+par .ge. 0) then
                        iat(na)=iat0(na)
                        jat(na)=jat0(na)
                        k_from(na)=jj
                        grid(iat(na),jat(na),k_from(na))=na
                        bind=1
                        exit         
!-1, -1, +1
                    elseif (grid(mod(nptg+iat0(na)-par,nptg),mod(nptg+jat0(na)-par,nptg),jj+par)>0 .and. jj+par .ge. 0) then
                        iat(na)=iat0(na)
                        jat(na)=jat0(na)
                        k_from(na)=jj
                        grid(iat(na),jat(na),k_from(na))=na
                        bind=1
                        exit         
!-1, +1, -1
                    elseif (grid(mod(nptg+iat0(na)-par,nptg),mod(nptg+jat0(na)+par,nptg),jj-par)>0  .and. jj-par .ge. 0) then
                        iat(na)=iat0(na)
                        jat(na)=jat0(na)
                        k_from(na)=jj
                        grid(iat(na),jat(na),k_from(na))=na
                        bind=1
                        exit        
! +1, -1, -1
                    elseif (grid(mod(nptg+iat0(na)+par,nptg),mod(nptg+jat0(na)-par,nptg),jj-par)>0  .and. jj-par .ge. 0) then
                        iat(na)=iat0(na)
                        jat(na)=jat0(na)
                        k_from(na)=jj
                        grid(iat(na),jat(na),k_from(na))=na
                        bind=1
                        exit  
!first layer
                    elseif(jj==0) then   
                        iat(na)=iat0(na)
                        jat(na)=jat0(na)
                        k_from(na)=jj
                        grid(iat(na),jat(na),k_from(na))=na
                        bind=1  
                        exit
                    endif

                elseif (grid(iat0(na),jat0(na),jj)>0) then                      ! if the place is occupied then get attached to this species
                    bind=-1
                    if (mod(jj,2)==0) then                                      ! even numbers +1,+1, +1 and -1, -1, +1 and -1, +1, -1 and +1, -1, -1
                        par=+1
                    else
                        par=-1                                                  ! odd numbers -1,-1, -1 and +1, +1, -1 and -1, +1, +1 and +1, -1, +1
                    endif 
                    
                    locations = (/ grid(mod(nptg+iat0(na)+par,nptg),mod(nptg+jat0(na)+par,nptg),jj+par),grid(mod(nptg+iat0(na)-par,nptg),mod(nptg+jat0(na)-par,nptg),jj+par),grid(mod(nptg+iat0(na)-par,nptg),mod(nptg+jat0(na)+par,nptg),jj-par),grid(mod(nptg+iat0(na)+par,nptg),mod(nptg+jat0(na)-par,nptg),jj-par) /)
                    i_loc = (/ mod(nptg+iat0(na)+par,nptg), mod(nptg+iat0(na)-par,nptg), mod(nptg+iat0(na)-par,nptg), mod(nptg+iat0(na)+par,nptg) /)
                    j_loc = (/ mod(nptg+jat0(na)+par,nptg), mod(nptg+jat0(na)-par,nptg), mod(nptg+jat0(na)+par,nptg), mod(nptg+jat0(na)-par,nptg)  /)
                    k_loc = (/ jj+par, jj+par, jj-par, jj-par /)
                    shuffle = (/ 1, 2, 3, 4 /)
                    
                    do i = 4, 2, -1
                        call ran2(idum,xrand)
                        j = INT(xrand * i) + 1
                        temp = shuffle(i)          ! Swap elements
                        shuffle(i) = shuffle(j)
                        shuffle(j) = temp
                    end do
                    do i=1,4
                        if (locations(shuffle(i)) == 0) then
                            if ((shuffle(i)<3 .and. jj+par .ge. 0) .OR. (shuffle(i)>2 .AND. jj-par .ge. 0)) then
                                iat(na)=i_loc(shuffle(i))
                                jat(na)=j_loc(shuffle(i))
                                k_from(na)=k_loc(shuffle(i))
                                grid(iat(na),jat(na),k_from(na))=na
                                bind=1
                                goto 67     
                            endif
                        endif
                    enddo        
                endif  
                if (bind==-1) then                                              ! if place is occupied and no place around it is found, restart loop                       
                    bind=0
                    goto 39
                endif
            enddo
67          continue
            nspec(kat(na))=nspec(kat(na))+1                                     ! Add one species to the count of species
            nmolstay(kat(na))=nmolstay(kat(na))+1                               ! Add one species to the count of species that stay
            
            do dii=-1,1                                                          !Update neighbours of new position
                do djj=-1,1
                    do dkk = -1,1
                        if ((k_from(na)+dkk)>=0) then
                            if (grid(mod(nptg+(iat(na)+dii),nptg),mod(nptg+(jat(na)+djj),nptg),(k_from(na)+dkk))>0) then
                                call add_one(grid(mod(nptg+(iat(na)+dii),nptg),mod(nptg+(jat(na)+djj),nptg),(k_from(na)+dkk)) ,tps,ev)
                            endif
                        endif
                    enddo
                enddo
            enddo 
        endif
   end subroutine bouge


    ! ####################################################################
    ! #################### 3.2. Subroutine: Add One ######################
    ! ####################################################################

    subroutine add_one(na, tps, ev)

        use phys_data
        use utilities

        implicit none
        
    ! ##################### Variable Definition ##########################
        
        integer                         :: na                                  ! Event number and atom identity
        double precision                :: tps                                  ! Current elapsed time
        integer, intent (out)           :: ev
        integer                         :: i, l, jj, di, dj, dk
        integer                         :: dii, djj, dkk                        ! Variables for hypothetical new position
        integer                         :: jjjj, jmove, wc, HE_n, HE_nf
        integer                         :: iatl, jatl, k_froml,nni,nnf,nniw,nniWd
        integer, allocatable     :: dfi(:,:), dfj(:,:), dfk(:,:)        
        integer                         :: nbposition
        integer, allocatable            :: nbourf(:,:),nbn(:,:)      
        integer                         :: highest_energy_neighbour, highest_energy_fn, p, q, revolve
        double precision                :: alpha, alpha_bulk, alpha_surface, evan, t00,deltaE, Ed, total_probability
        double precision, allocatable   :: kp(:), Efa(:,:), Eact(:),Es(:)
        double precision, allocatable   :: Ebin(:,:), increment(:)
        real                            :: result
 
        allocate(dfi(0:na_max,0:27), dfj(0:na_max,0:27), dfk(0:na_max,0:27))    
        allocate(nbourf(0:na_max,0:27),nbn(0:na_max,0:27))                      
        allocate(kp(27), Efa(0:na_max,0:27), Eact(27),Es(27))
        allocate(Ebin(0:nspecies,0:nspecies))      
        allocate(increment(0:nspecies))

 
    ! ################### Define Binding energies ########################

        alpha_surface = 0.43d0                                                    ! alpha is a coefficient that sets how the binding energy scales with diffusion
        alpha_bulk = 0.86d0
                
        increment(10) = 2550
        Ebin(10,10) = 2550                                                    ! H20-H20
        Ebin(10,36) = 2860
        Ebin(10,38) = 1695                                                 ! H20-Ar
        Ebin(10,39) = 2117                                                  ! H2O-Kr
        Ebin(10,40) = 2542                                                 ! H2O-xe
      
        Ebin(36,10) = 2860                                                    ! CO2-H20
        Ebin(36,36) = 2690                                                    ! CO2-CO2
        Ebin(36,38) = 1520                                                       ! CO2-Ar
        Ebin(36,39) = 1800                                                       ! CO2-Kr
        Ebin(36,40) = 2240                                                       ! CO2-Xe
        
        Ebin(38,10) = 1520!1695                                                      ! Ar-H20
        Ebin(38,36) = 1520                                                         ! Ar-CO2
        Ebin(38,38) = 875                                                       ! Ar-Ar
        Ebin(38,39) = 875                                                       ! Ar-Kr
        Ebin(38,40) = 875                                                       ! Ar-Xe

        Ebin(39,10) = 1800!2080                                                      ! Kr-H20
        Ebin(39,36) = 1800                                                         ! Kr-CO2
        Ebin(39,38) = 1160                                                      ! Kr-Ar
        Ebin(39,39) = 1160                                                      ! Kr-Kr
        Ebin(39,40) = 1160                                                      ! Kr-Xe
        
        Ebin(40,10) = 2240!2520                                                      ! Xe-H20
        Ebin(40,36) = 2240!0                                                         ! Xe-CO2
        Ebin(40,38) = 1600                                                      ! Xe-Ar
        Ebin(40,39) = 1600                                                      ! Xe-Kr
        Ebin(40,40) = 1600                                                      ! Xe-Xe

        nb(na,:)=0
       
        jmove=0
        bind=2
        nbposition=0
        
        do i=1,nspecies
            nu_p(10)=5.0d12 !H2O
            nu_p(36)=5.0d14 !CO2
            nu_p(38)=2.5d12 !Ar
            nu_p(39)=9.0d11 !Kr
            nu_p(40)=8.0d11 !Xe
        enddo

   !     if (iphase==0) then
  !          if (kat(na)==10) then
   !             Tg = 80
   !         else
  !              Tg = 20
  !          endif
  !      endif
    ! ############## Reset NA energies - Check position ##################
        if (k_from(na) .ge. 0) then                                             ! If -1 then it is evaporated,-2 is still of grid
            Eia(na)=0
            Efa(na,:)=0
      
    ! ###################### Check neighours #############################
            i=1                                                                 ! Dataset count
            l=1                                                                 ! Another count
            nni=0                                                               ! Number of neigbours
            nniW=0                                                              ! Number of water neighbours
            dfi(na,:)=0                                                         ! Reset neighbour positions
            dfj(na,:)=0
            dfk(na,:)=0

            do di=-1,1                                                          ! Check 
                do dj=-1,1      
                    do dk=-1,1
                        if (k_from(na)+dk.ge.0 .and. k_from(na)+dk < top-4 .and. (ABS(di)+ABS(dj)+ABS(dk)) .ne. 0) then	
                            nb(na,i)=grid(mod(nptg+iat(na)+di,nptg),mod(nptg+jat(na)+dj,nptg),k_from(na)+dk)
                        elseif (k_from(na)+dk < 0 .and. k_from(na) ==0) then
                            bind = 2
                            if(di==+1 .and. dj==+1 .and. dk== +1) then
                                bind=1
                            elseif (di==-1 .and. dj==-1 .and. dk== +1) then
                                bind=1
                            elseif (di==-1 .and. dj==+1 .and. dk== -1) then
                                bind=1
                            elseif (di==+1 .and. dj==-1 .and. dk== -1) then
                                bind=1
                            endif
                            if (bind==1) then
                                bind = 2
                                goto 22
                            endif
                        endif
                        if (nb(na,i)>0)  then                               ! If site is a species
                            nni=nni+1                                       ! Increase number of neighbours
                            bind = 0
                              
                            if (nni == 5) then
                                print*, 'error: nni>5'
                                print*, 'species : ', na, 'neighbours : ', nb(na,:)
                                print*, nb(na,1), ' : ', iat(nb(na,1)), jat(nb(na,1)), k_from(nb(na,1))
                                print*, nb(na,2), ' : ', iat(nb(na,2)), jat(nb(na,2)), k_from(nb(na,2))
                                print*, nb(na,3), ' : ', iat(nb(na,3)), jat(nb(na,3)), k_from(nb(na,3))
                                print*, nb(na,4), ' : ', iat(nb(na,4)), jat(nb(na,4)), k_from(nb(na,4))
                                print*, nb(na,5), ' : ', iat(nb(na,5)), jat(nb(na,5)), k_from(nb(na,5))
                                stop
                            endif

                            if (kat(nb(na,i))==10) then                     ! If water
                                nniW=nniW+1                                 ! Number of water neighbours +1
                            endif
                                                                            ! Check if it can revolve around its nieghbour to a more stable position
22                          continue                 
    
                            iatl=mod(nptg+iat(na)+di,nptg)              ! Define new location
                            jatl=mod(nptg+jat(na)+dj,nptg)
                            k_froml=k_from(na)+dk
                                                        
                            do dii=-1,1                                 ! If na would move there, would this place have neighbours?
                                do djj=-1,1 
                                    do dkk=-1,1 
                                        bind=0
                                        if (mod(k_froml,2)==0) then               ! check if position is a binding position, even numbers +1,+1, +1 and -1, -1, +1 and -1, +1, -1 and +1, -1, -1
                                            if(dii==+1 .and. djj==+1 .and. dkk== +1) then
                                                bind=1
                                            elseif (dii==-1 .and. djj==-1 .and. dkk== +1) then
                                                bind=1
                                            elseif (dii==-1 .and. djj==+1 .and. dkk== -1) then
                                                bind=1
                                            elseif (dii==+1 .and. djj==-1 .and. dkk== -1) then
                                                bind=1
                                            endif
                                        elseif (abs(mod(k_froml,2))==1) then           ! check if position is a binding position, odd numbers -1,-1, -1 and +1, +1, -1 and -1, +1, +1 and +1, -1, +1
                                            if(dii==-1 .and. djj==-1 .and. dkk== -1) then
                                                bind=1                          
                                            elseif (dii==+1 .and. djj==+1 .and. dkk== -1) then
                                                bind=1
                                            elseif (dii==+1 .and. djj==-1 .and. dkk== +1) then
                                                bind=1
                                            elseif (dii==-1 .and. djj==+1 .and. dkk== +1) then
                                                bind=1
                                            endif
                                        endif
                                        if (bind == 1) then
                                            if (k_froml+dkk .ge.0 .AND. k_froml + dkk < top-4  .and. (ABS(dii)+ABS(djj)+ABS(dkk)) .ne. 0) then  ! the new position is possible if there are some neighbors OR if it is on the surface
                                                if (grid(mod(nptg+iatl+dii,nptg),mod(nptg+jatl+djj,nptg),k_froml+dkk) == na) then       ! Original position will be empty so no neighbour @jazz if they were weakly bound, wouldnt they swap so there is somethingto be added here
                    
                                                elseif (grid(mod(nptg+iatl+dii,nptg),mod(nptg+jatl+djj,nptg),k_froml+dkk)==0) then
                                           !         if (ABS(dk+dkk)>0 .and. grid(mod(nptg+iat(na)-di,nptg),mod(nptg+jat(na)-dj,nptg),k_from(na)+dk)>0) then
                                            !        elseif (ABS(dk+dkk)==0 .AND. grid(mod(nptg+iat(na)-di,nptg),mod(nptg+jat(na)+dj,nptg),k_from(na)-dk)>0) then
                                             !       else
                                                        ipos(l)=mod(nptg+iatl+dii,nptg)   ! add possible position
                                                        jpos(l)=mod(nptg+jatl+djj,nptg)   ! add possible position
                                                        kpos(l)=k_froml+dkk               ! add possible position
                                                        dfi(na,l)=di+dii                        ! Store delta movement
                                                        dfj(na,l)=dj+djj                        ! Store delta movement
                                                        dfk(na,l)=dk+dkk                        ! Store delta movement
                                                        l=l+1                   ! Number of positions possible around na increased with 1 
                                             !       endif
                                                endif
                                            endif
                                        endif
                                    enddo
                                enddo
                            enddo
                                
                        elseif (nb(na,i)==0) then !.or. Ebin(kat(na),10)>Ebin(kat(nb(na,i)),10))  then                             ! there is no neighbour there, or it is more weakly bound. So it is a position where na could move
                            bind=0
                            if (mod(k_from(na),2)==0) then           ! check if position is a binding position, even numbers +1,+1, +1 and -1, -1, +1 and -1, +1, -1 and +1, -1, -1         
                                if(di==+1 .and. dj==+1 .and. dk== +1) then
                                    bind=1
                                elseif (di==-1 .and. dj==-1 .and. dk== +1) then
                                    bind=1
                                elseif (di==-1 .and. dj==+1 .and. dk== -1) then
                                    bind=1
                                elseif (di==+1 .and. dj==-1 .and. dk== -1) then
                                    bind=1
                                endif
                            elseif (mod(k_from(na),2)==1) then    ! check if position is a binding position, odd numbers -1,-1, -1 and +1, +1, -1 and -1, +1, +1 and +1, -1, +1
                                if(di==-1 .and. dj==-1 .and. dk== -1) then
                                    bind=1                          
                                elseif (di==+1 .and. dj==+1 .and. dk== -1) then
                                    bind=1
                                elseif (di==+1 .and. dj==-1 .and. dk== +1) then
                                    bind=1
                                elseif (di==-1 .and. dj==+1 .and. dk== +1) then
                                    bind=1
                                endif
                            endif
 
                            if (bind==1 .AND. k_from(na)<top-4) then                               ! here the species can go to this place di, dj, dk
                                iatl=mod(nptg+iat(na)+di,nptg)              ! Define new location
                                jatl=mod(nptg+jat(na)+dj,nptg)
                                k_froml=k_from(na)+dk
                                do dii=-1,1                                 ! If na would move there, would this place have neighbours?
                                    do djj=-1,1 
                                        do dkk=-1,1 
                                            if (k_froml+dkk .ge.0 .AND. k_froml + dkk < top-4  .and. (ABS(dii)+ABS(djj)+ABS(dkk)) .ne. 0) then  ! the new position is possible if there are some neighbors OR if it is on the surface 
                                                if (grid(mod(nptg+iatl+dii,nptg),mod(nptg+jatl+djj,nptg),k_froml+dkk) == na) then       ! Original position will be empty so no neighbour @jazz if they were weakly bound, wouldnt they swap so there is somethingto be added here
                                                        
                                                elseif (grid(mod(nptg+iatl+dii,nptg),mod(nptg+jatl+djj,nptg),k_froml+dkk)>0 .or. k_froml==0) then 
                                                    ipos(l)=mod(nptg+iat(na)+di,nptg)   ! add possible position
                                                    jpos(l)=mod(nptg+jat(na)+dj,nptg)   ! add possible position
                                                    kpos(l)=k_from(na)+dk               ! add possible position
                                                    dfi(na,l)=di                        ! Store delta movement
                                                    dfj(na,l)=dj                        ! Store delta movement
                                                    dfk(na,l)=dk                        ! Store delta movement
                                                    l=l+1                   ! Number of positions possible around na increased with 1 
                                                    goto 53 ! new position has neighbour and is stored, so no need to continue search
                                       
                                                endif
                                            endif
                                        enddo
                                    enddo
                                enddo
                            endif
                        endif
                        if (nb(na,i)>0) then
                                i=i+1 
                        endif
    53                      continue

                    enddo
                enddo
            enddo
        endif

        nnlist(na)=nni
        
        if (nni==1 .AND. k_from(na) > 0) then
            dimer_weak(na) = nb(na,1)
            do i=1,4
                if (dimer(dimer_weak(na),i) == na) then
                    goto 90
                endif
            enddo
            do i=1,4
                if (dimer(dimer_weak(na),i) == 0) then
                    dimer(dimer_weak(na),i) = na
                    exit
                endif
            enddo
90 continue
        elseif (nni>1 .OR. k_from(na)==0) then
            if (dimer_weak(na) > 0) then
                do i=1,4
                    if (dimer(dimer_weak(na),i) == na) then
                        dimer(dimer_weak(na),i) = 0.0d0
                        exit
                    endif
                enddo
                dimer_weak(na) = 0.0d0     
            endif
        endif
        


    ! ################# Calculate initial energies #######################
        highest_energy_neighbour = 0
        HE_n = 0
        do p=1,27
            if (nb(na,p)==0) then
                
            else
                if (Ebin(kat(na),kat((nb(na,p)))) > highest_energy_neighbour) then
                    highest_energy_neighbour = Ebin(kat(na),kat((nb(na,p))))
                    HE_n = kat((nb(na,p)))
                endif
            endif
        enddo
        if (k_from(na)==0) then                                             ! If species is at the surface then         
            Eia(na)=(Ebin(kat(na),10))                                        ! Bound to water
            if (kat(na)==10) then
                Eia(na)=Ebin(10,10)+increment(10)*0.8
                do i=1,nni
                    if ((kat(nb(na,i)) .ne. 10) .AND. (kat(nb(na,i)) > 0)) then
                        Eia(na) = Eia(na) + Ebin(kat(na),kat(nb(na,i)))
                    elseif (kat(nb(na,i))==10) then
                        Eia(na) = Eia(na) + increment(10)*0.8**(i+1)
                    endif
                enddo
            endif
            
        elseif (nni==0) then
            Eia(na)=0
            
        else
            if (kat(na) == 10) then
                if (nniW>0) then
                    wc = 0
                    do i=1,nni
                        if (wc > 0 .AND. kat(nb(na,i))==10) then
                            Eia(na) = Eia(na) + increment(10)*0.8**wc
                        else
                            Eia(na) = Eia(na) +  Ebin(kat(na),kat(nb(na,i)))
                        endif
                        if (kat(nb(na,i))==10) then
                            wc = wc + 1
                        endif
                    enddo
                else
                    Eia(na) = Ebin(HE_n,kat(na))
                endif
            else
                Eia(na) = Ebin(kat(na),HE_n)
            endif
        endif
        call sample_gaussian(1.0,0.05,2.0,result)
        eia(na) = eia(na)*result
                
        if (Eia(na) == 0) then                                              ! If the species is floating in space and has no energy, relocate     
            do di=1,na_max
                if (num(di)==na) then
                    exit
                endif
            enddo
            evt(di) = 5
            tlist1(di) = tlist1(1)
            tps = 0.0
            goto 44
        endif
        
    ! ############ Second loop -- Ef -- final binding energy ############# 
        nbposition=l-1                                                      ! number of positions possible around na
        do l=1,nbposition                                                   ! Loop through all possible positions
            jj=1 
            Efa(na,l)=0                                                     ! Initialise
            nniWd=0
            nnf=0
            nbourf(na,:)=0
     
            do dii=-1,1                                                     ! Loop through all neighbours of new position
                do djj=-1,1 
                    do dkk=-1,1
                        if (abs(dii)+abs(djj)+abs(dkk)==0) then             ! If self then do not consider
                            
                        elseif (kpos(l)+dkk.ge.0 .and. kpos(l)+dkk <top-1) then                  ! If on the grid         
                            nbourf(na,jj)=grid(mod(nptg+ipos(l)+dii,nptg),mod(nptg+jpos(l)+djj,nptg),kpos(l)+dkk)   ! neighbor of the possible future location of na at +dii,+djj,+dkk 
                            if (nbourf(na,jj)>0 .and. nbourf(na,jj).ne.na) then                                     ! Do I have neighbors in the new place-- > yes I can move there
                                nnf=nnf+1                                   ! Increase final nn
                                if (kat(nbourf(na,jj))==10) then            ! If neighbour is water count
                                    nniWd=nniWd+1
                                endif            
                                jj=jj+1                                             ! Increase loop to next neighbour
                            endif
                        endif
                    enddo
                enddo
            enddo
            nbn(na,l)=nnf
            highest_energy_fn = 0
            HE_nf = 0
            do q=1,27
                if (nbourf(na,q)==0) then
                
                else
                    if (Ebin(kat(na),kat(nbourf(na,q))) > highest_energy_fn) then
                        highest_energy_fn = Ebin(kat(na),kat(nbourf(na,q)))
                        HE_nf = kat(nbourf(na,q))
                    endif
                endif
            enddo
        
            if (kpos(l)==0) then                                            ! If new pos is at surface, bound to water
                Efa(na,l)=(Ebin(kat(na),10))                                ! assume two bounds with water, two free
                if (kat(na)==10) then
                    Efa(na,l)=Ebin(10,10)+increment(10)*0.8
                    do i=1,nnf
                        if ((kat(nbourf(na,i)) .ne. 10) .AND. (kat(nbourf(na,i)) > 0)) then
                            Efa(na,l) = Efa(na,l) + Ebin(kat(na),kat(nbourf(na,i)))
                        elseif (kat(nbourf(na,i))==10) then
                            Efa(na,l) = Efa(na,l) + increment(10)*0.8**(i+1)
                        endif
                    enddo
                endif

            else
                if (kat(na)==10) then                                       ! If species is water
                    if (nniWd>0) then  
                        wc = 0
                        do i=1,nnf
                            if (wc > 0 .AND. kat(nbourf(na,i))==10) then
                                 Efa(na,l) = Efa(na,l) + increment(10)*0.8**wc
                            else
                                Efa(na,l) = Efa(na,l) +  Ebin(kat(na),kat(nbourf(na,i)))
                            endif
                            if (kat(nbourf(na,i))==10) then
                                wc = wc + 1
                            endif
                        enddo
                    else
                        Efa(na,l) = Ebin(kat(na), HE_nf)
                    endif
                else
                    efa(na,l) = Ebin(kat(na),HE_nf)
                endif
            endif
            call sample_gaussian(1.0,0.05,2.0,result)
            efa(na,l) = efa(na,l)*result
        enddo
        
        kp=1d-99
        
        if (nni .ge. 3) then
            alpha=alpha_bulk
        else
            alpha=alpha_surface
        endif

        do jjjj=1,nbposition                                                ! Loop through all neighbours of position l  
            
            if (Efa(na,jjjj)==0) then                                       ! If final energy is 0 make it near 0
                Efa(na,jjjj)=1d-99
            endif
            deltaE=(max(Eia(na),Efa(na,jjjj))-min(Eia(na),Efa(na,jjjj)))    ! Calculate the energy difference between initial and final energy
            Ed=alpha*min(Eia(na),efa(na,jjjj))                             ! here diffusion barrier is alpha*Eb

            if (Efa(na,jjjj) > Eia(na)) then                                ! If final energy is larger
                Eact(jjjj)=Ed                                               ! Activation barrier is the diffusion energy
                Es(jjjj)=Eia(na)-Ed                                         ! Saddle point energy is initial - diffusion
            else
                Eact(jjjj)=Ed+deltaE                                        ! Activation energy is diffusion plus the increase
                Es(jjjj)=Efa(na,jjjj)-Ed                                    ! Saddle point energy is final - diffusion
            endif
            if (Efa(na,jjjj)-Es(jjjj) >0) then                              ! If the final energy is more than the saddle point energy 
                kp(jjjj)=4*(sqrt((Efa(na,jjjj)-Es(jjjj))/(Eia(na)-Es(jjjj)))/(1+sqrt((Efa(na,jjjj)-Es(jjjj))/(Eia(na)-Es(jjjj))))**2)*nu_p(kat(na))*exp(-(Eact(jjjj))/(Tg))
            else
                kp(jjjj)=1d-99                                              ! Otherwise near zero
            endif
            if (kp(jjjj)<1d-99) then                                        ! set to near zero again @jazz can this happen? think it is not used probably
                kp(jjjj)=1d-99
            endif  
        enddo

        kp(nbposition+1)=nu_p(kat(na))*exp(-Eia(na)/(Tg))   
   
        evan_list(na,1) = 1/kp(nbposition+1)
        t00=1d99

        total_probability = SUM(kp(1:nbposition+1))
        call ran2(idum,xrand)
        
        do jj=1,nbposition+1                                               ! Find first event  to occur
            if (sum(kp(1:jj))<total_probability*xrand) then
            
            else
                jmove=jj
                call ran2(idum,xrand)
                t00=1 / maxval(kp(1:nbposition+1))
                exit
            endif       
        enddo
        
        tps=t00
                
        if (nni == 4) then
            tps = tps + 1d6
        endif
   
        if ((jmove .le. nbposition).and.(jmove>0)) then                        !   
            ev=1
            dina(na)=dfi(na,jmove)
            djna(na)=dfj(na,jmove)
            dkna(na)=dfk(na,jmove)
            if (jmove>27) then                            ! If not within 1-27 something went wrong so stop
                print*, 'jmove>27'
                stop
            endif 
        elseif (jmove==nbposition+1) then                               ! evaporation event
            ev = 5 
            dina(na)=0
            djna(na)=0
            dkna(na)=0
        elseif (jmove==0) then                                          ! If no move then time is infinite, go to end of list
            dina(na)=0
            djna(na)=0
            dkna(na)=0
            ev=5 
            tps=1d95
        endif
        
        if (ev==1) then
            if (ABS(dina(na)+last_move(na,1))+ABS(djna(na)+last_move(na,2))+ABS(dkna(na)+last_move(na,3))==0) then
                attempt_rate(na,3) = attempt_rate(na,3) + 1
            else
                attempt_rate(na,4) = 1
            endif     
        endif
        
        if (attempt_rate(na,3) > 3) then
            ev = 5 
            dina(na)=0
            djna(na)=0
            dkna(na)=0
            tps = 1 / kp(nbposition+1)
        endif
   
        do di=1,na_max
            if (num(di)==na) then
                exit
            endif
        enddo
                
        tlist1(di) = tlist1(1)+tps  

        evt(di) = ev     

44      continue  
 !   if (iphase==0) then
 !       Tg = 20
 !   endif
    end subroutine add_one


    ! ####################################################################
    ! ################### 3.3. Subroutine: Accretion #####################
    ! ####################################################################

    subroutine accretion

        use phys_data
        use utilities

        implicit none
        
        integer                                     :: j, i			            ! Counters for loops : i->species; j->atoms
        integer                                     :: minpos(1)		        ! Position of the minimal time    
        double precision                            :: tacc    	                ! Total time of accretion
        double precision, dimension(:), allocatable :: tmps                     ! Times
  
    
        allocate (tmps(2*nspecies))	
        allocate (taccreti(0:na_max))	                                        ! Individual accretion time allocation
        allocate (indice(na_max))		                                        ! Final index allocation
        allocate (ind(na_max))		                                            ! Initial index allocation

        tacc=0.d0				                                                ! Init of tacc and taccreti
        taccreti(0)=0.d0			                                            !

        print*,"JUMP> filling taccret.out ..."
    
	    do j=1, na_max					                                        ! Loop over all arriving atoms
		    do i=1,nspecies				                                        ! Loop over each species
                if (k_accP(i)>1d90) then 
                    tmps(2*i-1)=1d99
                else
                    call ran2(idum,xrand)			                            ! Create a random double prec
                    tmps(2*i-1)=-log(1.0_dp-xrand)/k_accP(i)	                ! Time to accrete in a physisorbed site
                    ind(2*i-1)=2*i-1                                            ! Index is odd for physisorption
                endif
       	  
                if (k_accC(i)>1d90) then
                    tmps(2*i)=1d99
                else
                    call ran2(idum,xrand)			                            ! Create a random double prec
                    tmps(2*i)=-log(1.0_dp-xrand)/k_accC(i)	                    ! Time to accrete in a chemisorbed site
                    ind(2*i)=2*i                                                ! Index is even for chemi
                endif
		    enddo
	        minpos = minloc(tmps)				                                ! Minpos is equal to the position of the minimal time
	        indice(j) = ind(minpos(1))			                                ! Setting the corresponding index
	        taccreti(j) = taccreti(j-1)+tmps(minpos(1))	                        ! Accretion time
	        write(45,*)j,ind(minpos(1)),tmps(minpos(1)) 	                    ! Writing to taccret.out
	    enddo 

        close(45)
        tacc=taccreti(na_max)				                                    ! Total accretion time
    
        print*,"JUMP> Total accretion time : ",tacc,'number',na_max
        deallocate(tmps)
    end subroutine accretion

    ! ####################################################################
    ! ##################### 3.4. Subroutine: piksr3 ######################
    ! ####################################################################

    subroutine piksr3(n,arr,brr,crr)                                        ! Subroutine to order list
        integer         :: n
        REAL (kind=dp)  :: arr(n), a
        integer         :: brr(n), b, crr(n), c
        integer         :: i,j
      
        do 12 j=2,n
            a=arr(j)
            b=brr(j)
            c=crr(j)
            do 11 i=j-1,1,-1
                if(arr(i).le.a) goto 10
                arr(i+1)=arr(i)
                brr(i+1)=brr(i)
	            crr(i+1)=crr(i)
11          continue
            i=0
10          arr(i+1)=a
            brr(i+1)=b
	        crr(i+1)=c
12      continue
        return
    end subroutine piksr3

    ! ####################################################################
    ! ##################### 3.5. Subroutine: remove ######################
    ! ####################################################################

	subroutine remove(na,tlist1)                                     ! Subroutine to remove species
        
        integer, intent(in)		                        :: na
        integer                                         :: i,j 
        double precision,intent(inout), allocatable    :: tlist1(:)
        
        grid(iat(na),jat(na),k_from(na))=0

        iat(na)=-1
        jat(na)=-1
        k_from(na)=-2

        tlist1(1)=1d99
        do i=1,4
            if (nb(na,i)>0) then
                do j=1,4
                    if (nb(nb(na,i),j)==na) then
                        nb(nb(na,i),j)=0
                        nn(nb(na,i))=nn(nb(na,i))-1
                    endif
                enddo
            endif
        enddo
       
	end subroutine remove
         
    ! ####################################################################
    ! ################### 3.7. Subroutine: Neighbours ####################
    ! ####################################################################

    subroutine neighbors(na,nnb1,nn1,nbnb,nnnb)

        integer, intent (in)                :: na
        integer, intent (out), allocatable  :: nnb1(:,:),nbnb(:,:),nnnb(:)
        integer                             :: n0,na1,na2,na3,na4,i,j,par,nn1
        allocate(nnb1(0:na_max,0:4),nbnb(0:na_max,0:4),nnnb(0:4))
        na1=0
        na2=0
        na3=0
        na4=0

          if (mod(k_from(na),2)==0) then                                        
             par=+1                                                             ! even numbers +1,+1, +1 and -1, -1, +1 and -1, +1, -1 and +1, -1, -1
          else
             par=-1                                                             ! odd numbers -1,-1, -1 and +1, +1, -1 and -1, +1, +1 and +1, -1, +1
          endif

        if (k_from(na)+par .ge.0) then
           na1=grid(mod(nptg+iat(na)+par,nptg),mod(nptg+jat(na)+par,nptg),k_from(na)+par)
           na2=grid(mod(nptg+iat(na)-par,nptg),mod(nptg+jat(na)-par,nptg),k_from(na)+par)
        else
           na1=0
           na2=0
        endif
        if (k_from(na)-par .ge. 0) then
           na3=grid(mod(nptg+iat(na)-par,nptg),mod(nptg+jat(na)+par,nptg),k_from(na)-par)
           na4=grid(mod(nptg+iat(na)+par,nptg),mod(nptg+jat(na)-par,nptg),k_from(na)-par)
        else
           na3=0
           na4=0
        endif
        
        nn(na)=0
        nb(na,1)=na1
        nb(na,2)=na2
        nb(na,3)=na3
        nb(na,4)=na4
        nnb1(na,1)=na1
        nnb1(na,2)=na2
        nnb1(na,3)=na3
        nnb1(na,4)=na4

        do i=1,4
           if (nb(na,i)>0) then
              nn(na)=nn(na)+1
           endif
        enddo
        nn1=nn(na)

        do i=1,4                                                                ! now we calculate the neoghbors of the neighbors
            if (nb(na,i)>0) then
               n0=nb(na,i)
               if (mod(k_from(nb(na,i)),2)==0) then
                  par=+1                                                        ! even numbers +1,+1, +1 and -1, -1, +1 and -1, +1, -1 and +1, -1, -1
               else
                  par=-1                                                        ! odd numbers -1,-1, -1 and +1, +1, -1 and -1, +1, +1 and +1, -1, +1
               endif

                if (k_from(n0)+par .ge.0) then
                   na1=grid(mod(nptg+iat(n0)+par,nptg),mod(nptg+jat(n0)+par,nptg),k_from(n0)+par)
                   na2=grid(mod(nptg+iat(n0)-par,nptg),mod(nptg+jat(n0)-par,nptg),k_from(n0)+par)
                else
                   na1=0
                   na2=0
                endif
                if (k_from(n0)-par .ge. 0) then
                   na3=grid(mod(nptg+iat(n0)-par,nptg),mod(nptg+jat(n0)+par,nptg),k_from(n0)-par)
                   na4=grid(mod(nptg+iat(n0)+par,nptg),mod(nptg+jat(n0)-par,nptg),k_from(n0)-par)
                else
                   na3=0
                   na4=0
                endif
            endif

            nn((nb(na,i)))=0
            nb((nb(na,i)),1)=na1
            nb((nb(na,i)),2)=na2
            nb((nb(na,i)),3)=na3
            nb((nb(na,i)),4)=na4        
        enddo

        do i=1,4
           do j=1,4
              if (nb((nb(na,i)),j)>0) then
                 nbnb(nb(na,i),j)=(nb((nb(na,i)),j))
                 nn((nb(na,i)))=nn((nb(na,i)))+1
              endif
           enddo
        enddo

        nnnb(1)=nn((nb(na,1)))
        nnnb(2)=nn((nb(na,2)))
        nnnb(3)=nn((nb(na,3)))
        nnnb(4)=nn((nb(na,4)))

    end subroutine neighbors

    ! ####################################################################
    ! #################### 3.8. Subroutine: COUNTING #####################
    ! ####################################################################
        
    subroutine counting(nempty,surftot,thick,nmol)
        integer, intent (out)           :: nempty,surftot,nmol
        double precision, intent (out)  :: thick
        integer, allocatable            :: nnb1(:,:),nbnb(:,:),nnnb(:) 
        integer                         :: ii0,jj0, ii, iikkk, kkk, kkk0, ij, il, jj, kk, nn1
        integer                         :: bind, ibind, par0, ipar, jpar, kpar, surf, surf0
        integer                         :: ne0, ne1,ne2, ne3, ne4, ne5, ne6, ne7, ne8, nb0, nb1, nb2, nb3, nb4
        allocate(nnb1(0:na_max,0:4),nbnb(0:na_max,0:4),nnnb(0:4))
        thick=0
        iikkk=0
        bind=0
        nb0=0
        nb1=0
        nb2=0
        nb3=0
        nb4=0
        ne0=0
        ne1=0
        ne2=0
        ne3=0
        ne4=0          
        ne5=0
        ne6=0
        ne7=0
        ne8=0
        ii0=0
        jj0=0
        il=0
        empty=0
        nempty=0
        surftot=0
        nmol=0
        nn=0

        do ii=0,nptg-1
        do jj=0,nptg-1
            kkk0=0
            do kk=3,top
                kkk=top-kk
                if (grid(ii,jj,kkk)>0 .and. ii+jj .ne. ii0+jj0) then  
                    thick=thick+kkk
                    kkk0=kkk
                    ii0=ii
                    jj0=jj
                    if (nn(grid(ii,jj,kkk))>4) then
                        print*,'more than 4 neighbours'
                        stop
                    endif
                endif
                
                if (grid(ii,jj,kkk)>0) then
                    call neighbors(grid(ii,jj,kkk),nnb1,nn1,nbnb,nnnb)                              
                    nn(grid(ii,jj,kkk))=nn1
                    nmol=nmol+1
                    if (nn(grid(ii,jj,kkk))==1) then
                        nb1=nb1+1
                    elseif (nn(grid(ii,jj,kkk))==2) then
                        nb2=nb2+1
                    elseif (nn(grid(ii,jj,kkk))==3) then
                        nb3=nb3+1
                    elseif (nn(grid(ii,jj,kkk))==4) then
                        nb4=nb4+1
                    elseif (nn(grid(ii,jj,kkk))==0 .and. kkk<kkk0) then
                        nb0=nb0+1
                    endif
                    
                elseif (grid(ii,jj,kkk)==0 .and. kkk<kkk0) then
                    ibind=0
                    if (mod(kkk,4)==0) then
                        if (mod(jj,4)==0 .and. mod(ii,4)==0) then
                            ibind=2
                        elseif (mod(jj,4)==2 .and. mod(ii,4)==2) then
                            ibind=2
                        endif
                    endif
                    
                    if (mod(kkk,4)==1) then
                        if (mod(jj,4)==1 .and. mod(ii,4)==1) then
                            ibind=2
                        elseif (mod(jj,4)==3 .and. mod(ii,4)==3) then
                            ibind=2
                        endif
                    endif 
                    
                    if (mod(kkk,4)==3) then
                        if (mod(jj,4)==1 .and. mod(ii,4)==3) then
                            ibind=2
                        elseif (mod(jj,4)==3 .and. mod(ii,4)==1) then
                            ibind=2
                        endif
                    endif
                    
                    if (mod(kkk,4)==2) then
                        if (mod(ii,4)==0 .and. mod(jj,4)==2) then
                            ibind=2
                        elseif (mod(jj,4)==0 .and. mod(ii,4)==2) then
                            ibind=2                             
                        endif
                    endif  
                    
                    if (ibind==2) then
                        nempty=nempty+1
                        par0=15
99                      continue
                        par0=par0-1
                        do  ipar=-par0,par0
                            do  jpar=-par0,par0
                                do kpar=-par0,par0
                                    if (kkk+kpar .ge.0 .and. kkk+kpar < top .AND. par0 .ge. 1) then
                                        if (grid(mod(nptg+ii+ipar,nptg),mod(nptg+jj+jpar,nptg),kkk+kpar)>0 .and. abs(ipar)+abs(jpar)+abs(kpar)>0) then   
                                            goto 99
                                        endif
                                   elseif(par0 < 1) then
                                      ibind=1
                                   endif
                                enddo
                            enddo
                        enddo
                        
                        empty(ii,jj,kkk)=par0
                        surf0=0
                        surf=0
                        ij=0
                        if (par0==0) then
                            ne0=ne0+1
                            do  ipar=-1,1
                                do  jpar=-1,1
                                    do kpar=-1,1
                                        if (grid(mod(nptg+ii+ipar,nptg),mod(nptg+jj+jpar,nptg),kkk+kpar)>0) then
                                            surf=surf+1
                                            surf0=surf
                                        endif
                                    enddo
                                enddo
                            enddo
                            surftot=surftot+surf
                            elseif (par0==1) then
                                ne1=ne1+1
                            elseif (par0==2) then
                                ne2=ne2+1
                            elseif (par0==3) then
                                ne3=ne3+1
                            elseif (par0==4) then
                                ne4=ne4+1
                            elseif (par0==5) then
                                ne5=ne5+1
                            elseif (par0==6) then
                                ne6=ne6+1
                            elseif (par0==7) then
                                ne7=ne7+1
                            elseif (par0==8) then
                                ne8=ne8+1
                
                            endif
                        endif
                    endif
                enddo
            enddo
        enddo
    end subroutine counting

    ! ####################################################################
    ! #################### 3.8. Subroutine: COUNTING #####################
    ! ####################################################################
    SUBROUTINE sample_gaussian(mean, stddev, maxdev, result)
        IMPLICIT NONE
        REAL, INTENT(IN)  :: mean      ! Mean of the Gaussian distribution
        REAL, INTENT(IN)  :: stddev    ! Standard deviation of the Gaussian distribution
        REAL, INTENT(IN)  :: maxdev    ! Maximum allowed multiple of standard deviation
        REAL, INTENT(OUT) :: result    ! Output sampled value
        REAL :: u1, u2, z
        REAL :: pi
        pi = 3.141592653589793

        ! Generate two uniform random numbers in (0,1)
        CALL ran2(idum, xrand)
        u1 = xrand
        CALL ran2(idum, xrand)
        u2 = xrand

        ! Apply Box-Muller transform
        z = SQRT(-2.0 * LOG(u1)) * COS(2.0 * pi * u2)
        
        if (ABS(z) > maxdev) then
            if (z>0) then
                z = maxdev
            else
                z = -maxdev
            endif
        endif

        ! Scale and shift to match the desired mean and standard deviation
        result = mean + stddev * z
        
    END SUBROUTINE sample_gaussian
    
end module marche
    