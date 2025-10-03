module const
  implicit none

  integer,       parameter:: xp      = selected_real_kind(15)
  real(kind=xp), parameter:: unit_l  = 0.227550518265197E+28_xp
  real(kind=xp), parameter:: G       = 6.673E-11_xp
  real(kind=xp), parameter:: c       = 2.998E+8_xp
  real(kind=xp), parameter:: mp      = 1.673E-27_xp
  real(kind=xp), parameter:: me      = 9.109E-31_xp
  real(kind=xp), parameter:: kb      = 1.38E-23_xp
  real(kind=xp), parameter:: pi      = atan(1.E0_xp)*4.E0_xp
  real(kind=xp), parameter:: H_0     = 67.4_xp
  real(kind=xp), parameter:: mu      = 0.6
  real(kind=xp), parameter:: m_sun   = 1.99E+30
  real(kind=xp), parameter:: omega_m = 0.307114988565445E+00
  real(kind=xp), parameter:: omega_b = 0.450000017881393E-01
  real(kind=xp), parameter:: omega_c = omega_m-omega_b
  real(kind=xp), parameter:: rho_c   = 8.626795609585417E-27
  real(kind=xp), parameter:: pc      = 3.086E+16
end module const

module func
  use const
  implicit none

contains

  function gal_r500() result(i)

    real(kind=xp), dimension(1887) :: xcen,ycen,zcen
    !real(kind=xp), dimension(:), allocatable ::xdm,ydm,zdm,mdm,r
    real(kind=xp), dimension(1510846) ::xdm,ydm,zdm,mdm,rdm
    !real(kind=xp), dimension(364614)  ::x,y,z,r,m
    real(kind=xp), dimension(6779748)  ::x,y,z,r,m
    real(kind=xp), dimension(3) :: gal !31
    real(kind=xp), dimension(400):: m_rad,n_rad
    real(kind=xp), dimension(7)::  b 
    real(kind=xp), dimension(11):: a
    real(kind=xp) :: mtot,l,den,norm,ntot,r500,norm500,n500,m500,r200,norm200,n200,m200
    integer :: i,j,test,n,k

    test=1

 


    if(test==1) then
       !n=346614
       n=6779748
       !open(1,file='rdr_00251_l17.hydro')
       open(1,file='rdr_00251_l17.hydro')
       
       i=0
       do
          read(1,*,end=1) a
          i=i+1
          !write(*,*) i
          !write(*,*) "coucou"
          !r(i)=sqrt((a(1)-cen_x)**2+(a(2)-cen_y)**2+(a(3)-cen_z)**2)
          !v(i)=sqrt(a(4)**2+a(5)**2+a(6)**2)
          !P(i)=a(11)
          !rho(i)=a(7)
          !T(i)=a(10)
          m(i)=a(9)
          x(i)=a(1)
          y(i)=a(2)
          z(i)=a(3)
          !vx(i)=a(4)
          !vy(i)=a(5)
          !vz(i)=a(6)
     
       enddo
     
1      write(*,*) i
       close(1)

    endif

    if(test==2) then
       
       n=1510846
       
       open(2,file='rdr_00251_part.part')
  

       i=0
       do
          read(1,*,end=1) b
          i=i+1
          !write(*,*) i
          xdm(i)=b(1)
          ydm(i)=b(2)
          zdm(i)=b(3)
          !vx(i)=b(4)
          !vy(i)=b(5)
          !vz(i)=b(6)
          mdm(i)=b(7)
          
       enddo
        
2      write(*,*) mdm(1)
  
       close(2)
       
    endif
  

    !open(3,file='list_gal_251.dat_js_nocontam')
    open(3,file='gal_in_5Mpc.txt')
    
    !read(1,*) l
    
    !l=int(l)

    !allocate(xdm(l))
    !allocate(ydm(l))
    !allocate(zdm(l))
    !allocate(mdm(l))
    !allocate(r(l))
  
    !read(1,*) xdm
    !read(1,*) ydm
    !read(1,*) zdm
    !read(1,*)
    !read(1,*)
    !read(1,*)
    !read(1,*) mdm
    i=0
    !n=346614
    do
       i=i+1
       read(3,*,end=3) xcen(i),ycen(i),zcen(i)
       
       
       !if (i==2) then
       write(*,*) xcen(i),ycen(i),zcen(i)!gal(1),gal(2),gal(3)
       !endif
       
       
       !xcen(i)=(gal(1)-0.5)*(unit_l/3.08567758128E21) !4 5 6
       !ycen(i)=(gal(2)-0.5)*(unit_l/3.08567758128E21)
       !zcen(i)=(gal(3)-0.5)*(unit_l/3.08567758128E21)

        
     
   

    enddo

 
3   write(*,*) "fin lecture fichiers"

    close(3)

    write(*,*) "test"
    
    write(*,*) 
    
    close(2)

    !open(4,file='r500-r200_gal.txt')

  
    do k=1,2
       write(*,*) k
       m_rad=0
       n_rad=0
       norm=0
       norm500=0
       m500=0
       n500=0
       r500=0
       norm200=0
       m200=0
       n200=0
       r200=0
     
       do i=1,n
          
          r(i)=sqrt((xcen(k)-x(i))**2+(ycen(k)-y(i))**2+(zcen(k)-z(i))**2)
     
       enddo

  
       !write(*,*) "r fini"
       !write(*,*) r(1)
       
       m_rad=0
       do i=1,400
          do j=1,n
             if(r(j)>(i*10) .and. r(j)<((i+1)*10)) then
                m_rad(i)=m_rad(i)+m(j)
                n_rad(i)=n_rad(i)+1
             endif
          enddo
       enddo
       
       !write(*,*) m_rad
       
       !write(*,*) "bin fini"
       !write(*,*) m_rad(1)
       
       do i=1,400
          mtot=sum(m_rad(1:i))
          ntot=sum(n_rad(1:i))
          den=(3*m_sun*mtot)/(4*pi*(i*10*pc*1E3)**3)
          norm=den/(rho_c*(omega_b/omega_m))
          
          if (norm>500) then
             r500=i*10
             m500=mtot
             n500=ntot
             norm500=norm
          endif

          if (norm>200) then
             r200=i*10
             m200=mtot
             n200=ntot
             norm200=norm
          endif

        

          !write(*,*) "gal numero",k, "norm",norm,"mass",mtot,"r",i*10,"n",ntot
         

       enddo
     
           
       write(*,*) "gal numero",k, "norm",norm500,"m500",m500,"r500",r500,"n",n500
       write(*,*) "gal numero",k, "norm",norm200,"m200",m200,"r200",r200,"n",n200
       !write(4,*) k,r500,n500,m500,r200,n200,m200
    enddo
    
    !close(4)

  end function gal_r500

  function r200_subhalos() result(i)

    real(kind=xp), dimension(31)::gal
    real(kind=xp), dimension(1510846)::x,y,z,m,r,vx,vy,vz
    integer::i,j,g,n,bin
    real(kind=xp):: xcen,ycen,zcen,norm,den,mtot,ntot,r200
    real(kind=xp), dimension(400)::m_rad,n_rad
    real(kind=xp), dimension(5000)::r200_gals

    !open(2,file='virgo_xyz_dm.dat',access='direct',form='unformatted',recl=15)

    open(2,file="rdr_00251_part.part")
    i=0
    do
       i=i+1
       read(2,*,end=2) x(i),y(i),z(i),vx(i),vy(i),vz(i),m(i)
    !read(2,*) y
    !read(2,*) z
    !read(2,*)
    !read(2,*)
    !read(2,*)
       !read(2,*) m
    enddo
    

2   close(2)

    write(*,*) "lecture rdr fini"

    open(1,file='list_gal_251.dat_js_nocontam')

    n=1510846
    g=0
    
    do 
       read(1,*,end=1) gal

       xcen=gal(21)
       ycen=gal(22)
       zcen=gal(23)

       

       !write(*,*) xcen

        !cond=np.logical_and(xcen>0.47783,np.logical_and(xcen<0.49139,np.logical_and(ycen>0.50131,np.logical_and(ycen<0.51487,np.logical_and(zcen>0.49009,zcen<0.50365)))))
        if (xcen>0.47783 .and. xcen<0.49139 .and. ycen>0.50131 .and. ycen<0.51487 .and. zcen>0.49009 .and. zcen<0.50365) then

           xcen=(gal(21)-0.5)*(unit_l/3.08567758128E21)
           ycen=(gal(22)-0.5)*(unit_l/3.08567758128E21)
           zcen=(gal(23)-0.5)*(unit_l/3.08567758128E21)

           r200=0
           m_rad=0
           n_rad=0
           r=0

           g=g+1
           write(*,*) g
           do i=1,n
              r(i)=sqrt((xcen-x(i))**2+(ycen-y(i))**2+(zcen-z(i))**2)
              bin=int(r(i)/10)
              !write(*,*) r(i),bin
              if(bin<401) then
                 m_rad(bin)=m_rad(bin)+m(i)
                 n_rad(bin)=n_rad(bin)+1
              endif
           enddo

           write(*,*) "sum n_rad",sum(n_rad(:))

       

           !m_rad=0
         
           !do j=1,n
             
           !   if(r(j)>(i*10) .and. r(j)<((i+1)*10)) then
           !      m_rad(i)=m_rad(i)+m(j)
           !      n_rad(i)=n_rad(i)+1
           !   endif
       
       !write(*,*) m_rad

           
       !write(*,*) "bin fini"
       !write(*,*) m_rad(1)
       
           do i=1,400
              mtot=sum(m_rad(1:i))
              ntot=sum(n_rad(1:i))
              den=(3*m_sun*mtot)/(4*pi*(i*10*pc*1E3)**3)
              norm=den/(rho_c*(omega_c/omega_m))
              
              if (norm>200) then
                 r200=i*10
                 !m200=mtot
                 !n200=ntot
                 !norm200=norm
              endif
           enddo
           
        write(*,*) "r200=",r200   
        r200_gals(g)=r200
         
        endif
     enddo

1     close(1)

     open(3,file='r200_gals.txt')
     write(3,*) r200_gals(1:g)

     close(3)
    

   end function r200_subhalos
