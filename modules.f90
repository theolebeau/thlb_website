module const
  implicit none

  integer,       parameter:: xp      = selected_real_kind(15)
  real(kind=xp), parameter:: unit_l  = 0.227550518265197E+28_xp
  real(kind=xp), parameter:: G       = 6.67430E-11_xp
  real(kind=xp), parameter:: c       = 2.998E+8_xp
  real(kind=xp), parameter:: mp      = 1.673E-27_xp
  real(kind=xp), parameter:: me      = 9.109E-31_xp
  real(kind=xp), parameter:: kb      = 1.38E-23_xp
  real(kind=xp), parameter:: pi      = atan(1.E0_xp)*4.E0_xp
  real(kind=xp), parameter:: H_0     = 67.4_xp
  real(kind=xp), parameter:: mu      = 0.6_xp
  real(kind=xp), parameter:: m_sun   = 1.99E+30_xp
  real(kind=xp), parameter:: omega_m = 0.307114988565445E+00_xp
  real(kind=xp), parameter:: omega_b = 0.450000017881393E-01_xp
  real(kind=xp), parameter:: omega_c = omega_m-omega_b
  real(kind=xp), parameter:: rho_c   = 8.626795609585417E-27_xp
  real(kind=xp), parameter:: pc      = 3.086E+16_xp
  real(kind=xp), parameter:: sigma_t = 6.6524587158E-29_xp
end module const

module func
  use const
  use mpi
  implicit none

contains

  function ascii_to_bin() result(i)

    real(kind=xp), dimension(11)      :: a
    real(kind=xp), dimension(7)       :: b 
    real(kind=xp), dimension(400000000)  :: rho,P,T,r,vx,vy,vz,m,x,y,z,mdm,lvl
    integer                           :: i,test
    real(kind=xp)                     :: cen_x,cen_y,cen_z

  
    !x=0.48461053 !M87
    !y=0.50809777
    !z=0.49687093
    
    !x=0.47215262 !gal1
    !y=0.51401520
    !z=0.49303299
  
    !x=0.48461068 !Virgo
    !y=0.50809848
    !z=0.49687076

    !x=0.48287964 !gal10
    !y=0.50985229
    !z=0.49664778

  

    !cen_x= (x-0.5)*(unit_l/3.08567758128E21)
    !cen_y= (y-0.5)*(unit_l/3.08567758128E21)
    !cen_z= (z-0.5)*(unit_l/3.08567758128E21)

    test=1
    
    if(test==1) then
       open(1,file='/data/cluster/tlebeau/virgo/rdr_00251_l16_high.hydro')
       open(2,file='virgo_xyz_hydro_l16_high.dat',form='unformatted')
       !open(1,file='virgo_xyz_hydro_l21.dat',form='unformatted')
       !open(2,file='virgo_xyz_hydro_l21_bin.dat',form='unformatted')
     
       i=0
       do
          i=i+1
          !read(1,*,end=1) rho(i),T(i),P(i),x(i),y(i),z(i),vx(i),vy(i),vz(i),m(i),lvl(i)
          read(1,*,end=1) x(i),y(i),z(i),vx(i),vy(i),vz(i),rho(i),lvl(i),m(i),T(i),P(i)
  
          if(modulo(i,1000000)==0) then
             write(*,*) "i",i
          endif
          !write(*,*) i
          !write(*,*) "coucou"
          !r(i)=sqrt((a(1)-cen_x)**2+(a(2)-cen_y)**2+(a(3)-cen_z)**2)
          !v(i)=sqrt(a(4)**2+a(5)**2+a(6)**2)
          !P(i)=a(11)
          !rho(i)=a(7)
          !T(i)=a(10)
          !m(i)=a(9)
          !x(i)=a(1)
          !y(i)=a(2)
          !z(i)=a(3)
          !vx(i)=a(4)
          !vy(i)=a(5)
          !vz(i)=a(6)
     
       enddo
     
1      write(*,*) i
       write(*,*) x(i-1)
       write(*,*) x(i)
       
       write(2) i-1
       write(2) rho(1:i)
       write(2) T(1:i)
       write(2) P(1:i)
       write(2) x(1:i)
       write(2) y(1:i)
       write(2) z(1:i)
       write(2) vx(1:i)
       write(2) vy(1:i)
       write(2) vz(1:i)
       write(2) m(1:i)
       write(2) lvl(1:i)
     
       close(1)
       close(2)
       write(*,*) "conversion done"

    endif
    
    if(test==2) then
       open(1,file='/data/cluster/tlebeau/virgo/rdr_00251_part_496.part')
       open(2,file='/data/cluster/tlebeau/virgo/virgo_xyz_dm_low_res.dat',form='unformatted')
       
       i=0
       do
          
          i=i+1
          read(1,*,end=2) x(i),y(i),z(i),vx(i),vy(i),vz(i),mdm(i)
          !write(*,*) i
          if(modulo(i,1000000)==0) then
             write(*,*) i
          endif
          !write(*,*) i
          !r(i)=sqrt((a(1)-cen_x)**2+(a(2)-cen_y)**2+(a(3)-cen_z)**2)
          !v(i)=sqrt(a(4)**2+a(5)**2+a(6)**2)
          !x(i)=b(1)
          !y(i)=b(2)
          !z(i)=b(3)
          !vx(i)=b(4)
          !vy(i)=b(5)
          !vz(i)=b(6)
          !mdm(i)=b(7)
     
       enddo
       
2      write(*,*) "npart", i

       stop

       !write(*,*) maxval(x(1:i)),minval(x(1:i)),maxval(y(1:i)),minval(y(1:i)),maxval(z(1:i)),minval(z(1:i))

       !stop "test"
     
       write(2) i-1
       write(2) x(1:i)
       write(2) y(1:i)
       write(2) z(1:i)
       write(2) vx(1:i)
       write(2) vy(1:i)
       write(2) vz(1:i)
       write(2) mdm(1:i)
   
       close(1)
       close(2)

    endif

  end function ascii_to_bin
  

  function create_map(proj,width,unit,n,weighting,filein,fileout) result(cen_x)

    real(kind=xp), dimension(11)::a
    real(kind=xp), dimension(:,:),allocatable :: map,weights
    real(kind=xp) :: cen_x,cen_y,cen_z,minx,maxx,miny,maxy,minz,maxz,lenx,leny,lenz,dx,dy,dz,dist,volume,ratio,start,finish
    
    real(kind=xp), dimension(:), allocatable :: x,y,z,P,n_e,T,m,vx,vy,vz,lvl
    integer :: i,j,k,l,distx,disty,xinf,xsup,yinf,ysup,ix,iy,test,width,n,nx,ny,nz,ncell
    character(len=*)::filein,fileout
    character(len=2)::weighting
    character(len=1)::proj,unit

    write(*,*) "projection: ",proj
    write(*,*)  "unit: ",unit

    test=2
    if (test==1) then
       open(1,file='hydro_remove_rvir.dat',form='unformatted')
       !open(1,file='rdr_00251_l17.hydro')
       i=0
       do
          !read(1,end=1) a
          read(1,end=1) n_e(i),T(i),P(i),x(i),y(i),z(i),vx(i),vy(i),vz(i),m(i),lvl(i)
          i=i+1
          !P(i)=a(11)/10
          !n_e(i)=a(7)*(1e6/0.864)
          !write(*,*) "ne",n_e(1:100)
          !stop
       
          !T(i)=a(10)
          !m(i)=a(9)
          !x(i)=a(1)
          !y(i)=a(2)
          !z(i)=a(3)
          !vx(i)=a(4)
          !vy(i)=a(5)
          !vz(i)=a(6)
          !lvl(i)=a(8)

       enddo
       
  
1      write(*,*) "fin lecture fichier"

       P(:)=P(:)/10
       n_e(:)=n_e(:)*(1e6/0.864)

       close(1)

    endif

    if (test==2) then
       !open(1,file='virgo_xyz_hydro_l21_gal_clean.dat',form='unformatted')
       open(1,file=filein,form='unformatted')
       !open(1,file='virgo_xyz_hydro_l17.dat',form='unformatted')
       !ncell=6779748
       !ncell=5756892
       !ncell=166999251
       !ncell=223149095
       read(1) ncell
       write(*,*) "ncell",ncell
       allocate(n_e(ncell),T(ncell),P(ncell),x(ncell),y(ncell),z(ncell),vx(ncell),vy(ncell),vz(ncell),m(ncell),lvl(ncell))
       read(1) n_e
       n_e=n_e/0.864
       write(*,*) "n_e",n_e(1:10)
       write(*,*) minval(n_e(:)),maxval(n_e(:))
       read(1) T
       T=T*(kb/1.602e-16)
       read(1) P
       P=P/10
       P=P/(1.602e-10)
       P=P*(0.76/0.864)
       write(*,*) minval(P(:)),maxval(P(:))
       read(1) x
       read(1) y
       read(1) z
       read(1) vx
       read(1) vy
       read(1) vz
       read(1) m
       read(1) lvl

       write(*,*) "min lvl", minval(lvl(:))
       write(*,*) "vz", vz
       stop

       close(1)
       
    endif 
  
  cen_x=0.48461068    !Virgo
  cen_y=0.50809848
  cen_z=0.49687076

  !cen_x=0.48287964   !gal 10
  !cen_y=0.50985229
  !cen_z=0.49664778
    
  !cen_x=0.48268557    !gal 33 
  !cen_y=0.50639194
  !cen_z=0.49427480

  !cen_x=0.48512837    !gal 14
  !cen_y=0.50834167
  !cen_z=0.49735767 
  
  cen_x= (cen_x-0.5)*(unit_l/3.08567758128E21)
  cen_y= (cen_y-0.5)*(unit_l/3.08567758128E21)
  cen_z= (cen_z-0.5)*(unit_l/3.08567758128E21)

  !width=300
  minx=cen_x-width
  maxx=cen_x+width
  miny=cen_y-width
  maxy=cen_y+width
  minz=cen_z-width
  maxz=cen_z+width
  
  nx=n
  ny=n
  nz=n

  allocate(map(nx,ny))
  allocate(weights(nx,ny))

  map=0
  weights=0
  
  lenx=maxx-minx
  leny=maxy-miny
  lenz=maxz-minz
   
  dx=lenx/nx
  write(*,*) "dx",dx
  dy=leny/ny
  dz=lenz/nz

  write(*,*) "test"

  !open(2,file='map_2000px_15Mpc_zrot2_ne_los.txt')
  open(2,file=fileout,form='unformatted')
  
  !do i=1,nx
  !   write(*,*) i 
  !   do j=1,ny
        !write(*,*) i,j
  !      do k=1,6779748
  !         if (x(k)>(minx+dx*i) .and. x(k)<(minx+dx*(i+1)) .and. y(k)>(miny+dy*j) .and. y(k)<(miny+dy*(j+1))) then
  !            map(i,j)=map(i,j)+rho(k)
  !         endif
  !      enddo
  !   enddo
  !   write(2,*) map(i,:)
  !enddo

  write(*,*) "k=0"
  write(*,*) "ncell",ncell
  write(*,*) "len m",size(m)
  do k=1,size(m)
    call cpu_time(start)

     if (proj=="x") then   
        i=int(((y(k)-miny)/leny)*ny)
        j=int(((z(k)-minz)/lenz)*nz)
        l=int(((x(k)-minx)/lenx)*nx)
     endif

     if (proj=="y") then   
        i=int(((x(k)-minx)/lenx)*nx)
        j=int(((z(k)-minz)/lenz)*nz)
        l=int(((y(k)-miny)/leny)*ny)
     endif

     if (proj=="z") then  
        i=int(((x(k)-minx)/lenx)*nx)
        j=int(((y(k)-miny)/leny)*ny)
        l=int(((z(k)-minz)/lenz)*nz)
     endif

     
     !write(*,*) dist
     dist=737441/2**(lvl(k))
     volume=(dist*(pc*1e3))**3

     !write(*,*) "lvl", lvl(k)
     !write(*,*) "dist",dist
     
     distx=int((dist)/(dx*2))
     disty=int((dist)/(dy*2))

     !if (distx<0) then
     !   distx=0
     !   disty=0

     !write(*,*) "distx",distx

     !stop

     if(i>0 .and. j>0 .and. i<nx+1 .and. j<ny+1 .and. T(k)>8.6142e-3) then
     !if(i>0 .and. j>0 .and. i<nx+1 .and. j<ny+1) then
     !if(i>0 .and. j>0 .and. i<nx+1 .and. j<ny+1 .and. T(k)>1e5) then! .and. l>0 .and. l<nz+1) then !.and. T(k)>1e5) then
        !if(lvl(k)<12) then
        !   write(*,*) dist,distx,disty,lvl(k)
        !   write(*,*) i,j
        !endif
        if((i-distx)>0) then
           xinf=i-distx
        else
           xinf=1
        endif

        if((i+distx)<nx+1) then
           xsup=i+distx
        else
           xsup=nx
        endif

         if((j-disty)>0) then
           yinf=j-disty
        else
           yinf=1
        endif

        if((j+disty)<ny+1) then
           ysup=j+disty
        else
           ysup=ny
        endif

        !write(*,*) "bornes",xinf,xsup,yinf,ysup
        do ix=xinf,xsup
           do iy=yinf,ysup
              if(unit=="T") then
                 if (weighting=="mw") then
                    map(ix,iy)=map(ix,iy)+T(k)*m(k) !n_e(k) !P(k)*m(k)
                    weights(ix,iy)=weights(ix,iy)+m(k)
                 else if (weighting=="sl") then
                    !map(ix,iy)=map(ix,iy)+T(k)*m(k)*n_e(k)*T(k)**(-0.75)
                    !weights(ix,iy)=weights(ix,iy)+m(k)*n_e(k)*T(k)**(-0.75)
                    map(ix,iy)=map(ix,iy)+T(k)*n_e(k)*n_e(k)*T(k)**(-0.75)
                    weights(ix,iy)=weights(ix,iy)+n_e(k)*n_e(k)*T(k)**(-0.75)

                 else if (weighting=="ew") then
                    !map(ix,iy)=map(ix,iy)+T(k)*m(k)*n_e(k)*T(k)**(0.5)
                    !weights(ix,iy)=weights(ix,iy)+m(k)*n_e(k)*T(k)**(0.5)
                    map(ix,iy)=map(ix,iy)+T(k)*n_e(k)*n_e(k)*T(k)**(0.5)
                    weights(ix,iy)=weights(ix,iy)+n_e(k)*n_e(k)*T(k)**(0.5)
                 endif
                    
                 !write(*,*) "T"
              else if(unit=="n") then
                 map(ix,iy)=map(ix,iy)+n_e(k)*dist*pc*1e3
              else if(unit=="P") then
                 map(ix,iy)=map(ix,iy)+P(k)*m(k)
                 weights(ix,iy)=weights(ix,iy)+m(k)
              else if(unit=="y") then
                  map(ix,iy)=map(ix,iy)+P(k)*dist*pc*1e3
              else if(unit=="v") then
                 map(ix,iy)=map(ix,iy)+vz(k)
                 weights(ix,iy)=weights(ix,iy)+1

              endif
                 
           enddo
        enddo
        
     endif

     ratio=real(k)/6432912.0
     write(*,*) "ratio", ratio, "k", k
        call cpu_time(finish)
    print '("Time =")',finish-start
    stop
  
  enddo

  if(unit=="P" .or. unit=="T" .or. unit=="v") then
     map(:,:)=map(:,:)/weights(:,:)
  endif

  if(unit=="y") then
      map(:,:)=map(:,:)*(sigma_t/(me*c**2))
  end if

  write(*,*) "rangement fini"

  !do i=1,nx
  write(2) map
  !enddo
  
  !write(*,*) "map", map(:,:)
  
  !write(*,*) map 
  
  !write(2,*) nx,ny
  
  close(2)

  deallocate(n_e,T,P,x,y,z,vx,vy,vz,m,lvl,map)


  end function create_map


  function create_map_dm() result(cen_x)

    real(kind=xp), dimension(7)::a
    real(kind=xp), dimension(:,:),allocatable :: map
    real(kind=xp) :: cen_x,cen_y,cen_z,minx,maxx,miny,maxy,minz,maxz,nx,ny,nz,lenx,leny,lenz,dx,dy,dz,dist,width
    
    real(kind=xp), dimension(10000000)  :: x,y,z,m
    integer :: i,j,k,l,distx,disty,xinf,xsup,yinf,ysup,ix,iy,test,n

    test=1
    if (test==1) then
       open(1,file='rdr_00251_part.part')

       
       i=0
       do
          read(1,*,end=1) a
          i=i+1
          m(i)=a(7)
          x(i)=a(1)
          y(i)=a(2)
          z(i)=a(3)
          if (i==1) then
             maxx=x(i)
             minx=x(i)
             maxy=y(i)
             miny=y(i)
             minz=z(i)
             maxz=z(i)
          else
             if(x(i)>maxx) then
                maxx=x(i)
             endif
             
             if(x(i)<minx) then
                minx=x(i)
             endif

             if(y(i)>maxy) then
                maxy=y(i)
             endif

             if(y(i)<miny) then
                miny=y(i)
             endif

             if(z(i)>maxz) then
                maxz=z(i)
             endif

             if(z(i)<minz) then
                minz=z(i)
             endif
          endif
          
             

          
       

       enddo
  
1      write(*,*) "fin lecture fichier"
       write(*,*) "minx :",minx,"maxx",maxx
       write(*,*) "miny :",miny,"maxy",maxy
       write(*,*) "minz :",minz,"maxz",maxz
       

       close(1)

    endif
    n=i
    cen_x=0.48461068
    cen_y=0.50809848
    cen_z=0.49687076
  
    cen_x= (cen_x-0.5)*(unit_l/3.08567758128E21)
    cen_y= (cen_y-0.5)*(unit_l/3.08567758128E21)
    cen_z= (cen_z-0.5)*(unit_l/3.08567758128E21)

    write(*,*) cen_x,cen_y,cen_z

    width=5000
    minx=cen_x-width
    maxx=cen_x+width
    miny=cen_y-width
    maxy=cen_y+width
    minz=cen_z-width
    maxz=cen_z+width
  
    nx=2000
    ny=2000
    nz=2000

    allocate(map(nx,ny))
    
    lenx=maxx-minx
    leny=maxy-miny
    lenz=maxz-minz
    
    dx=lenx/nx
    dy=leny/ny
    dz=lenz/nz
    
    write(*,*) "test"
    
    open(2,file='map_2000px_dm_z_profond.txt')
  
  !do i=1,nx
  !   write(*,*) i 
  !   do j=1,ny
        !write(*,*) i,j
  !      do k=1,6779748
  !         if (x(k)>(minx+dx*i) .and. x(k)<(minx+dx*(i+1)) .and. y(k)>(miny+dy*j) .and. y(k)<(miny+dy*(j+1))) then
  !            map(i,j)=map(i,j)+rho(k)
  !         endif
  !      enddo
  !   enddo
  !   write(2,*) map(i,:)
  !enddo

    do k=1,n
       i=int(((x(k)-minx)/lenx)*nx)
       j=int(((y(k)-miny)/leny)*ny)
       l=int(((z(k)-minz)/lenz)*nz)
       !dist=737441/2**(lvl(k))
       !write(*,*) dist
     
       !distx=int((dist)/(dx*2))
       !disty=int((dist)/(dy*2))
     
       if(i>0 .and. j>0 .and. i<nx+1 .and. j<ny+1 .and. l>0 .and. l<nz+1 ) then !.and. T(k)>1e7) then
        !if(lvl(k)<12) then
           !write(*,*) dist,distx,disty,lvl(k)
           !write(*,*) i,j
        !endif
        !if((i-distx)>0) then
        !   xinf=i-distx
        !else
        !   xinf=1
        !endif

        !if((i+distx)<nx+1) then
        !   xsup=i+distx
        !else
        !   xsup=nx
        !endif

        !if((j-disty)>0) then
        !  yinf=j-disty
        !else
        !   yinf=1
        !endif

        !if((j+disty)<ny+1) then
        !   ysup=j+disty
        !else
        !   ysup=ny
        !endif

        !write(*,*) "bornes",xinf,xsup,yinf,ysup
        !do ix=xinf,xsup
        !   do iy=yinf,ysup
          map(i,j)=map(i,j)+m(k)
        !   enddo
        !enddo
        
        !map(i,j)=map(i,j)+rho(k)
       endif
  
    enddo

    write(*,*) "rangement fini"
    
    do i=1,nx
       write(2,*) map(:,i)
    enddo
  
  
  
  !write(*,*) map 
  
  !write(2,*) nx,ny
  
    close(2)

    deallocate(map)


  end function create_map_dm

  function remove_gal() result(i)
     real(kind=xp), dimension(6000) :: xcen,ycen,zcen,rvir,mgal
     !real(kind=xp), dimension(:), allocatable ::xdm,ydm,zdm,mdm,r
     real(kind=xp), dimension(60000006) ::xdm,ydm,zdm,mdm,rdm
     !real(kind=xp), dimension(364614)  ::x,y,z,r,m
     !real(kind=xp), dimension(223149095)  ::x,y,z,r,m,ne,T,P,vx,vy,vz,lvl
     real(kind=xp), dimension(:),allocatable  :: x,y,z,r,m,ne,T,P,vx,vy,vz,lvl
     real(kind=xp), dimension(31) :: gal !31
     real(kind=xp), dimension(400):: m_rad,n_rad
     real(kind=xp), dimension(7)::  b 
     real(kind=xp), dimension(11):: a
     real(kind=xp), dimension(:), allocatable :: liste,liste2
     real(kind=xp) :: mtot,den,norm,ntot,r500,norm500,n500,m500,r200,norm200,n200,m200,cen_x,cen_y,cen_z
     integer :: i,j,test,ncell,k,index,l,ngal,n
  
     test=1

 


     open(1,file='/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l15_high.dat',form='unformatted')
     !open(1,file='/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l19_low.dat',form='unformatted')


     read(1) ncell
     write(*,*) ncell
     allocate(ne(ncell),T(ncell),P(ncell),x(ncell),y(ncell),z(ncell),vx(ncell),vy(ncell),vz(ncell),m(ncell),lvl(ncell),r(ncell),liste(ncell))
     !allocate(ne(ncell),T(ncell),P(ncell),x(ncell),y(ncell),z(ncell),vx(ncell),vy(ncell),vz(ncell),m(ncell),lvl(ncell),r(ncell),liste(ncell))

     read(1) ne
     read(1) T
     read(1) P
     read(1) x
     read(1) y
     read(1) z
     read(1) vx
     read(1) vy
     read(1) vz
     read(1) m
     read(1) lvl
        
     close(1)

     write(*,*) "fin lecture hydro"

     !write(*,*) "blabla"
  

     open(3,file='/data/cluster/tlebeau/virgo/list_gal_251.dat_js_nocontam_high_res')
     !open(3,file='/data/cluster/tlebeau/virgo/list_gal_251.dat_js_nocontam_low_res')
     
     cen_x=0.48461068
     cen_y=0.50809848
     cen_z=0.49687076
      
     cen_x= (cen_x-0.5)*(unit_l/3.08567758128E21)
     cen_y= (cen_y-0.5)*(unit_l/3.08567758128E21)
     cen_z= (cen_z-0.5)*(unit_l/3.08567758128E21)
      
     i=0
     do
        
        read(3,*,end=3) gal
        !if(gal(4)>0.47783 .and. gal(4)<0.49139 .and. gal(5)>0.50131 .and. gal(5)<0.51487 .and. gal(6)>0.49009 .and. gal(6)<0.50365) then
        if(gal(4)>0.46 .and. gal(4)<0.51 .and. gal(5)>0.48 .and. gal(5)<0.53 .and. gal(6)>0.47 .and. gal(6)<0.52) then
           i=i+1
           xcen(i)=(gal(4)-0.5)*(unit_l/3.08567758128E21)
           ycen(i)=(gal(5)-0.5)*(unit_l/3.08567758128E21)
           zcen(i)=(gal(6)-0.5)*(unit_l/3.08567758128E21)
           rvir(i)=gal(24)*(unit_l/3.08567758128E21)
           mgal(i)=gal(3)
           !write(*,*) "i",i,"x",xcen(i),"y",ycen(i),"z",zcen(i),"rvir",rvir(i),"m",mgal(i)



           !test to understand peak on 2D proj EM rad prof (splashback project), gal at x=-0.4179 kpc and y = -1.0230 kpc ??
           !if (xcen(i)<0 .and. xcen(i)>-1000 .and. ycen(i)<0 .and. ycen(i)>-2000) then
           !   write(*,*) "gal",i,"x",xcen(i),"y",ycen(i),"z",zcen(i)!,"rvir",gal(24),"m",gal(3)

           !end if
        endif 

     enddo



 
3    write(*,*) "fin lecture fichiers"

     !write(*,*) "end"

     !stop

     
     close(3)

     write(*,*) "nbr gal", i

     !stop

     ngal=i
      
     write(*,*) "test"
      
     !open(4,file='virgo_xyz_hydro_l19_gal_clean_m1e8.5.dat',form='unformatted')
     open(4,file='./virgo_xyz_files/virgo_xyz_hydro_l15_high_gal_clean_m1e8.5.dat',form='unformatted')

     write(*,*) "taille x", size(x)
     write(*,*) "test avant if"
      
     do i=1,size(x)
        liste(i)=i
     enddo
     
     index=size(x)
     
     do k=3,ngal
       
        if (mgal(k)>10**(8.5)) then

           allocate(liste2(index))

           index=0
            
           do i=1,size(liste)
           
              r(liste(i))=sqrt((xcen(k)-x(liste(i)))**2+(ycen(k)-y(liste(i)))**2+(zcen(k)-z(liste(i)))**2)
              if (r(liste(i))>rvir(k)) then
                 index=index+1
                 liste2(index)=liste(i)
              endif
           enddo

           deallocate(liste)
           allocate(liste(index)) 
           liste(:)=liste2(1:index)
           deallocate(liste2)
        endif
                                     
        write(*,*) "gal nbr",k," len liste", size(liste)
                                     
     enddo
     
     write(*,*) "len liste", size(liste)
     
     write(4) size(liste)
     write(4) ne(liste)
     write(4) T(liste)
     write(4) P(liste)
     write(4) x(liste)
     write(4) y(liste)
     write(4) z(liste)
     write(4) vx(liste)
     write(4) vy(liste)
     write(4) vz(liste)
     write(4) m(liste)
     write(4) lvl(liste)
      
     close(4)
     
     deallocate(ne,T,P,x,y,z,vx,vy,vz,m,lvl,liste)
  
   end function remove_gal

    function rotation() result(phi)
      real(kind=xp), dimension(:),allocatable:: ne,T,P,m,lvl
      real(kind=xp), dimension(:,:),allocatable::pos,v
      real(kind=xp) :: mw_x,mw_y,mw_z,phi,theta,s,c,zero,un,fil_x,fil_y,fil_z,vx_virgo_halo,vy_virgo_halo,vz_virgo_halo
      real(kind=xp), dimension(3,3):: rx,ry,rz
      integer:: ncell
      real(kind=xp), dimension(3):: mw,fil,rotaxis,virgo

      
      write(*,*) "open datafile"
      !open(1,file='virgo_xyz_hydro_l19_gal_clean_m1e8.5.dat',form='unformatted')
      open(1,file='./virgo_xyz_files/virgo_xyz_hydro_l15_high_gal_clean_m1e8.5.dat',form='unformatted')
      
      read(1) ncell
      allocate(ne(ncell),T(ncell),P(ncell),m(ncell),lvl(ncell),pos(3,ncell),v(3,ncell))
      
      read(1) ne
      !ne=ne*(1E6/0.864)
      read(1) T
      read(1) P
      !P=P/10
      read(1) pos(1,:)
      read(1) pos(2,:)
      read(1) pos(3,:)
      read(1) v(1,:)
      read(1) v(2,:)
      read(1) v(3,:)
      read(1) m
      read(1) lvl

      close(1)


      !stop



      
      virgo=(/0.48461068,0.50809848,0.49687076/)
      fil=(/0.497,0.4984,0.4996/)
      mw=(/0.5,0.5,0.5/)

      !pos(:,1)=(/0,0,0/)
      pos(:,1)=(fil(:)-mw(:))*(unit_l/3.08567758128E21)

      write(*,*) "pos 1",pos(:,1)

      rotaxis=mw

      vx_virgo_halo = -507.8579
      vy_virgo_halo = 229.4530
      vz_virgo_halo = -136.9451

      v(1,:) = v(1,:) - vx_virgo_halo
      v(2,:) = v(2,:) - vy_virgo_halo
      v(3,:) = v(3,:) - vz_virgo_halo
      

      pos(1,:)=pos(1,:)-(virgo(1)-mw(1))*(unit_l/3.08567758128E21)
      pos(2,:)=pos(2,:)-(virgo(2)-mw(2))*(unit_l/3.08567758128E21)
      pos(3,:)=pos(3,:)-(virgo(3)-mw(3))*(unit_l/3.08567758128E21)

      

      !pos(:,1)=(/0,0,0/)


      !mw_x=0.5-0.48461068
      !mw_y=0.5-0.50809848
      !mw_z=0.5-0.49687076

      !fil_x=0.5-0.497
      !fil_y=0.5-0.4984
      !fil_z=0.5-0.4996

      
      
      zero=0.0
      un=1.0
      
      !phi=atan(mw_x/mw_z)
      !phi=atan(mw_y/mw_z)
      !theta=atan(mw_x/mw_y)
      !theta=atan(fil_x/fil_y)
      
      theta=atan((rotaxis(1)-virgo(1))/(rotaxis(2)-virgo(2)))

      !write(*,*) "phi",phi
      write(*,*) "theta",theta

      !theta=pi/2
      !phi=0

      write(*,*) "mw",pos(:,1)

      c=cos(theta)
      s=sin(theta)

      write(*,*) "cos",c,"sin",s

      rz(:,1)=(/c,-s,zero/)
      rz(:,2)=(/s,c,zero/)
      rz(:,3)=(/zero,zero,un/)
      write(*,*) "rz",rz

      pos=transpose(matmul(transpose(pos),rz))
      v=transpose(matmul(transpose(v),rz))

      write(*,*) "mw",pos(:,1)

      phi=atan(pos(2,1)/pos(3,1))
      !phi=

      write(*,*) "phi",phi

      c=cos(phi)
      s=sin(phi)

      

      !ry(:,1)=(/c,zero,s/)
      !ry(:,2)=(/zero,un,zero/)
      !ry(:,3)=(/-s,zero,c/)

      rx(:,1)=(/un,zero,zero/)
      rx(:,2)=(/zero,c,-s/)
      rx(:,3)=(/zero,s,c/)

      pos=transpose(matmul(transpose(pos),rx))
      v=transpose(matmul(transpose(v),rx))

      phi=pi/2

      write(*,*) "phi",phi

      c=cos(phi)
      s=sin(phi)

      write(*,*) "cos",c,"sin",s

      

      ry(:,1)=(/c,zero,s/)
      ry(:,2)=(/zero,un,zero/)
      ry(:,3)=(/-s,zero,c/)

      rx(:,1)=(/un,zero,zero/)
      rx(:,2)=(/zero,c,-s/)
      rx(:,3)=(/zero,s,c/)

      !pos=transpose(matmul(transpose(pos),ry))

      write(*,*) "pos 1",pos(:,1)

      pos(1,:)=pos(1,:)+(virgo(1)-mw(1))*(unit_l/3.08567758128E21)
      pos(2,:)=pos(2,:)+(virgo(2)-mw(2))*(unit_l/3.08567758128E21)
      pos(3,:)=pos(3,:)+(virgo(3)-mw(3))*(unit_l/3.08567758128E21)

      write(*,*) "pos 1",pos(:,1)

      !stop

      write(*,*) "ecriture"
      
      open(2,file='./virgo_xyz_files/virgo_xyz_hydro_l15_gal_clean_m1e8.5_MW_los.dat',form='unformatted')
      write(2) ncell
      write(2) ne
      write(2) T
      write(2) P
      write(2) pos(1,:)
      write(2) pos(2,:)
      write(2) pos(3,:)
      write(2) v(1,:)
      write(2) v(2,:)
      write(2) v(3,:)
      write(2) m
      write(2) lvl
      close(2)

      write(*,*) "data saved"

    end function rotation

    function rotation_DM() result(phi)
      real(kind=xp), dimension(:),allocatable::m
      real(kind=xp), dimension(:,:),allocatable::pos,v
      real(kind=xp) :: mw_x,mw_y,mw_z,phi,theta,s,c,zero,un,fil_x,fil_y,fil_z,vx_virgo_halo,vy_virgo_halo,vz_virgo_halo
      real(kind=xp), dimension(3,3):: rx,ry,rz
      integer:: ncell
      real(kind=xp), dimension(3):: mw,fil,rotaxis,virgo


      write(*,*) "open datafile"
      !open(1,file='virgo_xyz_hydro_l19_gal_clean_m1e8.5.dat',form='unformatted')
      open(1,file='virgo_xyz_dm_high_res.dat',form='unformatted')

      read(1) ncell
      allocate(pos(3,ncell),v(3,ncell),m(ncell))

      read(1) pos(1,:)
      read(1) pos(2,:)
      read(1) pos(3,:)
      read(1) v(1,:)
      read(1) v(2,:)
      read(1) v(3,:)
      read(1) m

      close(1)


      !stop




      virgo=(/0.48461068,0.50809848,0.49687076/)
      fil=(/0.497,0.4984,0.4996/)
      mw=(/0.5,0.5,0.5/)

      !pos(:,1)=(/0,0,0/)
      pos(:,1)=(fil(:)-mw(:))*(unit_l/3.08567758128E21)

      write(*,*) "pos 1",pos(:,1)

      rotaxis=fil

      vx_virgo_halo = -507.8579
      vy_virgo_halo = 229.4530
      vz_virgo_halo = -136.9451

      v(1,:) = v(1,:) - vx_virgo_halo
      v(2,:) = v(2,:) - vy_virgo_halo
      v(3,:) = v(3,:) - vz_virgo_halo


      pos(1,:)=pos(1,:)-(virgo(1)-mw(1))*(unit_l/3.08567758128E21)
      pos(2,:)=pos(2,:)-(virgo(2)-mw(2))*(unit_l/3.08567758128E21)
      pos(3,:)=pos(3,:)-(virgo(3)-mw(3))*(unit_l/3.08567758128E21)



      !pos(:,1)=(/0,0,0/)


      !mw_x=0.5-0.48461068
      !mw_y=0.5-0.50809848
      !mw_z=0.5-0.49687076

      !fil_x=0.5-0.497
      !fil_y=0.5-0.4984
      !fil_z=0.5-0.4996



      zero=0.0
      un=1.0

      !phi=atan(mw_x/mw_z)
      !phi=atan(mw_y/mw_z)
      !theta=atan(mw_x/mw_y)
      !theta=atan(fil_x/fil_y)

      theta=atan((rotaxis(1)-virgo(1))/(rotaxis(2)-virgo(2)))

      !write(*,*) "phi",phi
      write(*,*) "theta",theta

      !theta=pi/2
      !phi=0

      write(*,*) "mw",pos(:,1)

      c=cos(theta)
      s=sin(theta)

      write(*,*) "cos",c,"sin",s

      rz(:,1)=(/c,-s,zero/)
      rz(:,2)=(/s,c,zero/)
      rz(:,3)=(/zero,zero,un/)
      write(*,*) "rz",rz

      pos=transpose(matmul(transpose(pos),rz))
      v=transpose(matmul(transpose(v),rz))

      write(*,*) "mw",pos(:,1)

      phi=atan(pos(2,1)/pos(3,1))
      !phi=

      write(*,*) "phi",phi

      c=cos(phi)
      s=sin(phi)



      !ry(:,1)=(/c,zero,s/)
      !ry(:,2)=(/zero,un,zero/)
      !ry(:,3)=(/-s,zero,c/)

      rx(:,1)=(/un,zero,zero/)
      rx(:,2)=(/zero,c,-s/)
      rx(:,3)=(/zero,s,c/)

      pos=transpose(matmul(transpose(pos),rx))
      v=transpose(matmul(transpose(v),rx))

      phi=pi/2

      write(*,*) "phi",phi

      c=cos(phi)
      s=sin(phi)

      write(*,*) "cos",c,"sin",s



      ry(:,1)=(/c,zero,s/)
      ry(:,2)=(/zero,un,zero/)
      ry(:,3)=(/-s,zero,c/)

      rx(:,1)=(/un,zero,zero/)
      rx(:,2)=(/zero,c,-s/)
      rx(:,3)=(/zero,s,c/)

      !pos=transpose(matmul(transpose(pos),ry))

      write(*,*) "pos 1",pos(:,1)

      pos(1,:)=pos(1,:)+(virgo(1)-mw(1))*(unit_l/3.08567758128E21)
      pos(2,:)=pos(2,:)+(virgo(2)-mw(2))*(unit_l/3.08567758128E21)
      pos(3,:)=pos(3,:)+(virgo(3)-mw(3))*(unit_l/3.08567758128E21)

      write(*,*) "pos 1",pos(:,1)

      !stop

      write(*,*) "ecriture"

      open(2,file='virgo_xyz_dm_high_res_fil_los.dat',form='unformatted')
      write(2) ncell
      !write(2) ne
      !write(2) T
      !write(2) P
      write(2) pos(1,:)
      write(2) pos(2,:)
      write(2) pos(3,:)
      write(2) v(1,:)
      write(2) v(2,:)
      write(2) v(3,:)
      write(2) m
      !write(2) lvl
      close(2)

      write(*,*) "fin ecriture"

    end function rotation_DM

    recursive subroutine cell_division(x,y,z,P,n_e,T,m,vx,vy,vz,lvl,rescale_lvl,ndiv,xn,yn,zn,Pn,nn,Tn,mn,vxn,vyn,vzn,lvln)
        real(kind=xp),intent(in) :: x,y,z,P,n_e,T,m,vx,vy,vz,lvl
        real(kind=xp) :: len
        integer, intent(in) :: rescale_lvl,ndiv
        real(kind=xp), intent(out), dimension(ndiv) :: xn,yn,zn,Pn,nn,Tn,mn,vxn,vyn,vzn,lvln
        real(kind=xp), dimension(8) :: xi,yi,zi,Pi,ni,Ti,mi,vxi,vyi,vzi,lvli
        integer :: ncell,i,ndivn



        ncell=8
        !allocate(xi(ncell),yi(ncell),zi(ncell),Pi(ncell),ni(ncell),Ti(ncell),mi(ncell),vxi(ncell),vyi(ncell),vzi(ncell),lvli(ncell))
        !allocate(xn(ncell),yn(ncell),zn(ncell),Pn(ncell),nn(ncell),Tn(ncell),mn(ncell),vxn(ncell),vyn(ncell),vzn(ncell),lvln(ncell))
            len=0.25*(737441/(2**lvl))
            !len = 0.25*(1/2**(lvl-1))
            !write(*,*) "len",len
            xn(1)=x-len
            xn(2)=x+len
            xn(3)=x-len
            xn(4)=x+len
            xn(5)=x-len
            xn(6)=x+len
            xn(7)=x-len
            xn(8)=x+len

            yn(1)=y+len
            yn(2)=y+len
            yn(3)=y-len
            yn(4)=y-len
            yn(5)=y+len
            yn(6)=y+len
            yn(7)=y-len
            yn(8)=y-len

            zn(1)=z-len
            zn(2)=z-len
            zn(3)=z-len
            zn(4)=z-len
            zn(5)=z+len
            zn(6)=z+len
            zn(7)=z+len
            zn(8)=z+len

            Pn(:)=P
            Tn(:)=T
            nn(:)=n_e
            mn(:)=m/8
            vxn(:)=vx
            vyn(:)=vy
            vzn(:)=vz
            lvln(:)=lvl+1
            !write(*,*) "x",xn(1:8)
            !write(*,*) "ytour 1",yn(1:8)
            !write(*,*) "ztour 1",zn(1:8)
            !write(*,*) "lvl check",lvln(1)
            !write(*,*) "rescale_level", rescale_lvl
        if (lvln(1)<rescale_lvl) then
            !write(*,*) "lvl check",lvln(1)
            !write(*,*) "rescale_level", rescale_lvl
            xi(:)=xn(1:8)
            yi(:)=yn(1:8)
            zi(:)=zn(1:8)
            Pi(:)=Pn(1:8)
            ni(:)=nn(1:8)
            Ti(:)=Tn(1:8)
            mi(:)=mn(1:8)
            vxi(:)=vxn(1:8)
            vyi(:)=vyn(1:8)
            vzi(:)=vzn(1:8)
            lvli(:)=lvln(1:8)
            !deallocate(xn,yn,zn,Pn,nn,Tn,mn,vxn,vyn,vzn,lvln)
            ncell=ncell*(8**(rescale_lvl-lvln(1)))
            !allocate(xn(ncell),yn(ncell),zn(ncell),Pn(ncell),nn(ncell),Tn(ncell),mn(ncell),vxn(ncell),vyn(ncell),vzn(ncell),lvln(ncell))
            ndivn=(2**(rescale_lvl-lvli(1)))**3
            !write(*,*) "ndivn",ndivn
               do i=1,8
                   !write(*,*) "i",i

                   call cell_division(xi(i),yi(i),zi(i),Pi(i),ni(i),Ti(i),mi(i),vxi(i),vyi(i),vzi(i),lvli(i),rescale_lvl,ndivn,xn(1+ndivn*(i-1):ndivn*i),yn(1+ndivn*(i-1):ndivn*i),zn(1+ndivn*(i-1):ndivn*i),Pn(1+ndivn*(i-1):ndivn*i),nn(1+ndivn*(i-1):ndivn*i),Tn(1+ndivn*(i-1):ndivn*i),mn(1+ndivn*(i-1):ndivn*i),vxn(1+ndivn*(i-1):ndivn*i),vyn(1+ndivn*(i-1):ndivn*i),vzn(1+ndivn*(i-1):ndivn*i),lvln(1+ndivn*(i-1):ndivn*i))
                   !write(*,*) "xn",xn(:)
               enddo
        endif

    end subroutine cell_division


    function rescale_data(filein,fileout,rescale_lvl,savetype) result(ncell)

        character(len=*)::filein,fileout,savetype
        integer :: rescale_lvl,ncell,i,ndiv,n
        integer(8) :: nnew,nline
        integer::line_test
        real(kind=xp) :: nnew_r,ratio,ratiobis
        real(kind=xp), dimension(:), allocatable :: x,y,z,P,n_e,T,m,vx,vy,vz,lvl,xn,yn,zn,Pn,nn,Tn,mn,vxn,vyn,vzn,lvln,xi,yi,zi,Pi,ni,Ti,mi,vxi,vyi,vzi,lvli
        integer, dimension(:), allocatable::liminf,limsup

        open(1,file=filein,form='unformatted')
        read(1) ncell
        write(*,*) "ncell",ncell
        allocate(n_e(ncell),T(ncell),P(ncell),x(ncell),y(ncell),z(ncell),vx(ncell),vy(ncell),vz(ncell),m(ncell),lvl(ncell))
        write(*,*) "n_e"
        read(1) n_e
        write(*,*) "T"
        read(1) T
        write(*,*) "P"
        read(1) P
        write(*,*) "x"
        read(1) x
        write(*,*) "y"
        read(1) y
        write(*,*) "z"
        read(1) z
        write(*,*) "vx"
        read(1) vx
        write(*,*) "vy"
        read(1) vy
        write(*,*) "vz"
        read(1) vz
        write(*,*) "m"
        read(1) m
        write(*,*) "lvl"
        read(1) lvl

        close(1)

        n=0

        !do i=1,ncell
        !    if (lvl(i)<rescale_lvl) then
        !        n=n+1
        !    end if
        !end do

        !write(*,*) "n",n, "ncell",ncell,"ratio",real(n)/real(ncell)

        !stop
        !!!!test

        !ncell = 1
        !allocate(n_e(ncell),T(ncell),P(ncell),x(ncell),y(ncell),z(ncell),vx(ncell),vy(ncell),vz(ncell),m(ncell),lvl(ncell))

        !n_e(1) = 1.0
        !T(1) = 1.0
        !P(1) = 1.0
        !x(1) = 0.5
        !y(1) = 0.5
        !z(1) = 0.5
        !vx(1) = 0.0
        !vy(1) = 0.0
        !vz(1) = 0.0
        !m(1) = 1.0
        !lvl(1) = 1

        !nnew=int((0.03*2**rescale_lvl)**3)
        nnew_r=1e10
        write(*,*) "nnewr",nnew_r
        nnew=int8(nnew_r)
        !nnew=10**10
        write(*,*) "nnew",nnew

        allocate(nn(nnew),Tn(nnew),Pn(nnew),xn(nnew),yn(nnew),zn(nnew),vxn(nnew),vyn(nnew),vzn(nnew),mn(nnew),lvln(nnew))

        write(*,*) "starting rescaling"

        n=0
        do i=1,ncell
            if (lvl(i)==rescale_lvl) then
                nn(n)=n_e(i)
                Tn(n)=T(i)
                Pn(n)=P(i)
                xn(n)=x(i)
                yn(n)=y(i)
                zn(n)=z(i)
                vxn(n)=vx(i)
                vyn(n)=vy(i)
                vzn(n)=vz(i)
                mn(n)=m(i)
                lvln(n)=lvl(i)
                n=n+1
            else
                !write(*,*) "in"
                !write(*,*) "x(i)", x(i)
                !write(*,*) "y(i)", y(i)
                !write(*,*) "z(i)", z(i)
                ndiv=(2**(rescale_lvl-lvl(i)))**3
                !write(*,*) lvl(i),ndiv
                !write(*,*) "lvl",lvl(i)
                !write(*,*) "ndiv",ndiv
                call cell_division(x(i),y(i),z(i),P(i),n_e(i),T(i),m(i),vx(i),vy(i),vz(i),lvl(i),rescale_lvl,ndiv,xn(n:n+ndiv),yn(n:n+ndiv),zn(n:n+ndiv),Pn(n:n+ndiv),nn(n:n+ndiv),Tn(n:n+ndiv),mn(n:n+ndiv),vxn(n:n+ndiv),vyn(n:n+ndiv),vzn(n:n+ndiv),lvln(n:n+ndiv))
                !write(*,*) "xnew",xn
                !write(*,*) "ynew",yn
                !write(*,*) "znew",zn
                !write(*,*) xn(i:i+ndiv-1)
                !write(*,*) yn(i:i+ndiv-1)
                !write(*,*) zn(i:i+ndiv-1)
                n=n+ndiv
                !write(*,*) "n",n
                !stop

                !stop
            end if

        ratio=real(i)/real(ncell)
        !ratiobis=real(n)/949978046.0
        !write(*,*) "ne",nn(1:n)
        !stop

         if(modulo(i,1000000)==0) then
            write(*,*) "ratio", ratio, "i", i !, "lvl",lvl(k)
            !write(*,*) "max map", maxval(map), "min map", minval(map)
            !stop
         end if



        end do





        write(*,*) "nbr of cells", n-1

        nline=int(n/2.5e8)+1

        write(*,*) "nline", nline





        write(*,*) "xn(n-2:n+2)",xn(n-2:n+2)
        write(*,*) "nn(n-2:n+2)",nn(n-2:n+2)
        !write(*,*) "xn(n)",xn(n)
        !write(*,*) "xn(n+1)",xn(n+1)

        !stop


        write(*,*) "line_test",line_test
        !line_test=0

         open(2,file=fileout,form='unformatted')

        if (line_test==1) then
            write(2) nline
            allocate(liminf(nline),limsup(nline))
            do i=1,nline
            liminf(i)=1+(i-1)*250000000
            limsup(i)=250000000*i
        end do
            !write(*,*) "inf",liminf,"sup",limsup
            !stop
            do i=1,nline

            if (i<nline) then
                write(2) nn(liminf(i):limsup(i))
                !write(*,*) "i",i,"inf",int(1.0+(i-1)*2.5e8),"sup",i*2.5e8
            else
                write(2) nn(liminf(i):n-1)
                !write(*,*) "i",i,"inf",int(1.0+(i-1)*2.5e8),"sup",n-1
            end if
        end do
            !stop
            do i=1,nline
            if (i<nline) then
                write(2) Tn(liminf(i):limsup(i))
            else
                write(2) Tn(liminf(i):n-1)
            end if
        end do
            do i=1,nline
            if (i<nline) then
                write(2) Pn(liminf(i):limsup(i))
            else
                write(2) Pn(liminf(i):n)
            end if
        end do
            do i=1,nline
            if (i<nline) then
                write(2) xn(liminf(i):limsup(i))
            else
                write(2) xn(liminf(i):n)
            end if
        end do
            do i=1,nline
            if (i<nline) then
                write(2) yn(liminf(i):limsup(i))
            else
                write(2) yn(liminf(i):n)
            end if
        end do
            do i=1,nline
            if (i<nline) then
                write(2) zn(liminf(i):limsup(i))
            else
                write(2) zn(liminf(i):n)
            end if
        end do
            do i=1,nline
            if (i<nline) then
                write(2) vxn(liminf(i):limsup(i))
            else
                write(2) vxn(liminf(i):n)
            end if
        end do
            do i=1,nline
            if (i<nline) then
                write(2) vyn(liminf(i):limsup(i))
            else
                write(2) vyn(liminf(i):n)
            end if
        end do
            do i=1,nline
            if (i<nline) then
                write(2) vzn(liminf(i):limsup(i))
            else
                write(2) vzn(liminf(i):n)
            end if
        end do
            do i=1,nline
            if (i<nline) then
                write(2) mn(liminf(i):limsup(i))
            else
                write(2) mn(liminf(i):n)
            end if
        end do
            do i=1,nline
            if (i<nline) then
                write(2) lvln(liminf(i):limsup(i))
            else
                write(2) lvln(liminf(i):n)
            end if
        end do
        else
            if (savetype=='all') then
                write(*,*) "writing rescaled data"

                write(2) n-1
                write(2) nn(1:n-1)
                write(*,*) "n_e saved"
                write(2) Tn(1:n-1)
                write(*,*) "T saved"
                write(2) Pn(1:n-1)
                write(*,*) "P saved"
                write(2) xn(1:n-1)
                write(*,*) "x saved"
                write(2) yn(1:n-1)
                write(*,*) "y saved"
                write(2) zn(1:n-1)
                write(*,*) "z saved"
                write(2) vxn(1:n-1)
                write(*,*) "vx saved"
                write(2) vyn(1:n-1)
                write(*,*) "vy saved"
                write(2) vzn(1:n-1)
                write(*,*) "vz saved"
                write(2) mn(1:n-1)
                write(*,*) "m saved"
                write(2) lvln(1:n-1)
                write(*,*) "lvl saved"

            else if (savetype=='vel') then

                write(2) n-1
                write(2) xn(1:n-1)
                write(*,*) "x saved"
                write(2) yn(1:n-1)
                write(*,*) "y saved"
                write(2) zn(1:n-1)
                write(*,*) "z saved"
                write(2) vxn(1:n-1)
                write(*,*) "vx saved"
                write(2) vyn(1:n-1)
                write(*,*) "vy saved"
                write(2) vzn(1:n-1)
                write(*,*) "vz saved"

            !else if (savetype=='velmat') then

            end if

        end if

        close(2)

        write(*,*) "data saved"



    end function rescale_data

    function create_map_bis(proj,pxlvl,unit,weighting,filein,fileout) result(cen_x)


        real(kind=xp), dimension(11)::a
        real(kind=xp), dimension(:,:),allocatable :: map,weights
        real(kind=xp) :: cen_x,cen_y,cen_z,minx,maxx,miny,maxy,minz,maxz,lenx,leny,lenz,dx,dy,dz,dist,volume,ratio,i,j,l,start,finish,size_box,x_cen,y_cen,z_cen

        real(kind=xp), dimension(:), allocatable :: x,y,z,P,n_e,T,m,vx,vy,vz,lvl,r
        integer ::k,distx,disty,xinf,xsup,yinf,ysup,ix,iy,test,nx,ny,nz,ncell,n,pxlvl,ierror
        integer :: xtest,ytest,ztest,bl,br,tl,tr
        character(len=*)::filein,fileout
        character(len=2)::weighting
        character(len=1)::proj,unit
        real(kind=xp), dimension(21) :: m_rad,n_rad



        write(*,*) "projection: ",proj
        write(*,*)  "unit: ",unit


        open(1,file=filein,form='unformatted')
        read(1) ncell
        write(*,*) "ncell",ncell
        allocate(n_e(ncell),T(ncell),P(ncell),x(ncell),y(ncell),z(ncell),vx(ncell),vy(ncell),vz(ncell),m(ncell),lvl(ncell),r(ncell))
        read(1) n_e
        n_e=n_e/0.864 !conversion to electron density
        write(*,*) "n_e",n_e(1:10)
        read(1) T
        T=T*(kb/1.602e-16) ! conversion to kev instead of kelvin
        read(1) P
        P=P/10
        P=P/(1.602e-10)
        P=P*(0.76/0.864) ! conversion from Pa in cgs system to kev.cm-3 in SI units + conversion to electron pressure (this line)
        read(1) x
        read(1) y
        read(1) z
        read(1) vx
        read(1) vy
        read(1) vz
        read(1) m
        read(1) lvl

        write(*,*) "min lvl", minval(lvl(:))
        !write(*,*) "sum(m)",sum(m)
        close(1)

        !x_cen = 0.48461068
        !y_cen = 0.50809848
        !z_cen = 0.49687076

        !x_cen = (x_cen - 0.5) * (unit_l / 3.08567758128E21)
        !y_cen = (y_cen - 0.5) * (unit_l / 3.08567758128E21)
        !z_cen = (z_cen - 0.5) * (unit_l / 3.08567758128E21)


        !r = log10(sqrt((x - x_cen) ** 2 + (y - y_cen) ** 2 + (z - z_cen) ** 2))
        !m_rad=0
        !n_rad=0

        !do i=1,21
        !   write(*,*) i
  !   do j=1,size(m)
  !      if(r(j)>(i+17.5)*0.1 .and. r(j)<(i+18.5)*0.1) then
  !         m_rad(i)=m_rad(i)+m(j)
  !         n_rad(i)=n_rad(i)+1
  !      endif
  !   enddo
  !enddo

  !open(3,file='/data/cluster/tlebeau/virgo/rad_log_baryons_21.dat',form='unformatted')
  !write(3) m_rad
  !write(3) n_rad

  !write(*,*) m_rad
  !write(*,*) n_rad

  !  close(3)
  !  stop


        cen_x=0.48461068    !Virgo
        cen_y=0.50809848
        cen_z=0.49687076

        cen_x= (cen_x-0.5)*(unit_l/3.08567758128E21)
        cen_y= (cen_y-0.5)*(unit_l/3.08567758128E21)
        cen_z= (cen_z-0.5)*(unit_l/3.08567758128E21)

        write(*,*) "cen_x",cen_x

    !dx=737441.0/2**pxlvl
    !lenx=737441.0*0.03
        size_box=737441 !737433.333
        dx=size_box/2**pxlvl
        lenx=size_box*0.03
        nx=int(lenx/dx)
        dy=dx
        dz=dx
        leny=lenx
        lenz=lenx
        ny=nx
        nz=nx
        write(*,*) "dx",dx,"lenx",lenx,"nx",nx
    !stop

        minx=cen_x-lenx/2
        miny=cen_y-leny/2
        minz=cen_z-lenz/2

        allocate(map(nx,ny))
        allocate(weights(nx,ny))

        map=0
        weights=0

        open(2,file=fileout,form='unformatted')

        write(*,*) "ncell",ncell
        write(*,*) "len m",size(m)

        bl=0
        br=0
        tl=0
        tr=0

      !call MPI_INIT(ierror)

    do k=1,size(m)
      xtest=0
      ytest=0
      ztest=0
      !call cpu_time(start)
     if (proj=="x") then
        i=(y(k)-miny)/dy
        j=(z(k)-minz)/dz
        l=(x(k)-minx)/dx
     endif

     if (proj=="y") then
        i=(x(k)-minx)/dx
        j=(z(k)-minz)/dz
        l=(y(k)-miny)/dy
     endif

     if (proj=="z") then
        i=(x(k)-minx)/dx
        j=(y(k)-miny)/dy
        l=(z(k)-minz)/dz
        !write(*,*) "i real",(x(k)-minx)/dx
        !write(*,*) "j real",(y(k)-miny)/dy
     endif





     if(i>0 .and. j>0 .and. i<nx+1 .and. j<ny+1) then! .and. T(k)>8.6142e-3) then !.and. l>int(nx/2)-71 .and. l<int(nx/2)+71
         distx=2**(pxlvl-lvl(k)-1)
         disty=distx
         dist=size_box/2**(lvl(k))
         !write(*,*) "test"
         if (int(lvl(k))==pxlvl) then
             !write(*,*) "yes"
             xinf=int(i)
             xsup=xinf
             yinf=int(j)
             ysup=yinf
         else
             !if((i-int(i))>0.5) then
             !    xinf=int(i)-distx
             !    xsup=int(i)+(distx-1)
             !    xtest=1
             !else
             !    xinf=int(i)-(distx-1)
             !    xsup=int(i)+distx
             !end if
             !if((j-int(j))>0.5) then
             !    yinf=int(j)-disty
             !    ysup=int(j)+(disty-1)
             !    ztest=xtest+2
             !else
             !    yinf=int(j)-(disty-1)
             !    ysup=int(j)+disty
             !    ztest=xtest
             !end if
             xinf=int(i)-(distx-1)
             xsup=int(i)+distx
             yinf=int(j)-(disty-1)
             ysup=int(j)+disty

         end if
         !if (ztest==0) then
         !    bl=bl+1
         !else if (ztest==1) then
         !    br=br+1
         !else if (ztest==2) then
         !    tl=tl+1
         !else
         !    tr=tr+1
         !end if
         xinf=max(xinf,1)
         xsup=min(xsup,nx)
         yinf=max(yinf,1)
         ysup=min(ysup,ny)

         ratio=real(k)/real(size(m))
         if(modulo(k,1000000)==0) then
            write(*,*) "ratio", ratio, "k", k, "lvl",lvl(k)
         end if
         !write(*,*) map(xinf:xsup,yinf:ysup)
         !map(xinf:xsup,yinf:ysup)=map(xinf:xsup,yinf:ysup)+P(k)*m(k)
         !weights(xinf:xsup,yinf:ysup)=weights(xinf:xsup,yinf:ysup)+m(k)

         !write(*,*) "test"
         !weights(ix,iy)=weights(ix,iy)+m(k)
         !malwrite(*,*) map(xinf:xsup,yinf:ysup)
         !weights(ix,iy)=weights(ix,iy)+m(k)
         !stop
         do ix=xinf,xsup
             do iy=yinf,ysup
                 if(unit=="T") then
                 if (weighting=="mw") then
                    map(ix,iy)=map(ix,iy)+T(k)*m(k) !n_e(k) !P(k)*m(k)
                    weights(ix,iy)=weights(ix,iy)+m(k)
                 else if (weighting=="sl") then
                    !map(ix,iy)=map(ix,iy)+T(k)*m(k)*n_e(k)*T(k)**(-0.75)
                    !weights(ix,iy)=weights(ix,iy)+m(k)*n_e(k)*T(k)**(-0.75)
                    map(ix,iy)=map(ix,iy)+T(k)*n_e(k)*n_e(k)*T(k)**(-0.75)
                    weights(ix,iy)=weights(ix,iy)+n_e(k)*n_e(k)*T(k)**(-0.75)

                 else if (weighting=="ew") then
                    !map(ix,iy)=map(ix,iy)+T(k)*m(k)*n_e(k)*T(k)**(0.5)
                    !weights(ix,iy)=weights(ix,iy)+m(k)*n_e(k)*T(k)**(0.5)
                    map(ix,iy)=map(ix,iy)+T(k)*n_e(k)*n_e(k)*T(k)**(0.5)
                    weights(ix,iy)=weights(ix,iy)+n_e(k)*n_e(k)*T(k)**(0.5)
                 endif

                else if(unit=="n") then
                 map(ix,iy)=map(ix,iy)+n_e(k)*dist*pc*1e3
                else if(unit=="P") then
                 !write(*,*) "in P loop"
                 map(ix,iy)=map(ix,iy)+P(k)*m(k)
                 !if (weights(ix,iy)>real(2**31)) then
                     !write(*,*) "map(ix,iy) superior to 2**31"
                     !stop
                 !end if
                 !if (k==42373123) then
                     !write(*,*) "ix",ix,"iy",iy,"m(k)",m(k)
                     !write(*,*) "weights(ix,iy)",weights(ix,iy)
                 !end if
                 weights(ix,iy)=weights(ix,iy)+m(k)
                else if(unit=="y") then
                  map(ix,iy)=map(ix,iy)+P(k)*dist*pc*1e3*1.602e-10 ! multiplication by pc*1e3 to convert the distance from kpc to m, multiplication by 1.602e-10 to convert the electron pressure from kev.cm-3 to Pa (still in SI units), the compton parameter is dimensionless
                else if(unit=="v") then
                 map(ix,iy)=map(ix,iy)+vz(k)
                 weights(ix,iy)=weights(ix,iy)+1
                else if(unit=="vd") then
                 map(ix,iy)=map(ix,iy)+sqrt(vx(k)**2+vy(k)**2)
                 weights(ix,iy)=weights(ix,iy)+1
                else if(unit=="m") then ! emission measure
                  !write(*,*) "in EM loop"
                  map(ix,iy)=map(ix,iy)+((n_e(k)**2)/0.864)*(dist*pc*1e3*1e2)**3 ! multiplication by pc*1e3 to convert the distance from kpc to m and *1e2 to convert in cm
                  !write(*,*) "ne**2", (n_e(k)**2)/0.864
                  !write(*,*) "vol", (dist*pc*1e3*1e2)**3
                  !write(*,*) "map",map(ix,iy)
                  !stop
                else if(unit=="s") then !surface density
                  map(ix,iy)=map(ix,iy)+ m(k)/(dist)**2




                endif

           enddo
        enddo

     endif

      !call cpu_time(finish)
      !print '("Time = ",f6.3," seconds.")',finish-start
      !stop
    enddo

  !call MPI_FINALIZE(ierror)



    if(unit=="P" .or. unit=="T" .or. unit=="v") then
      !write(*,*) "max map", maxval(map(:,:))
      !write(*,*) "max weights", maxval(weights(:,:))
      map(:,:)=map(:,:)/weights(:,:)
    endif

    if(unit=="y") then
      map(:,:)=map(:,:)*(sigma_t/(me*c**2))
    end if

  !write(*,*) "bl",bl,"br",br,"tl",tl,"tr",tr
  !stop
    write(*,*) "rangement fini"
    write(*,*) "max map", maxval(map)
  !stop

  !do i=1,nx
    write(2) map
  !enddo

  !write(*,*) "map", map(:,:)

  !write(*,*) map

  !write(2,*) nx,ny

    close(2)

    deallocate(n_e,T,P,x,y,z,vx,vy,vz,m,lvl,map)


    end function create_map_bis

    function create_map_ter(proj,pxlvl,unit,weighting,filein,fileout,datatype,filament,map_width_in_mpc,tcut) result(cen_x)


        real(kind=xp), dimension(11)::a
        real(kind=xp), dimension(:,:),allocatable :: map,weights,map2
        real(kind=xp) :: cen_x,cen_y,cen_z,minx,maxx,miny,maxy,minz,maxz,lenx,leny,lenz,dx,dy,dz,dist,volume,ratio,i,j,l,start,finish,size_box,x_cen,y_cen,z_cen,factor,mean,wmean,vx_virgo_halo,vy_virgo_halo,vz_virgo_halo,tcut

        real(kind=xp), dimension(:), allocatable :: x,y,z,P,n_e,T,m,vx,vy,vz,lvl,r,rho
        integer ::k,distx,disty,xinf,xsup,yinf,ysup,ix,iy,test,nx,ny,nz,ncell,n,pxlvl,ierror,ni,nj,nl,zsup,zinf,distz
        integer :: xtest,ytest,ztest,bl,br,tl,tr,filament,in_map,nbr_cells_in_map,nmean,nnan,nslice,map_width_in_mpc
        character(len=*)::filein,fileout
        character(len=*)::weighting
        character(len=2)::datatype
        character(len=1)::proj,unit
        real(kind=xp), dimension(21) :: m_rad,n_rad



        write(*,*) "projection: ",proj
        write(*,*)  "unit: ",unit
        !write(*,*) "datatype:",datatype
        write(*,*) "weigthing:",weighting
        !stop


        open(1,file=filein,form='unformatted')
        read(1) ncell
        !write(*,*) "ncell",ncell

        !ncell = ncell-1

        if (datatype=="hy") then

        write(*,*) "datatype:",datatype

        allocate(n_e(ncell),T(ncell),P(ncell),x(ncell),y(ncell),z(ncell),vx(ncell),vy(ncell),vz(ncell),m(ncell),lvl(ncell),r(ncell))
        read(1) n_e
        n_e=n_e/0.864 !conversion to electron density
        !write(*,*) "n_e",n_e(1:10)
        read(1) T
        T=T*(kb/1.602e-16)! conversion to kev instead of kelvin
        !write(*,*) "max",maxval(T),"min",minval(T)
        !stop
        read(1) P
        P=P/10
        P=P/(1.602e-10)
        P=P*(0.76/0.864) ! conversion from Pa in cgs system to kev.cm-3 in SI units + conversion to electron pressure (this line)
        read(1) x
        read(1) y
        read(1) z
        read(1) vx
        read(1) vy
        read(1) vz
        read(1) m
        read(1) lvl

        write(*,*) "data read"

        vx_virgo_halo = -507.8579
        vy_virgo_halo = 229.4530
        vz_virgo_halo = -136.9451

        vx = vx - vx_virgo_halo
        vy = vy - vy_virgo_halo
        vz = vz - vz_virgo_halo



        if (unit=="f") then

            allocate(rho(ncell))
            rho = n_e*0.864 * 1.66D-24 * 0.76 !conversion from electron density in cm-3 to mass density in g.cm-3
            rho = rho * 1e3 !conversion from g.cm-3 to kg.m-3
            !vx_virgo = -509.1301
            !vy_virgo = 228.9488
            !vz_virgo = -131.9249
            !vx = vx - vx_virgo
            !vy = vy - vy_virgo
            !vz = vz - vz_virgo
            !vx = vx * 1e3 !conversion from km.s-1 to m.s-1
            !vy = vy * 1e3 !conversion from km.s-1 to m.s-1
            !vz = vz * 1e3 !conversion from km.s-1 to m.s-1
        end if

        !write(*,*) "min lvl", minval(lvl(:))
        !write(*,*) "sum(m)",sum(m)

        !write(*,*) "x",x(1:10)

        close(1)

        else if (datatype=="dm") then

        write(*,*) "datatype:",datatype

        allocate(x(ncell),y(ncell),z(ncell),vx(ncell),vy(ncell),vz(ncell),m(ncell),lvl(ncell))
        read(1) x
        read(1) y
        read(1) z
        read(1) vx
        read(1) vy
        read(1) vz
        read(1) m
        lvl = pxlvl

        !write(*,*) "m",m(1:20)

        endif

        !stop


        !x_cen = 0.48461068
        !y_cen = 0.50809848
        !z_cen = 0.49687076

        !x_cen = (x_cen - 0.5) * (unit_l / 3.08567758128E21)
        !y_cen = (y_cen - 0.5) * (unit_l / 3.08567758128E21)
        !z_cen = (z_cen - 0.5) * (unit_l / 3.08567758128E21)


        !r = log10(sqrt((x - x_cen) ** 2 + (y - y_cen) ** 2 + (z - z_cen) ** 2))
        !m_rad=0
        !n_rad=0

        !do i=1,21
        !   write(*,*) i
  !   do j=1,size(m)
  !      if(r(j)>(i+17.5)*0.1 .and. r(j)<(i+18.5)*0.1) then
  !         m_rad(i)=m_rad(i)+m(j)
  !         n_rad(i)=n_rad(i)+1
  !      endif
  !   enddo
  !enddo

  !open(3,file='/data/cluster/tlebeau/virgo/rad_log_baryons_21.dat',form='unformatted')
  !write(3) m_rad
  !write(3) n_rad

  !write(*,*) m_rad
  !write(*,*) n_rad

  !  close(3)
  !  stop


        cen_x=0.48461068    !Virgo
        cen_y=0.50809848
        cen_z=0.49687076

        !cen_x=0.48461053 !M87
        !cen_y=0.50809777
        !cen_z=0.49687093

        cen_x= (cen_x-0.5)*(unit_l/3.08567758128E21)
        cen_y= (cen_y-0.5)*(unit_l/3.08567758128E21)
        cen_z= (cen_z-0.5)*(unit_l/3.08567758128E21)

        !write(*,*) "cen_x",cen_x, "cen_y",cen_y, "cen_z",cen_z

    !dx=737441.0/2**pxlvl
    !lenx=737441.0*0.03
        size_box=737441 !737433.333
        dx=size_box/2**pxlvl
        lenx=size_box*0.03
        nx=int(lenx/dx)
        dy=dx
        dz=dx
        leny=lenx
        lenz=lenx
        ny=nx
        nz=nx
        !write(*,*) "dx",dx,"lenx",lenx,"nx",nx
    !stop

        !minx=cen_x-lenx/2
        !miny=cen_y-leny/2
        !minz=cen_z-lenz/2

        !write(*,*) "minx",minx,"miny",miny,"minz",minz

        if (filament==0) then


            if (proj=="x") then
                !lenx = lenx/221.23
                leny = leny/(22.123/map_width_in_mpc)
                lenz = lenz/(22.123/map_width_in_mpc)
            else if (proj=="y") then
                lenx = lenx/(22.123/map_width_in_mpc)
                !leny = leny/221.23
                lenz = lenz/(22.123/map_width_in_mpc)
            else if (proj=="z") then
                lenx = lenx/(22.123/map_width_in_mpc)
                leny = leny/(22.123/map_width_in_mpc)
                !lenz = lenz/2212.3
            end if


            minx=cen_x-lenx/2
            miny=cen_y-leny/2
            minz=cen_z-lenz/2

            nx=int(lenx/dx)
            ny=int(leny/dy)
            nz=int(lenz/dz)

            write(*,*) "dx",dx !,"lenx",lenx,"leny",leny,"lenz",lenz
            !write(*,*) "nx",nx,"ny",ny,"nz",nz
            !stop

            !allocate(map(nx,ny))
            !allocate(weights(nx,ny))

            if (proj=="y") then
                allocate(map(nz,nx))
                allocate(weights(nz,nx))
                if (unit=="d") then
                    allocate(map2(nz,nx))
                end if
            else if (proj=="z") then
                allocate(map(nx,ny))
                allocate(weights(nx,ny))
                if (unit=="d") then
                    allocate(map2(nx,ny))
                end if
            else
                allocate(map(nz,ny))
                allocate(weights(nz,ny))
                if (unit=="d") then
                    allocate(map2(nz,ny))
                end if
            end if

            write(*,*) "lenx",lenx,"leny",leny,"lenz",lenz

            write(*,*) "cen_x",cen_x,"cen_y",cen_y,"cen_z",cen_z

            write(*,*) "nx",nx,"ny",ny,"nz",nz

            !stop

        end if

        !!!!Modification for filament visualisation
        if (filament==1) then



            !lenx = lenx/2 ! = 5.53 Mpc
            !leny = leny/2 ! = 5.53 Mpc
            !lenz = lenz/256 ! = 11.06 Mpc
            !lenz = lenz/512
            !lenz = lenz/1024

            !!transverse slices

            !lenx = 1000
            !leny = leny/8
            !lenx = leny

            !lenz = lenz/8
            lenx = lenx/1024
            cen_x = cen_x - 2000


            write(*,*) "lenx",lenx,"leny",leny,"lenz",lenz

            !cen_z = cen_z - lenz * nslice
            !cen_x = cen_x + 0.75*lenx
            !cen_y = cen_y - 0.5*leny

            !cen_x = cen_x - lenx/2 ! * (nslice-20)*5

            write(*,*) "cen_x",cen_x,"cen_y",cen_y,"cen_z",cen_z

            minx = cen_x - lenx/2
            miny = cen_y - leny/2
            minz = cen_z - lenz/2

            write(*,*) "minx",minx,"miny",miny,"minz",minz

            nx=int(lenx/dx)
            ny=int(leny/dy)
            nz=int(lenz/dz)

            write(*,*) "nx",nx,"ny",ny,"nz",nz

            !stop

            if (proj=="y") then
                allocate(map(nz,nx))
                allocate(weights(nz,nx))
            else if (proj=="z") then
                allocate(map(nx,ny))
                allocate(weights(nx,ny))
            else
                allocate(map(nz,ny))
                allocate(weights(nz,ny))
            end if

        end if


        !cen_x = -5817 !-16878.75
        !cen_y = 441.25
        !cen_z = -2307

        !lenx = (size_box*0.03) / 2

        !nx=int(lenx/dx)
        !dy=dx
        !dz=dx
        !leny=lenx
        !lenz=lenx
        !ny=nx
        !nz=nx

        !minx=cen_x-lenx/2
        !miny=cen_y-leny/2
        !minz=cen_z-lenz/2

        map=0
        weights=0

        open(2,file=fileout,form='unformatted')

        write(*,*) "ncell",ncell
        write(*,*) "len m",size(m)

        bl=0
        br=0
        tl=0
        tr=0

      !call MPI_INIT(ierror)
    nbr_cells_in_map = 0

    do k=1,size(m)
      factor = 1
      in_map = 0
      xtest=0
      ytest=0
      ztest=0
      !call cpu_time(start)
     if (proj=="x") then
        i=(z(k)-minz)/dz
        j=(y(k)-miny)/dy
        l=(x(k)-minx)/dx

        ni=nz
        nj=ny
        nl=nx

     endif

     if (proj=="y") then

        i=(z(k)-minz)/dz
        j=(x(k)-minx)/dx
        l=(y(k)-miny)/dy

        ni=nz
        nj=nx
        nl=ny

     endif

     if (proj=="z") then
        i=(x(k)-minx)/dx
        j=(y(k)-miny)/dy
        l=(z(k)-minz)/dz
        !write(*,*) "i real",(x(k)-minx)/dx
        !write(*,*) "j real",(y(k)-miny)/dy

        ni=nx
        nj=ny
        nl=nz

     endif

     distx=2**(pxlvl-lvl(k)-1)
     disty=dist
     distz=distx
     dist=size_box/2**(lvl(k))
     !write(*,*) "test"
     if (int(lvl(k))>=pxlvl) then
             !write(*,*) "yes"
             xinf=int(i)
             xsup=xinf
             yinf=int(j)
             ysup=yinf
             zinf=int(l)
             zsup=zinf

     else
             !if((i-int(i))>0.5) then
             !    xinf=int(i)-distx
             !    xsup=int(i)+(distx-1)
             !    xtest=1
             !else
             !    xinf=int(i)-(distx-1)
             !    xsup=int(i)+distx
             !end if
             !if((j-int(j))>0.5) then
             !    yinf=int(j)-disty
             !    ysup=int(j)+(disty-1)
             !    ztest=xtest+2
             !else
             !    yinf=int(j)-(disty-1)
             !    ysup=int(j)+disty
             !    ztest=xtest
             !end if
             xinf=int(i)-(distx-1)
             xsup=int(i)+distx
             yinf=int(j)-(disty-1)
             ysup=int(j)+disty
             zinf=int(l)-(distz-1)
             zsup=int(l)+distz

     end if





     if(i>0 .and. j>0 .and. i<ni+1 .and. j<nj+1 .and. l>0 .and. l<nl+1) then !.and. l>int(nx/2)-71 .and. l<int(nx/2)+71

         !if (ztest==0) then
         !    bl=bl+1
         !else if (ztest==1) then
         !    br=br+1
         !else if (ztest==2) then
         !    tl=tl+1
         !else
         !    tr=tr+1
         !end if
         xinf=max(xinf,1)
         xsup=min(xsup,ni)
         yinf=max(yinf,1)
         ysup=min(ysup,nj)

         if (zinf<0) then
            zinf = 1
            factor = (real(zsup)-real(zinf))/(2*real(distz))
         end if

         if (zsup>nl) then
            zsup = nl
            factor = (real(zsup)-real(zinf))/(2*real(distz))
         end if

         !zinf=max(zinf,1)
         !zsup=min(zsup,nl)
         !factor = (real(zsup)-real(zinf))/(2*real(distz))

         !write(*,*) map(xinf:xsup,yinf:ysup)
         !map(xinf:xsup,yinf:ysup)=map(xinf:xsup,yinf:ysup)+P(k)*m(k)
         !weights(xinf:xsup,yinf:ysup)=weights(xinf:xsup,yinf:ysup)+m(k)

         !write(*,*) "test"
         !weights(ix,iy)=weights(ix,iy)+m(k)
         !malwrite(*,*) map(xinf:xsup,yinf:ysup)
         !weights(ix,iy)=weights(ix,iy)+m(k)
         !stop
         in_map = 1
     else
        if (xsup>0 .and. xinf<1) then
           xinf = 1
           !in_map = 1
        end if
        if (xinf<ni+1 .and. xsup>ni) then
           xsup = ni
           !in_map = 1
        end if
        if (ysup>0 .and. yinf<1) then
           yinf = 1
           !in_map = 1
        end if
        if (yinf<nj+1 .and. ysup>nj) then
           ysup = nj
           !in_map = 1
        end if
        if (zsup>0 .and. zinf<1) then
           zinf = 1
           !in_map = 1
           factor = (real(zsup)-real(zinf))/(2*real(distz))
           !write(*,*) "zsup",zsup,"zinf",zinf,"distx",distx,"factor",factor
           !stop
           !write(*,*) "factor",factor
        end if
        if (zinf<nl+1 .and. zsup>nl) then
           zsup = nl
           factor = (real(zsup)-real(zinf))/(2*real(distz))
           !in_map = 1
           !write(*,*) "factor",factor
           !write(*,*) "zsup",zsup,"zinf",zinf,"distx",distx,"factor",factor
           !stop
        end if

        if(xinf>0 .and. yinf>0 .and. xsup<ni+1 .and. ysup<nj+1 .and. zinf>0 .and. zsup<nl+1) then
           in_map = 1
        !   write(*,*) "in map"
           !write(*,*) "zsup",zsup,"zinf",zinf,"distz",distz,"factor",factor,"lvl",lvl(k),"i",i,"j",j,"l",l
           !stop
        !else
        !   in_map = 0
        end if

     end if

      if (in_map==1 .and. T(k)>tcut) then !.and. T(k)>1.7 .and. T(k)<12) then !
        !write(*,*) "T",T(k)
        !stop

        nbr_cells_in_map = nbr_cells_in_map + 1

        !write(*,*) "xinf",xinf,"xsup",xsup,"yinf",yinf,"ysup",ysup,"zinf",zinf,"zsup",zsup,"factor",factor

        ratio=real(k)/real(size(m))
         if(modulo(k,1000000)==0) then
            write(*,*) "ratio", ratio, "k", k, "lvl",lvl(k)
            !write(*,*) "max map", maxval(map), "min map", minval(map)
            !stop
         end if

        do ix=xinf,xsup
             do iy=yinf,ysup
                 if(unit=="T") then
                    if (weighting=="mw") then
                    map(ix,iy)=map(ix,iy)+T(k)*m(k)*factor !n_e(k) !P(k)*m(k)
                    weights(ix,iy)=weights(ix,iy)+m(k)*factor
                    !if (ix>27710 .and. ix<27720 .and. iy>29385 .and. iy<29395) then
                    !    write(*,*) "i",i,"j",j,"l",l,"xinf",xinf,"xsup",xsup,"yinf",yinf,"ysup",ysup,"factor",factor,'lvl',lvl(k)
                    !end if

                    else if (weighting=="sl") then
                    !map(ix,iy)=map(ix,iy)+T(k)*m(k)*n_e(k)*T(k)**(-0.75)
                    !weights(ix,iy)=weights(ix,iy)+m(k)*n_e(k)*T(k)**(-0.75)
                    map(ix,iy)=map(ix,iy)+T(k)*n_e(k)*n_e(k)*T(k)**(-0.75)*factor
                    weights(ix,iy)=weights(ix,iy)+n_e(k)*n_e(k)*T(k)**(-0.75)*factor

                    else if (weighting=="ew") then
                    !map(ix,iy)=map(ix,iy)+T(k)*m(k)*n_e(k)*T(k)**(0.5)
                    !weights(ix,iy)=weights(ix,iy)+m(k)*n_e(k)*T(k)**(0.5)
                    map(ix,iy)=map(ix,iy)+T(k)*n_e(k)*n_e(k)*T(k)**(0.5)
                    weights(ix,iy)=weights(ix,iy)+n_e(k)*n_e(k)*T(k)**(0.5)
                    endif

                 else if(unit=="n") then
                 map(ix,iy)=map(ix,iy)+n_e(k)*dist*pc*1e3*factor*1e2
                 !map(ix,iy)=map(ix,iy)+n_e(k)*m(k)*factor
                 !test for filament project maps
                 !map(ix,iy)=map(ix,iy)+n_e(k)*factor !n_e(k) !P(k)*m(k)
                 !weights(ix,iy)=weights(ix,iy)+factor*m(k)
                 weights(ix,iy)=weights(ix,iy)+factor

                else if(unit=="P") then
                 !write(*,*) "in P loop"
                 map(ix,iy)=map(ix,iy)+P(k)*m(k)*factor
                 !if (weights(ix,iy)>real(2**31)) then
                     !write(*,*) "map(ix,iy) superior to 2**31"
                     !stop
                 !end if
                 !if (k==42373123) then
                     !write(*,*) "ix",ix,"iy",iy,"m(k)",m(k)
                     !write(*,*) "weights(ix,iy)",weights(ix,iy)
                 !end if
                 weights(ix,iy)=weights(ix,iy)+m(k)*factor
                else if(unit=="y") then
                  map(ix,iy)=map(ix,iy)+P(k)*dist*pc*1e3*1.602e-10*factor ! multiplication by pc*1e3 to convert the distance from kpc to m, multiplication by 1.602e-10 to convert the electron pressure from kev.cm-3 to Pa (still in SI units), the compton parameter is dimensionless
                else if(unit=="v") then
                     if (weighting=="vxew") then
                         map(ix,iy)=map(ix,iy)+vx(k)*factor*n_e(k)**2
                         weights(ix,iy)=weights(ix,iy)+factor*n_e(k)**2
                         !map(ix,iy)=map(ix,iy)+vx(k)*factor*m(k)
                     else if (weighting=="vyew") then
                         map(ix,iy)=map(ix,iy)+vy(k)*factor*n_e(k)**2
                         weights(ix,iy)=weights(ix,iy)+factor*n_e(k)**2
                     else if (weighting=="vzew") then
                         map(ix,iy)=map(ix,iy)+vz(k)*factor*n_e(k)**2
                         weights(ix,iy)=weights(ix,iy)+factor*n_e(k)**2
                     else if (weighting=="vew") then
                         map(ix,iy)=map(ix,iy)+sqrt(vx(k)**2+vy(k)**2+vz(k)**2)*factor*n_e(k)**2
                         weights(ix,iy)=weights(ix,iy)+factor*n_e(k)**2

                     else if (weighting=="vxmw") then
                         map(ix,iy)=map(ix,iy)+vx(k)*factor*m(k)
                         weights(ix,iy)=weights(ix,iy)+factor*m(k)
                     else if (weighting=="vymw") then
                         map(ix,iy)=map(ix,iy)+vy(k)*factor*m(k)
                         weights(ix,iy)=weights(ix,iy)+factor*m(k)
                     else if (weighting=="vzmw") then
                         map(ix,iy)=map(ix,iy)+vz(k)*factor*m(k)
                         weights(ix,iy)=weights(ix,iy)+factor*m(k)
                     else if (weighting=="vmw") then
                         map(ix,iy)=map(ix,iy)+sqrt(vx(k)**2+vy(k)**2+vz(k)**2)*factor*m(k)
                         weights(ix,iy)=weights(ix,iy)+factor*m(k)
                     end if
                 !map(ix,iy)=map(ix,iy)+vx(k)*factor
                     !if (weighting=="vz") then
                     !    weights(ix,iy)=weights(ix,iy)+factor*n_e(k)**2
                     !else
                     !    weights(ix,iy)=weights(ix,iy)+factor
                     !end if

                 !weights(ix,iy)=weights(ix,iy)+factor*m(k)
                else if(unit=="d") then !velocity dispersion
                     !write(*,*) "test"
                     !stop
                    if (weighting=="vxew") then
                        map(ix,iy)=map(ix,iy)+vx(k)**2*factor*n_e(k)**2
                        map2(ix,iy)=map2(ix,iy)+vx(k)*factor*n_e(k)**2
                        weights(ix,iy)=weights(ix,iy)+factor*n_e(k)**2
                    else if (weighting=="vyew") then
                        map(ix,iy)=map(ix,iy)+vy(k)**2*factor*n_e(k)**2
                        map2(ix,iy)=map2(ix,iy)+vy(k)*factor*n_e(k)**2
                        weights(ix,iy)=weights(ix,iy)+factor*n_e(k)**2
                    else if (weighting=="vzew") then
                        map(ix,iy)=map(ix,iy)+vz(k)**2*factor*n_e(k)**2
                        map2(ix,iy)=map2(ix,iy)+vz(k)*factor*n_e(k)**2
                        weights(ix,iy)=weights(ix,iy)+factor*n_e(k)**2
                    else if (weighting=="vew") then
                        map(ix,iy)=map(ix,iy)+sqrt(vx(k)**2+vy(k)**2+vz(k)**2)**2*factor*n_e(k)**2
                        map2(ix,iy)=map2(ix,iy)+sqrt(vx(k)**2+vy(k)**2+vz(k)**2)*factor*n_e(k)**2
                        weights(ix,iy)=weights(ix,iy)+factor*n_e(k)**2
                        
                    else if (weighting=="vxmw") then
                        map(ix,iy)=map(ix,iy)+vx(k)**2*factor*m(k)
                        map2(ix,iy)=map2(ix,iy)+vx(k)*factor*m(k)
                        weights(ix,iy)=weights(ix,iy)+factor*m(k)
                    else if (weighting=="vymw") then
                        !write(*,*) "test"
                        !stop
                        map(ix,iy)=map(ix,iy)+vy(k)**2*factor*m(k)
                        map2(ix,iy)=map2(ix,iy)+vy(k)*factor*m(k)
                        weights(ix,iy)=weights(ix,iy)+factor*m(k)
                    else if (weighting=="vzmw") then
                        map(ix,iy)=map(ix,iy)+vz(k)**2*factor*m(k)
                        map2(ix,iy)=map2(ix,iy)+vz(k)*factor*m(k)
                        weights(ix,iy)=weights(ix,iy)+factor*m(k)
                    else if (weighting=="vmw") then
                        map(ix,iy)=map(ix,iy)+sqrt(vx(k)**2+vy(k)**2+vz(k)**2)**2*factor*m(k)
                        map2(ix,iy)=map2(ix,iy)+sqrt(vx(k)**2+vy(k)**2+vz(k)**2)*factor*m(k)
                        weights(ix,iy)=weights(ix,iy)+factor*m(k)
                    end if 
                    


                else if(unit=="m") then ! emission measure
                  !write(*,*) "in EM loop"
                  map(ix,iy)=map(ix,iy)+((n_e(k)**2)/0.864)*(dist*1e-3)*factor ! multiplication by 1e-3 to convert in cm-6*Mpc like Eckert et al. 2012 (before: multiplication by pc*1e3 to convert the distance from kpc to m and *1e2 to convert in cm)
                  !map(ix,iy)=map(ix,iy)+((n_e(k)**2)/0.864)*(dist*pc*1e3*1e2)**3 ! multiplication by 1e-3 to convert in cm-6*Mpc like Eckert et al. 2012 (before: multiplication by pc*1e3 to convert the distance from kpc to m and *1e2 to convert in cm)
                  !write(*,*) "ne**2", (n_e(k)**2)/0.864
                  !write(*,*) "vol", (dist*pc*1e3*1e2)**3
                  !write(*,*) "map",map(ix,iy)
                  !stop
                else if(unit=="s") then !surface density
                  map(ix,iy)=map(ix,iy)+ (m(k)*factor)/(2**((19-lvl(k))*2))

                else if(unit=="k") then
                     map(ix,iy)=map(ix,iy)+(T(k)/n_e(k)**(2/3))*m(k)*factor !n_e(k) !P(k)*m(k)
                     weights(ix,iy)=weights(ix,iy)+m(k)*factor

                else if(unit=="f") then
                     if (weighting=="vx") then
                         map(ix,iy)=map(ix,iy)+vx(k)*rho(k)*factor
                     else if (weighting=="vy") then
                         map(ix,iy)=map(ix,iy)+vy(k)*rho(k)*factor
                     else if (weighting=="vz") then
                         map(ix,iy)=map(ix,iy)+vz(k)*rho(k)*factor

                         end if
                     weights(ix,iy)=weights(ix,iy)+factor

                else if(unit=="a") then  ! compute mach number map
                     map(ix,iy)=map(ix,iy)+(sqrt(vx(k)**2+vy(k)**2+vz(k)**2)/sqrt(((5 / 3) * kb * T(k)) / (mu * mp)))*factor
                     weights(ix,iy)=weights(ix,iy)+factor

                else if(unit=="l") then !mean resolution in slice
                     map(ix,iy)=map(ix,iy)+lvl(k)
                     weights(ix,iy)=weights(ix,iy)+1




                endif

           enddo
        enddo
      end if


    end do

    if (filament==1) then

      write(*,*) 'calculating the mean neighbourhood value for Nan pixels'

      nnan = 0

      do i=1,ni
          do j=1,nj
              if (map(i,j)==0) then
                  nnan = nnan + 1
                  !write(*,*) 'pixel is zero'
                  mean = 0
                  nmean = 0
                  wmean = 0
                  do ix=-2,2
                      do iy=-2,2
                          if (map(i+ix,j+iy)/=0 .and. i+ix>0 .and. j+iy>0) then
                              mean = mean + map(i+ix,j+iy)
                              nmean = nmean + 1
                              wmean = wmean + weights(i+ix,j+iy)
                          end if
                      end do
                  end do
                  map(i,j) = mean/nmean
                  weights(i,j) = wmean/nmean
                  !write(*,*) 'map(i,j)',map(i,j),'weights(i,j)',weights(i,j)

              end if
          end do
      end do


    end if

      !write(*,*) 'number of nan pixels',nnan

      !stop



      !call cpu_time(finish)
      !print '("Time = ",f6.3," seconds.")',finish-start
      !stop


  !call MPI_FINALIZE(ierror)

    !write(*,*) "map",map(27710:27720,29391)



    if(unit=="P" .or. unit=="T" .or. unit=="v" .or. unit=="k" .or. unit=="f" .or. unit=="a" .or. unit=="l") then
      !write(*,*) "max map", maxval(map(:,:))
      !write(*,*) "max weights", maxval(weights(:,:))
      !write(*,*) "max map", maxval(map),'min map',minval(map)
      map(:,:)=map(:,:)/weights(:,:)
      !write(*,*) "max map", maxval(map),'min map',minval(map)

    !if(unit=="n") then !test for filament project maps
    !  map(:,:)=map(:,:)/weights(:,:)
    !end if
    endif

    if(unit=="y") then
      map(:,:)=map(:,:)*(sigma_t/(me*c**2))
    end if

    if (unit=="s") then
        map(:,:)=map(:,:)/(size_box/2**(pxlvl))
    end if

    if (unit=="d") then
        map2(:,:)=map2(:,:)/weights(:,:)
        map(:,:)=map(:,:)/weights(:,:) - map2(:,:)**2
        map(:,:)=sqrt(map(:,:))
    end if

    if (unit=="n") then
        map(:,:)=map(:,:)/weights(:,:)
    end if

  !write(*,*) "bl",bl,"br",br,"tl",tl,"tr",tr
  !stop
    !write(*,*) "map",map(27710:27720,29391)
    write(*,*) "rangement fini"
    write(*,*) "max map", maxval(map), "min map", minval(map)
    write(*,*) "nbr_cells_in_map",nbr_cells_in_map
  !stop

  !do i=1,nx
    write(2) nx,ny,nz
    write(2) cen_x,cen_y,cen_z
    write(2) map
  !enddo

  !write(*,*) "map", map(:,:)

  !write(*,*) map

  !write(2,*) nx,ny

    close(2)

    !if (datatype=="hy") then
    !    deallocate(n_e,T,P,x,y,z,vx,vy,vz,m,lvl,map)
    !else if (datatype=="dm") then
    !    deallocate(x,y,z,vx,vy,vz,m,lvl,map)
    !end if

    !stop


    end function create_map_ter

    function cut_lines_in_file(filein,fileout) result(i)
        character(len=*)::filein,fileout
        real(kind=xp), dimension(:), allocatable :: map
        integer(8) :: nline,i,n,npx
        integer(8), dimension(:), allocatable :: liminf,limsup

        open(1,file=filein,form='unformatted')
        npx=3958171396
        !allocate(map(3958171397))
        !npx=(15728*4)**2
        write(*,*) "npx",npx
        allocate(map(npx))
        write(*,*) "map reading"
        !map=0.0
        read(1) map

        write(*,*) "file read, map:"
        write(*,*) "size",size(map)
        write(*,*) "end", map(npx)

        !stop

        !close(1)

        !do i=1,size(map)
        !    if(map(i)<0.0) then
        !        n=i
        !    end if
        !end do

        !write(*,*) "n",n

        n=3958171396
        write(*,*) "n",n

        !write(*,*) map
        !write(*,*) "size map", size(map)
        nline=int(n/2.5e8)+1
        write(*,*) "nline",nline

        open(2,file=fileout,form='unformatted')

        write(2) nline
        allocate(liminf(nline),limsup(nline))

        do i=1,nline
            liminf(i)=1+(i-1)*250000000
            limsup(i)=250000000*i
        end do
            !write(*,*) "inf",liminf,"sup",limsup
            !stop
        do i=1,nline

            if (i<nline) then
                write(*,*) "i",i,"liminf",liminf(i),"limsup",limsup(i)
                write(*,*) "map test" , map(limsup(i):limsup(i)+2)
                !write(*,*) map(liminf(i):limsup(i))
                !stop
                write(2) map(liminf(i):limsup(i))
                !write(*,*) "i",i,"inf",int(1.0+(i-1)*2.5e8),"sup",i*2.5e8
            else
                write(*,*) "last line"
                write(*,*) "liminf",liminf(i),"end",n
                write(2) map(liminf(i):n)
                !write(*,*) "i",i,"inf",int(1.0+(i-1)*2.5e8),"sup",n-1
            end if
        end do

        close(1)
        close(2)

    end function cut_lines_in_file

    function gas_in_gal() result(i)
     real(kind=xp), dimension(6000) :: xcen,ycen,zcen,rvir,mgal
     !real(kind=xp), dimension(:), allocatable ::xdm,ydm,zdm,mdm,r
     real(kind=xp), dimension(60000006) ::xdm,ydm,zdm,mdm,rdm
     !real(kind=xp), dimension(364614)  ::x,y,z,r,m
     !real(kind=xp), dimension(223149095)  ::x,y,z,r,m,ne,T,P,vx,vy,vz,lvl
     real(kind=xp), dimension(:),allocatable  :: x,y,z,r,m,ne,T,P,vx,vy,vz,lvl
     real(kind=xp), dimension(31) :: gal !31
     real(kind=xp), dimension(400):: m_rad,n_rad
     real(kind=xp), dimension(7)::  b
     real(kind=xp), dimension(11):: a
     real(kind=xp), dimension(:), allocatable :: liste,liste2
     real(kind=xp) :: mtot,den,norm,ntot,r500,norm500,n500,m500,r200,norm200,n200,m200,cen_x,cen_y,cen_z,dist_gal
     integer :: i,j,test,ncell,k,index,l,ngal,n,nsub

     test=1




     open(1,file='/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21.dat',form='unformatted')
     !open(1,file='/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l15_low.dat',form='unformatted')


     read(1) ncell
     write(*,*) ncell
     allocate(ne(ncell),T(ncell),P(ncell),x(ncell),y(ncell),z(ncell),vx(ncell),vy(ncell),vz(ncell),m(ncell),lvl(ncell),r(ncell),liste(ncell))

     read(1) ne
     read(1) T
     read(1) P
     read(1) x
     read(1) y
     read(1) z
     read(1) vx
     read(1) vy
     read(1) vz
     read(1) m
     read(1) lvl

     close(1)

     write(*,*) "fin lecture hydro"


     open(3,file='/data/cluster/tlebeau/virgo/HighRes/list_gal_251.dat_js_nocontam_high_res')
     !open(3,file='/data/cluster/tlebeau/virgo/list_gal_251.dat_js_nocontam_low_res')

     cen_x=0.48461068
     cen_y=0.50809848
     cen_z=0.49687076

     cen_x= (cen_x-0.5)*(unit_l/3.08567758128E21)
     cen_y= (cen_y-0.5)*(unit_l/3.08567758128E21)
     cen_z= (cen_z-0.5)*(unit_l/3.08567758128E21)

     i=0
     nsub=0
     do

        read(3,*,end=3) gal
        !if(gal(4)>0.47783 .and. gal(4)<0.49139 .and. gal(5)>0.50131 .and. gal(5)<0.51487 .and. gal(6)>0.49009 .and. gal(6)<0.50365) then
        !if(gal(4)>0.46 .and. gal(4)<0.51 .and. gal(5)>0.48 .and. gal(5)<0.53 .and. gal(6)>0.47 .and. gal(6)<0.52) then
           i=i+1
           xcen(i)=(gal(4)-0.5)*(unit_l/3.08567758128E21)
           ycen(i)=(gal(5)-0.5)*(unit_l/3.08567758128E21)
           zcen(i)=(gal(6)-0.5)*(unit_l/3.08567758128E21)
           rvir(i)=gal(24)*(unit_l/3.08567758128E21)
           mgal(i)=gal(3)
           dist_gal=sqrt((xcen(i)-cen_x)**2+(ycen(i)-cen_y)**2+(zcen(i)-cen_z)**2)
           !if (i<10) then
           ! write(*,*) "dist_gal", dist_gal
           !end if
            if (dist_gal<2113) then
                nsub=nsub+1
           end if
        !endif

     enddo


3    write(*,*) "fin lecture fichiers"


     close(3)

     write(*,*) "nsub", nsub

     write(*,*) "nbr gal", i

     stop

     ngal=i

     write(*,*) "test"

     !open(4,file='virgo_xyz_hydro_l19_gal_clean_m1e8.5.dat',form='unformatted')
     open(4,file='virgo_xyz_hydro_l17_high_gal_clean_m1e8.5.dat',form='unformatted')

     write(*,*) "taille x", size(x)
     write(*,*) "test avant if"

     do i=1,size(x)
        liste(i)=i
     enddo

     index=size(x)

     do k=3,ngal

        if (mgal(k)>10**(8.5)) then

           allocate(liste2(index))

           index=0

           do i=1,size(liste)

              r(liste(i))=sqrt((xcen(k)-x(liste(i)))**2+(ycen(k)-y(liste(i)))**2+(zcen(k)-z(liste(i)))**2)
              if (r(liste(i))>rvir(k)) then
                 index=index+1
                 liste2(index)=liste(i)
              endif
           enddo

           deallocate(liste)
           allocate(liste(index))
           liste(:)=liste2(1:index)
           deallocate(liste2)
        endif

        write(*,*) "gal nbr",k," len liste", size(liste)

     enddo

     write(*,*) "len liste", size(liste)

     write(4) size(liste)
     write(4) ne(liste)
     write(4) T(liste)
     write(4) P(liste)
     write(4) x(liste)
     write(4) y(liste)
     write(4) z(liste)
     write(4) vx(liste)
     write(4) vy(liste)
     write(4) vz(liste)
     write(4) m(liste)
     write(4) lvl(liste)

     close(4)

     deallocate(ne,T,P,x,y,z,vx,vy,vz,m,lvl,liste)

   end function gas_in_gal

    function test_mpi() result(i)
        !include 'mpif.h
        !#include <mpif.h>
        !use mpi

        integer rank, size, ierror,i

        call MPI_INIT(ierror)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)


        write(*,*) 'Hello World from process: ', rank, 'of ', size

        call MPI_FINALIZE(ierror)


    end function test_mpi

    function create_map_3D(proj,pxlvl,unit,weighting,filein,fileout) result(cen_x)


        real(kind=xp), dimension(11)::a
        real(kind=xp), dimension(:,:,:),allocatable :: map,weights
        real(kind=xp) :: cen_x,cen_y,cen_z,minx,maxx,miny,maxy,minz,maxz,lenx,leny,lenz,dx,dy,dz,dist,volume,ratio,i,j,l,start,finish,size_box,x_cen,y_cen,z_cen

        real(kind=xp), dimension(:), allocatable :: x,y,z,P,n_e,T,m,vx,vy,vz,lvl,r
        integer ::k,distx,disty,xinf,xsup,yinf,ysup,ix,iy,test,nx,ny,nz,ncell,n,pxlvl,ierror,distz,zinf,zsup,iz
        integer :: xtest,ytest,ztest,bl,br,tl,tr
        character(len=*)::filein,fileout
        character(len=2)::weighting
        character(len=1)::proj,unit
        real(kind=xp), dimension(21) :: m_rad,n_rad



        write(*,*) "projection: ",proj
        write(*,*)  "unit: ",unit


        open(1,file=filein,form='unformatted')
        read(1) ncell
        write(*,*) "ncell",ncell
        allocate(n_e(ncell),T(ncell),P(ncell),x(ncell),y(ncell),z(ncell),vx(ncell),vy(ncell),vz(ncell),m(ncell),lvl(ncell),r(ncell))
        read(1) n_e
        n_e=n_e/0.864 !conversion to electron density
        write(*,*) "n_e",n_e(1:10)
        read(1) T
        T=T*(kb/1.602e-16) ! conversion to kev instead of kelvin
        read(1) P
        P=P/10
        P=P/(1.602e-10)
        P=P*(0.76/0.864) ! conversion from Pa in cgs system to kev.cm-3 in SI units + conversion to electron pressure (this line)
        read(1) x
        read(1) y
        read(1) z
        read(1) vx
        read(1) vy
        read(1) vz
        read(1) m
        read(1) lvl

        write(*,*) "min lvl", minval(lvl(:))
        !write(*,*) "sum(m)",sum(m)
        close(1)

        !x_cen = 0.48461068
        !y_cen = 0.50809848
        !z_cen = 0.49687076

        !x_cen = (x_cen - 0.5) * (unit_l / 3.08567758128E21)
        !y_cen = (y_cen - 0.5) * (unit_l / 3.08567758128E21)
        !z_cen = (z_cen - 0.5) * (unit_l / 3.08567758128E21)


        !r = log10(sqrt((x - x_cen) ** 2 + (y - y_cen) ** 2 + (z - z_cen) ** 2))
        !m_rad=0
        !n_rad=0

        !do i=1,21
        !   write(*,*) i
  !   do j=1,size(m)
  !      if(r(j)>(i+17.5)*0.1 .and. r(j)<(i+18.5)*0.1) then
  !         m_rad(i)=m_rad(i)+m(j)
  !         n_rad(i)=n_rad(i)+1
  !      endif
  !   enddo
  !enddo

  !open(3,file='/data/cluster/tlebeau/virgo/rad_log_baryons_21.dat',form='unformatted')
  !write(3) m_rad
  !write(3) n_rad

  !write(*,*) m_rad
  !write(*,*) n_rad

  !  close(3)
  !  stop


        cen_x=0.48461068    !Virgo
        cen_y=0.50809848
        cen_z=0.49687076

        cen_x= (cen_x-0.5)*(unit_l/3.08567758128E21)
        cen_y= (cen_y-0.5)*(unit_l/3.08567758128E21)
        cen_z= (cen_z-0.5)*(unit_l/3.08567758128E21)

        write(*,*) "cen_x",cen_x
        write(*,*) "cen_y",cen_y
        write(*,*) "cen_z",cen_z

        !dx=737441.0/2**pxlvl
        !lenx=737441.0*0.03
        size_box=737441 !737433.333
        dx=size_box/2**pxlvl
        lenx=size_box*0.03
        !!Modif lenx pour zoomer sur centre ammas
        !lenx=10000
        nx=int(lenx/dx)
        dy=dx
        dz=dx
        leny=lenx
        lenz=lenx
        ny=nx
        nz=nx
        write(*,*) "dx",dx,"lenx",lenx,"nx",nx
    !stop

        minx=cen_x-lenx/2
        miny=cen_y-leny/2
        minz=cen_z-lenz/2

        write(*,*) "minx",minx
        write(*,*) "miny",miny
        write(*,*) "minz",minz

        !stop

        !!!!Modification for filament visualisation

        !cen_x = -5817
        !cen_y = 441.25
        !cen_z = -2307

        !lenx = (size_box*0.03) / 2

        !nx=int(lenx/dx)
        !dy=dx
        !dz=dx
        !leny=lenx
        !lenz=lenx
        !ny=nx
        !nz=nx

        !minx=cen_x-lenx/2
        !miny=cen_y-leny/2
        !minz=cen_z-lenz/2



        allocate(map(nx,ny,nz))
        allocate(weights(nx,ny,nz))

        map=0
        weights=0

        open(2,file=fileout,form='unformatted')

        write(*,*) "ncell",ncell
        write(*,*) "len m",size(m)

        bl=0
        br=0
        tl=0
        tr=0

      !call MPI_INIT(ierror)

    do k=1,size(m)
      xtest=0
      ytest=0
      ztest=0
      !call cpu_time(start)
     if (proj=="x") then
        i=(y(k)-miny)/dy
        j=(z(k)-minz)/dz
        l=(x(k)-minx)/dx
     endif

     if (proj=="y") then
        i=(x(k)-minx)/dx
        j=(z(k)-minz)/dz
        l=(y(k)-miny)/dy
     endif

     if (proj=="z") then
        i=(x(k)-minx)/dx
        j=(y(k)-miny)/dy
        l=(z(k)-minz)/dz
        !write(*,*) "i real",(x(k)-minx)/dx
        !write(*,*) "j real",(y(k)-miny)/dy
     endif





     if(i>0 .and. j>0 .and. i<nx+1 .and. j<ny+1  .and. l>0 .and. l<nz+1) then !.and. l>int(nx/2)-71 .and. l<int(nx/2)+71 .and. T(k)>8.6142e-3
         !write(*,*) "in"
         distx=2**(pxlvl-lvl(k)-1)
         disty=distx
         distz=distx
         dist=size_box/2**(lvl(k))
         !write(*,*) "test"
         if (int(lvl(k))>=pxlvl) then
             !write(*,*) "yes"
             xinf=int(i)
             xsup=xinf
             yinf=int(j)
             ysup=yinf
             zinf=int(l)
             zsup=zinf
         else
             !if((i-int(i))>0.5) then
             !    xinf=int(i)-distx
             !    xsup=int(i)+(distx-1)
             !    xtest=1
             !else
             !    xinf=int(i)-(distx-1)
             !    xsup=int(i)+distx
             !end if
             !if((j-int(j))>0.5) then
             !    yinf=int(j)-disty
             !    ysup=int(j)+(disty-1)
             !    ztest=xtest+2
             !else
             !    yinf=int(j)-(disty-1)
             !    ysup=int(j)+disty
             !    ztest=xtest
             !end if
             xinf=int(i)-(distx-1)
             xsup=int(i)+distx
             yinf=int(j)-(disty-1)
             ysup=int(j)+disty
             zinf=int(l)-(distz-1)
             zsup=int(l)+distz

         end if
         !if (ztest==0) then
         !    bl=bl+1
         !else if (ztest==1) then
         !    br=br+1
         !else if (ztest==2) then
         !    tl=tl+1
         !else
         !    tr=tr+1
         !end if
         xinf=max(xinf,1)
         xsup=min(xsup,nx)
         yinf=max(yinf,1)
         ysup=min(ysup,ny)
         zinf=max(zinf,1)
         zsup=min(zsup,nz)

         ratio=real(k)/real(size(m))
         if(modulo(k,1000000)==0) then
            write(*,*) "ratio", ratio, "k", k, "lvl",lvl(k)
         end if
         !write(*,*) map(xinf:xsup,yinf:ysup)
         !map(xinf:xsup,yinf:ysup)=map(xinf:xsup,yinf:ysup)+P(k)*m(k)
         !weights(xinf:xsup,yinf:ysup)=weights(xinf:xsup,yinf:ysup)+m(k)

         !write(*,*) "test"
         !weights(ix,iy)=weights(ix,iy)+m(k)
         !malwrite(*,*) map(xinf:xsup,yinf:ysup)
         !weights(ix,iy)=weights(ix,iy)+m(k)
         !stop
         do ix=xinf,xsup
             do iy=yinf,ysup
                 do iz=zinf,zsup
                    if(unit=="T") then
                        if (weighting=="mw") then
                            map(ix,iy,iz)=map(ix,iy,iz)+T(k)*m(k) !n_e(k) !P(k)*m(k)
                            weights(ix,iy,iz)=weights(ix,iy,iz)+m(k)
                        else if (weighting=="sl") then
                            !map(ix,iy)=map(ix,iy)+T(k)*m(k)*n_e(k)*T(k)**(-0.75)
                            !weights(ix,iy)=weights(ix,iy)+m(k)*n_e(k)*T(k)**(-0.75)
                            map(ix,iy,iz)=map(ix,iy,iz)+T(k)*n_e(k)*n_e(k)*T(k)**(-0.75)
                            weights(ix,iy,iz)=weights(ix,iy,iz)+n_e(k)*n_e(k)*T(k)**(-0.75)

                        else if (weighting=="ew") then
                            !map(ix,iy)=map(ix,iy)+T(k)*m(k)*n_e(k)*T(k)**(0.5)
                            !weights(ix,iy)=weights(ix,iy)+m(k)*n_e(k)*T(k)**(0.5)
                            map(ix,iy,iz)=map(ix,iy,iz)+T(k)*n_e(k)*n_e(k)*T(k)**(0.5)
                            weights(ix,iy,iz)=weights(ix,iy,iz)+n_e(k)*n_e(k)*T(k)**(0.5)
                    endif

                    else if(unit=="n") then
                        map(ix,iy,iz)=map(ix,iy,iz)+n_e(k)*dist*pc*1e3
                    else if(unit=="P") then
                        map(ix,iy,iz)=map(ix,iy,iz)+P(k)*m(k)
                    !if (weights(ix,iy)>real(2**31)) then
                        !write(*,*) "map(ix,iy) superior to 2**31"
                        !stop
                    !end if
                        !if (k==42373123) then
                            !write(*,*) "ix",ix,"iy",iy,"m(k)",m(k)
                            !write(*,*) "weights(ix,iy)",weights(ix,iy)
                    !end if
                        weights(ix,iy,iz)=weights(ix,iy,iz)+m(k)
                    else if(unit=="y") then
                        map(ix,iy,iz)=map(ix,iy,iz)+P(k)*dist*pc*1e3*1.602e-10 ! multiplication by pc*1e3 to convert the distance from kpc to m, multiplication by 1.602e-10 to convert the electron pressure from kev.cm-3 to Pa (still in SI units), the compton parameter is dimensionless
                    else if(unit=="v") then
                        map(ix,iy,iz)=map(ix,iy,iz)+vz(k)
                        weights(ix,iy,iz)=weights(ix,iy,iz)+1
                    else if(unit=="vd") then
                        map(ix,iy,iz)=map(ix,iy,iz)+sqrt(vx(k)**2+vy(k)**2)
                        weights(ix,iy,iz)=weights(ix,iy,iz)+1
                    endif
                 enddo

           enddo
        enddo

     endif

      !call cpu_time(finish)
      !print '("Time = ",f6.3," seconds.")',finish-start
      !stop
    enddo

  !call MPI_FINALIZE(ierror)



    if(unit=="P" .or. unit=="T" .or. unit=="v") then
      !write(*,*) "max map", maxval(map(:,:))
      !write(*,*) "max weights", maxval(weights(:,:))
      map(:,:,:)=map(:,:,:)/weights(:,:,:)
    endif

    if(unit=="y") then
      map(:,:,:)=map(:,:,:)*(sigma_t/(me*c**2))
    end if

  !write(*,*) "bl",bl,"br",br,"tl",tl,"tr",tr
  !stop
    write(*,*) "rangement fini"
    write(*,*) "max map", maxval(map)
  !stop

  !do i=1,nx
    write(2) map
  !enddo

  !write(*,*) "map", map(:,:)

  !write(*,*) map

  !write(2,*) nx,ny

    close(2)

    write(*,*) "file saved"

    deallocate(n_e,T,P,x,y,z,vx,vy,vz,m,lvl,map)


    end function create_map_3D

    function create_map_3D_bis(proj,pxlvl,unit,weighting,filein,fileout,datatype,filament) result(cen_x)


        real(kind=xp), dimension(11)::a
        real(kind=xp), dimension(:,:,:),allocatable :: map,weights
        real(kind=xp) :: cen_x,cen_y,cen_z,minx,maxx,miny,maxy,minz,maxz,lenx,leny,lenz,dx,dy,dz,dist,volume,ratio,i,j,l,start,finish,size_box,x_cen,y_cen,z_cen,factor,mean,wmean,vx_virgo_halo, vy_virgo_halo, vz_virgo_halo

        real(kind=xp), dimension(:), allocatable :: x,y,z,P,n_e,T,m,vx,vy,vz,lvl,r
        integer ::k,distx,disty,xinf,xsup,yinf,ysup,ix,iy,test,nx,ny,nz,ncell,n,pxlvl,ierror,ni,nj,nl,zsup,zinf,distz,iz
        integer :: xtest,ytest,ztest,bl,br,tl,tr,filament,in_map,nbr_cells_in_map,nmean,nnan
        character(len=*)::filein,fileout
        character(len=2)::weighting,datatype
        character(len=1)::proj,unit
        real(kind=xp), dimension(21) :: m_rad,n_rad



        write(*,*) "projection: ",proj
        write(*,*)  "unit: ",unit


        open(1,file=filein,form='unformatted')
        read(1) ncell
        write(*,*) "ncell",ncell

        if (datatype=="hy") then

        write(*,*) "datatype:",datatype

        allocate(n_e(ncell),T(ncell),P(ncell),x(ncell),y(ncell),z(ncell),vx(ncell),vy(ncell),vz(ncell),m(ncell),lvl(ncell),r(ncell))
        read(1) n_e
        n_e=n_e/0.864 !conversion to electron density
        write(*,*) "n_e",n_e(1:10)
        read(1) T
        T=T*(kb/1.602e-16) ! conversion to kev instead of kelvin
        read(1) P
        P=P/10
        P=P/(1.602e-10)
        P=P*(0.76/0.864) ! conversion from Pa in cgs system to kev.cm-3 in SI units + conversion to electron pressure (this line)
        read(1) x
        read(1) y
        read(1) z
        read(1) vx
        read(1) vy
        read(1) vz
        read(1) m
        read(1) lvl

        write(*,*) "min lvl", minval(lvl(:))
        !write(*,*) "sum(m)",sum(m)

        write(*,*) "x",x(1:10)

        close(1)

        else if (datatype=="dm") then

        write(*,*) "datatype:",datatype

        allocate(x(ncell),y(ncell),z(ncell),vx(ncell),vy(ncell),vz(ncell),m(ncell),lvl(ncell))
        read(1) x
        read(1) y
        read(1) z
        read(1) vx
        read(1) vy
        read(1) vz
        read(1) m
        lvl = pxlvl

        !write(*,*) "m",m(1:20)

        endif

        !stop


        !x_cen = 0.48461068
        !y_cen = 0.50809848
        !z_cen = 0.49687076

        !x_cen = (x_cen - 0.5) * (unit_l / 3.08567758128E21)
        !y_cen = (y_cen - 0.5) * (unit_l / 3.08567758128E21)
        !z_cen = (z_cen - 0.5) * (unit_l / 3.08567758128E21)


        !r = log10(sqrt((x - x_cen) ** 2 + (y - y_cen) ** 2 + (z - z_cen) ** 2))
        !m_rad=0
        !n_rad=0

        !do i=1,21
        !   write(*,*) i
  !   do j=1,size(m)
  !      if(r(j)>(i+17.5)*0.1 .and. r(j)<(i+18.5)*0.1) then
  !         m_rad(i)=m_rad(i)+m(j)
  !         n_rad(i)=n_rad(i)+1
  !      endif
  !   enddo
  !enddo

  !open(3,file='/data/cluster/tlebeau/virgo/rad_log_baryons_21.dat',form='unformatted')
  !write(3) m_rad
  !write(3) n_rad

  !write(*,*) m_rad
  !write(*,*) n_rad

  !  close(3)
  !  stop


        cen_x=0.48461068    !Virgo
        cen_y=0.50809848
        cen_z=0.49687076

        cen_x= (cen_x-0.5)*(unit_l/3.08567758128E21)
        cen_y= (cen_y-0.5)*(unit_l/3.08567758128E21)
        cen_z= (cen_z-0.5)*(unit_l/3.08567758128E21)

        write(*,*) "cen_x",cen_x, "cen_y",cen_y, "cen_z",cen_z

    !dx=737441.0/2**pxlvl
    !lenx=737441.0*0.03
        size_box=737441 !737433.333
        dx=size_box/2**pxlvl
        lenx=size_box*0.03
        nx=int(lenx/dx)
        dy=dx
        dz=dx
        leny=lenx
        lenz=lenx
        ny=nx
        nz=nx
        write(*,*) "dx",dx,"lenx",lenx,"nx",nx
    !stop

        minx=cen_x-lenx/2
        miny=cen_y-leny/2
        minz=cen_z-lenz/2

        if (filament==0) then

            allocate(map(nx,ny,nz))
            allocate(weights(nx,ny,nz))

        end if

        !!!!Modification for filament visualisation
        if (filament==1) then

            vx_virgo_halo = -507.8579
            vy_virgo_halo = 229.4530
            vz_virgo_halo = -136.9451

            vx = vx - vx_virgo_halo
            vy = vy - vy_virgo_halo
            vz = vz - vz_virgo_halo

            !lenx = lenx/8 ! = 5.53 Mpc
            !leny = leny/8 ! = 5.53 Mpc
            !lenz = lenz/8 !256 ! = 11.06 Mpc

            lenx = 20000
            leny = 10000
            lenz = leny
            !cen_x = cen_x + 6000
            !cen_y = cen_y - 2500

            write(*,*) "lenx",lenx,"leny",leny,"lenz",lenz

            !cen_z = cen_z + 0.5*lenz
            !cen_x = cen_x + 0.75*lenx
            !cen_y = cen_y - 0.5*leny

            !cen_x = cen_x - 0.5*lenx

            write(*,*) "cen_x",cen_x,"cen_y",cen_y,"cen_z",cen_z

            !stop

            minx = cen_x - lenx/2
            miny = cen_y - leny/2
            minz = cen_z - lenz/2

            nx=int(lenx/dx)
            ny=int(leny/dy)
            nz=int(lenz/dz)

            write(*,*) "nx",nx,"ny",ny,"nz",nz

            !stop

            if (proj=="y") then
                allocate(map(nz,nx,nz))
                allocate(weights(nz,nx,nz))
            else if (proj=="z") then
                allocate(map(nx,ny,nz))
                allocate(weights(nx,ny,nz))
            else
                allocate(map(nx,ny,nz))
                allocate(weights(nx,ny,nz))
            end if

        end if


        !cen_x = -5817 !-16878.75
        !cen_y = 441.25
        !cen_z = -2307

        !lenx = (size_box*0.03) / 2

        !nx=int(lenx/dx)
        !dy=dx
        !dz=dx
        !leny=lenx
        !lenz=lenx
        !ny=nx
        !nz=nx

        !minx=cen_x-lenx/2
        !miny=cen_y-leny/2
        !minz=cen_z-lenz/2

        !map=0
        !weights=0

        open(2,file=fileout,form='unformatted')

        write(*,*) "ncell",ncell
        write(*,*) "len m",size(m)

        bl=0
        br=0
        tl=0
        tr=0

      !call MPI_INIT(ierror)
    nbr_cells_in_map = 0

    do k=1,size(m)
      factor = 1
      in_map = 0
      xtest=0
      ytest=0
      ztest=0
      !call cpu_time(start)
     if (proj=="x") then
        i=(y(k)-miny)/dy
        j=(z(k)-minz)/dz
        l=(x(k)-minx)/dx

        ni=ny
        nj=nz
        nl=nx

     endif

     if (proj=="y") then

        i=(z(k)-minz)/dz
        j=(x(k)-minx)/dx
        l=(y(k)-miny)/dy

        ni=nz
        nj=nx
        nl=ny

     endif

     if (proj=="z") then
        i=(x(k)-minx)/dx
        j=(y(k)-miny)/dy
        l=(z(k)-minz)/dz
        !write(*,*) "i real",(x(k)-minx)/dx
        !write(*,*) "j real",(y(k)-miny)/dy

        ni=nx
        nj=ny
        nl=nz

     endif

     distx=2**(pxlvl-lvl(k)-1)
     disty=distx
     distz=distx
     dist=size_box/2**(lvl(k))
     !write(*,*) "test"
     if (int(lvl(k))>=pxlvl) then
             !write(*,*) "yes"
             xinf=int(i)
             xsup=xinf
             yinf=int(j)
             ysup=yinf
             zinf=int(l)
             zsup=zinf

     else
             !if((i-int(i))>0.5) then
             !    xinf=int(i)-distx
             !    xsup=int(i)+(distx-1)
             !    xtest=1
             !else
             !    xinf=int(i)-(distx-1)
             !    xsup=int(i)+distx
             !end if
             !if((j-int(j))>0.5) then
             !    yinf=int(j)-disty
             !    ysup=int(j)+(disty-1)
             !    ztest=xtest+2
             !else
             !    yinf=int(j)-(disty-1)
             !    ysup=int(j)+disty
             !    ztest=xtest
             !end if
             xinf=int(i)-(distx-1)
             xsup=int(i)+distx
             yinf=int(j)-(disty-1)
             ysup=int(j)+disty
             zinf=int(l)-(distz-1)
             zsup=int(l)+dist

     end if





     if(i>0 .and. j>0 .and. i<ni+1 .and. j<nj+1 .and. l>0 .and. l<nl+1) then ! .and. T(k)>8.6142e-1) then !.and. l>int(nx/2)-71 .and. l<int(nx/2)+71

         !if (ztest==0) then
         !    bl=bl+1
         !else if (ztest==1) then
         !    br=br+1
         !else if (ztest==2) then
         !    tl=tl+1
         !else
         !    tr=tr+1
         !end if
         xinf=max(xinf,1)
         xsup=min(xsup,ni)
         yinf=max(yinf,1)
         ysup=min(ysup,nj)

         if (zinf<1) then
            zinf = 1
            factor = (real(zsup)-real(zinf))/(2*real(distz))
         end if

         if (zsup>nl) then
            zsup = nl
            factor = (real(zsup)-real(zinf))/(2*real(distz))
         end if




         !write(*,*) map(xinf:xsup,yinf:ysup)
         !map(xinf:xsup,yinf:ysup)=map(xinf:xsup,yinf:ysup)+P(k)*m(k)
         !weights(xinf:xsup,yinf:ysup)=weights(xinf:xsup,yinf:ysup)+m(k)

         !write(*,*) "test"
         !weights(ix,iy)=weights(ix,iy)+m(k)
         !malwrite(*,*) map(xinf:xsup,yinf:ysup)
         !weights(ix,iy)=weights(ix,iy)+m(k)
         !stop
         in_map = 1

     else
        if (xsup>0 .and. xinf<1) then
           xinf = 1
           !in_map = 1
        end if
        if (xinf<ni+1 .and. xsup>ni) then
           xsup = ni
           !in_map = 1
        end if
        if (ysup>0 .and. yinf<1) then
           yinf = 1
           !in_map = 1
        end if
        if (yinf<nj+1 .and. ysup>nj) then
           ysup = nj
           !in_map = 1
        end if
        if (zsup>0 .and. zinf<1) then
           zinf = 1
           !in_map = 1
           factor = (real(zsup)-real(zinf))/(2*real(distz))
           !write(*,*) "zsup",zsup,"zinf",zinf,"distx",distx,"factor",factor
           !stop
           !write(*,*) "factor",factor
        end if
        if (zinf<nl+1 .and. zsup>nl) then
           zsup = nl
           factor = (real(zsup)-real(zinf))/(2*real(distz))
           !in_map = 1
           !write(*,*) "factor",factor
           !write(*,*) "zsup",zsup,"zinf",zinf,"distx",distx,"factor",factor
           !stop
        end if

        if(xinf>0 .and. yinf>0 .and. xsup<ni+1 .and. ysup<nj+1 .and. zinf>0 .and. zsup<nl+1) then
           in_map = 1
        !   write(*,*) "in map"
           !write(*,*) "zsup",zsup,"zinf",zinf,"distz",distz,"factor",factor,"lvl",lvl(k),"i",i,"j",j,"l",l
           !stop
        !else
        !   in_map = 0
        end if

     end if



      if (in_map==1) then
        nbr_cells_in_map = nbr_cells_in_map + 1

        !write(*,*) "xinf",xinf,"xsup",xsup,"yinf",yinf,"ysup",ysup,"zinf",zinf,"zsup",zsup,"factor",factor

        ratio=real(k)/real(size(m))
         if(modulo(k,1000000)==0) then
            write(*,*) "ratio", ratio, "k", k, "lvl",lvl(k)
            !stop
         end if

        do ix=xinf,xsup
             do iy=yinf,ysup
                 do iz=zinf,zsup
                    if(unit=="T") then
                 if (weighting=="mw") then
                    map(ix,iy,iz)=map(ix,iy,iz)+T(k)*m(k)*factor !n_e(k) !P(k)*m(k)
                    weights(ix,iy,iz)=weights(ix,iy,iz)+m(k)*factor
                 else if (weighting=="sl") then
                    !map(ix,iy)=map(ix,iy)+T(k)*m(k)*n_e(k)*T(k)**(-0.75)
                    !weights(ix,iy)=weights(ix,iy)+m(k)*n_e(k)*T(k)**(-0.75)
                    map(ix,iy,iz)=map(ix,iy,iz)+T(k)*n_e(k)*n_e(k)*T(k)**(-0.75)
                    weights(ix,iy,iz)=weights(ix,iy,iz)+n_e(k)*n_e(k)*T(k)**(-0.75)

                 else if (weighting=="ew") then
                    !map(ix,iy)=map(ix,iy)+T(k)*m(k)*n_e(k)*T(k)**(0.5)
                    !weights(ix,iy)=weights(ix,iy)+m(k)*n_e(k)*T(k)**(0.5)
                    map(ix,iy,iz)=map(ix,iy,iz)+T(k)*n_e(k)*n_e(k)*T(k)**(0.5)
                    weights(ix,iy,iz)=weights(ix,iy,iz)+n_e(k)*n_e(k)*T(k)**(0.5)
                 endif

                    else if(unit=="n") then
                 map(ix,iy,iz)=map(ix,iy,iz)+n_e(k)*dist*pc*1e3*factor
                    else if(unit=="P") then
                 !write(*,*) "in P loop"
                 map(ix,iy,iz)=map(ix,iy,iz)+P(k)*m(k)*factor
                 !if (weights(ix,iy)>real(2**31)) then
                     !write(*,*) "map(ix,iy) superior to 2**31"
                     !stop
                 !end if
                 !if (k==42373123) then
                     !write(*,*) "ix",ix,"iy",iy,"m(k)",m(k)
                     !write(*,*) "weights(ix,iy)",weights(ix,iy)
                 !end if
                 weights(ix,iy,iz)=weights(ix,iy,iz)+m(k)*factor
                    else if(unit=="y") then
                  map(ix,iy,iz)=map(ix,iy,iz)+P(k)*dist*pc*1e3*1.602e-10*factor ! multiplication by pc*1e3 to convert the distance from kpc to m, multiplication by 1.602e-10 to convert the electron pressure from kev.cm-3 to Pa (still in SI units), the compton parameter is dimensionless
                    else if(unit=="v") then
                        if (weighting=="vx") then
                         map(ix,iy,iz)=map(ix,iy,iz)+vx(k)*factor
                     else if (weighting=="vy") then
                         map(ix,iy,iz)=map(ix,iy,iz)+vy(k)*factor
                     else if (weighting=="vz") then
                         !if (T(k)>4) then
                         !    map(ix,iy)=map(ix,iy)+vz(k)*factor
                         !end if
                         map(ix,iy,iz)=map(ix,iy,iz)+vz(k)*factor
                     else if (weighting=="v") then
                         map(ix,iy,iz)=map(ix,iy,iz)+sqrt(vx(k)**2+vy(k)**2+vz(k)**2)*factor
                     end if
                 !map(ix,iy,iz)=map(ix,iy,iz)+sqrt(vx(k)**2+vy(k)**2+vx(k)**2)
                 weights(ix,iy,iz)=weights(ix,iy,iz)+factor
                    else if(unit=="d") then
                 !write(*,*) "ix",ix,"iy",iy,"iz",iz
                 !map(ix,iy,iz)=map(ix,iy,iz)+sqrt(vx(k)**2+vy(k)**2)*factor*m(k)
                 map(ix,iy,iz)=map(ix,iy,iz)+vx(k)*factor*m(k)
                 weights(ix,iy,iz)=weights(ix,iy,iz)+factor*m(k)
                    else if(unit=="m") then ! emission measure
                  !write(*,*) "in EM loop"
                  map(ix,iy,iz)=map(ix,iy,iz)+((n_e(k)**2)/0.864)*(dist*1e-3)*factor ! multiplication by 1e-3 to convert in cm-6*Mpc like Eckert et al. 2012 (before: multiplication by pc*1e3 to convert the distance from kpc to m and *1e2 to convert in cm)
                    !map(ix,iy)=map(ix,iy)+((n_e(k)**2)/0.864)*(dist*pc*1e3*1e2)**3 ! multiplication by 1e-3 to convert in cm-6*Mpc like Eckert et al. 2012 (before: multiplication by pc*1e3 to convert the distance from kpc to m and *1e2 to convert in cm)
                    !write(*,*) "ne**2", (n_e(k)**2)/0.864
                    !write(*,*) "vol", (dist*pc*1e3*1e2)**3
                    !write(*,*) "map",map(ix,iy)
                    !stop
                    else if(unit=="s") then !surface density
                  map(ix,iy,iz)=map(ix,iy,iz)+ (m(k)*factor)/(2**((19-lvl(k))*2))




                    endif
                 enddo
           enddo
        enddo

      end if



      !call cpu_time(finish)
      !print '("Time = ",f6.3," seconds.")',finish-start
      !stop
    enddo

  !call MPI_FINALIZE(ierror)

    write(*,*) 'calculating the mean neighbourhood value for Nan pixels'

      nnan = 0

      do i=1,ni
          do j=1,nj
              do l=1,nl
                if (map(i,j,l)==0) then
                  nnan = nnan + 1
                  !write(*,*) 'pixel is zero'
                  mean = 0
                  nmean = 0
                  wmean = 0
                  do ix=-2,2
                      do iy=-2,2
                          do iz=-2,2
                            if (map(i+ix,j+iy,l+iz)/=0 .and. i+ix>0 .and. j+iy>0 .and. l+iz>0 .and. i+ix<ni .and. j+iy<nj .and. l+iz<nl ) then
                              mean = mean + map(i+ix,j+iy,l+iz)
                              nmean = nmean + 1
                              wmean = wmean + weights(i+ix,j+iy,l+iz)
                          end if
                          end do
                      end do
                  end do
                  map(i,j,l) = mean/nmean
                  weights(i,j,l) = wmean/nmean
                  !write(*,*) 'map(i,j)',map(i,j),'weights(i,j)',weights(i,j)

                end if
              end do
          end do
      end do





    if(unit=="P" .or. unit=="T" .or. unit=="v" .or. unit=="d") then
      !write(*,*) "max map", maxval(map(:,:))
      !write(*,*) "max weights", maxval(weights(:,:))
      map(:,:,:)=map(:,:,:)/weights(:,:,:)
    endif

    if(unit=="y") then
      map(:,:,:)=map(:,:,:)*(sigma_t/(me*c**2))
    end if

    if (unit=="s") then
        map(:,:,:)=map(:,:,:)/(size_box/2**(pxlvl))
    end if

  !write(*,*) "bl",bl,"br",br,"tl",tl,"tr",tr
  !stop
    write(*,*) "rangement fini"
    write(*,*) "max map", maxval(map)
    write(*,*) "nbr_cells_in_map",nbr_cells_in_map
  !stop

  !do i=1,nx
    write(2) nx,ny,nz
    write(2) map
  !enddo

  !write(*,*) "map", map(:,:)

  !write(*,*) map

  !write(2,*) nx,ny

    close(2)

    !if (datatype=="hy") then
    !    deallocate(n_e,T,P,x,y,z,vx,vy,vz,m,lvl,map)
    !else if (datatype=="dm") then
    !    deallocate(x,y,z,vx,vy,vz,m,lvl,map)
    !end if

    !stop


    end function create_map_3D_bis



    !!!!!!!!!!!!  Module to do random projections for splashback paper !!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine create_map_from_random(proj,pxlvl,unit,weighting,fileout,ncell,pos,v,n_e,T,P,m,lvl,datatype)


        real(kind=xp), dimension(11)::a
        real(kind=xp), dimension(:,:),allocatable :: map,weights
        real(kind=xp) :: cen_x,cen_y,cen_z,minx,maxx,miny,maxy,minz,maxz,lenx,leny,lenz,dx,dy,dz,dist,volume,ratio,i,j,l,start,finish,size_box,x_cen,y_cen,z_cen

        !real(kind=xp), dimension(ncell) :: x,y,z,P,n_e,T,m,vx,vy,vz,lvl,r
        real(kind=xp), dimension(ncell) :: x,y,z,vx,vy,vz
        integer ::k,distx,disty,xinf,xsup,yinf,ysup,ix,iy,test,nx,ny,nz,n,ierror
        integer :: xtest,ytest,ztest,bl,br,tl,tr
        real(kind=xp), dimension(21) :: m_rad,n_rad

        character(*), intent(in)::fileout
        character(len=2), intent(in)::weighting,datatype
        character(len=1), intent(in)::proj,unit
        integer, intent(in)::pxlvl,ncell
        real(kind=xp), dimension(ncell), intent(inout) :: P,n_e,T,m,lvl
        real(kind=xp), dimension(3,ncell), intent(in) :: pos,v


        x(:) = pos(1,:)
        y(:) = pos(2,:)
        z(:) = pos(3,:)

        vx(:) = v(1,:)
        vy(:) = v(2,:)
        vz(:) = v(3,:)

        if (datatype=='dm') then
            lvl=pxlvl
        end if



        !write(*,*) "projection: ",proj
        !write(*,*)  "unit: ",unit


        !open(1,file=filein,form='unformatted')
        !read(1) ncell
        !write(*,*) "ncell",ncell

        !if (datatype=="hy") then

        !write(*,*) "datatype:",datatype

        !allocate(n_e(ncell),T(ncell),P(ncell),x(ncell),y(ncell),z(ncell),vx(ncell),vy(ncell),vz(ncell),m(ncell),lvl(ncell),r(ncell))
        !read(1) n_e
        n_e=n_e/0.864 !conversion to electron density
        !write(*,*) "n_e",n_e(1:10)
        !read(1) T
        T=T*(kb/1.602e-16) ! conversion to kev instead of kelvin
        !read(1) P
        P=P/10
        P=P/(1.602e-10)
        P=P*(0.76/0.864) ! conversion from Pa in cgs system to kev.cm-3 in SI units + conversion to electron pressure (this line)
        !read(1) x
        !read(1) y
        !read(1) z
        !read(1) vx
        !read(1) vy
        !read(1) vz
        !read(1) m
        !read(1) lvl

        !write(*,*) "min lvl", minval(lvl(:))
        !write(*,*) "sum(m)",sum(m)

        !write(*,*) "x",x(1:10)

        !close(1)

        !else if (datatype=="dm") then

        !write(*,*) "datatype:",datatype

        !allocate(x(ncell),y(ncell),z(ncell),vx(ncell),vy(ncell),vz(ncell),m(ncell),lvl(ncell))
        !read(1) x
        !read(1) y
        !read(1) z
        !read(1) vx
        !read(1) vy
        !read(1) vz
        !read(1) m
        !lvl = pxlvl

        !write(*,*) "m",m(1:20)

        !endif

        !stop


        !x_cen = 0.48461068
        !y_cen = 0.50809848
        !z_cen = 0.49687076

        !x_cen = (x_cen - 0.5) * (unit_l / 3.08567758128E21)
        !y_cen = (y_cen - 0.5) * (unit_l / 3.08567758128E21)
        !z_cen = (z_cen - 0.5) * (unit_l / 3.08567758128E21)


        !r = log10(sqrt((x - x_cen) ** 2 + (y - y_cen) ** 2 + (z - z_cen) ** 2))
        !m_rad=0
        !n_rad=0

        !do i=1,21
        !   write(*,*) i
  !   do j=1,size(m)
  !      if(r(j)>(i+17.5)*0.1 .and. r(j)<(i+18.5)*0.1) then
  !         m_rad(i)=m_rad(i)+m(j)
  !         n_rad(i)=n_rad(i)+1
  !      endif
  !   enddo
  !enddo

  !open(3,file='/data/cluster/tlebeau/virgo/rad_log_baryons_21.dat',form='unformatted')
  !write(3) m_rad
  !write(3) n_rad

  !write(*,*) m_rad
  !write(*,*) n_rad

  !  close(3)
  !  stop


        cen_x=0.48461068    !Virgo
        cen_y=0.50809848
        cen_z=0.49687076

        cen_x= (cen_x-0.5)*(unit_l/3.08567758128E21)
        cen_y= (cen_y-0.5)*(unit_l/3.08567758128E21)
        cen_z= (cen_z-0.5)*(unit_l/3.08567758128E21)

        !write(*,*) "cen_x",cen_x

    !dx=737441.0/2**pxlvl
    !lenx=737441.0*0.03
        size_box=737441 !737433.333
        dx=size_box/2**pxlvl
        lenx=size_box*0.03
        nx=int(lenx/dx)
        dy=dx
        dz=dx
        leny=lenx
        lenz=lenx
        ny=nx
        nz=nx
        !write(*,*) "dx",dx,"lenx",lenx,"nx",nx
    !stop

        minx=cen_x-lenx/2
        miny=cen_y-leny/2
        minz=cen_z-lenz/2

        allocate(map(nx,ny))
        allocate(weights(nx,ny))

        map=0
        weights=0



        !write(*,*) "ncell",ncell
        !write(*,*) "len m",size(m)

        bl=0
        br=0
        tl=0
        tr=0

      !call MPI_INIT(ierror)

    do k=1,size(m)
      xtest=0
      ytest=0
      ztest=0
      !call cpu_time(start)
     if (proj=="x") then
        i=(y(k)-miny)/dy
        j=(z(k)-minz)/dz
        l=(x(k)-minx)/dx
     endif

     if (proj=="y") then
        i=(x(k)-minx)/dx
        j=(z(k)-minz)/dz
        l=(y(k)-miny)/dy
     endif

     if (proj=="z") then
        i=(x(k)-minx)/dx
        j=(y(k)-miny)/dy
        l=(z(k)-minz)/dz
        !write(*,*) "i real",(x(k)-minx)/dx
        !write(*,*) "j real",(y(k)-miny)/dy
     endif





     if(i>0 .and. j>0 .and. i<nx+1 .and. j<ny+1) then ! .and. T(k)>8.6142e-1) then !.and. l>int(nx/2)-71 .and. l<int(nx/2)+71
         distx=2**(pxlvl-lvl(k)-1)
         disty=distx
         dist=size_box/2**(lvl(k))
         !write(*,*) "test"
         if (int(lvl(k))==pxlvl) then
             !write(*,*) "yes"
             xinf=int(i)
             xsup=xinf
             yinf=int(j)
             ysup=yinf
         else
             !if((i-int(i))>0.5) then
             !    xinf=int(i)-distx
             !    xsup=int(i)+(distx-1)
             !    xtest=1
             !else
             !    xinf=int(i)-(distx-1)
             !    xsup=int(i)+distx
             !end if
             !if((j-int(j))>0.5) then
             !    yinf=int(j)-disty
             !    ysup=int(j)+(disty-1)
             !    ztest=xtest+2
             !else
             !    yinf=int(j)-(disty-1)
             !    ysup=int(j)+disty
             !    ztest=xtest
             !end if
             xinf=int(i)-(distx-1)
             xsup=int(i)+distx
             yinf=int(j)-(disty-1)
             ysup=int(j)+disty

         end if
         !if (ztest==0) then
         !    bl=bl+1
         !else if (ztest==1) then
         !    br=br+1
         !else if (ztest==2) then
         !    tl=tl+1
         !else
         !    tr=tr+1
         !end if
         xinf=max(xinf,1)
         xsup=min(xsup,nx)
         yinf=max(yinf,1)
         ysup=min(ysup,ny)

         ratio=real(k)/real(size(m))
         !if(modulo(k,1000000)==0) then
         !   write(*,*) "ratio", ratio, "k", k, "lvl",lvl(k)
            !stop
         !end if
         !write(*,*) map(xinf:xsup,yinf:ysup)
         !map(xinf:xsup,yinf:ysup)=map(xinf:xsup,yinf:ysup)+P(k)*m(k)
         !weights(xinf:xsup,yinf:ysup)=weights(xinf:xsup,yinf:ysup)+m(k)

         !write(*,*) "test"
         !weights(ix,iy)=weights(ix,iy)+m(k)
         !malwrite(*,*) map(xinf:xsup,yinf:ysup)
         !weights(ix,iy)=weights(ix,iy)+m(k)
         !stop
         do ix=xinf,xsup
             do iy=yinf,ysup
                 if(unit=="T") then
                 if (weighting=="mw") then
                    map(ix,iy)=map(ix,iy)+T(k)*m(k) !n_e(k) !P(k)*m(k)
                    weights(ix,iy)=weights(ix,iy)+m(k)
                 else if (weighting=="sl") then
                    !map(ix,iy)=map(ix,iy)+T(k)*m(k)*n_e(k)*T(k)**(-0.75)
                    !weights(ix,iy)=weights(ix,iy)+m(k)*n_e(k)*T(k)**(-0.75)
                    map(ix,iy)=map(ix,iy)+T(k)*n_e(k)*n_e(k)*T(k)**(-0.75)
                    weights(ix,iy)=weights(ix,iy)+n_e(k)*n_e(k)*T(k)**(-0.75)

                 else if (weighting=="ew") then
                    !map(ix,iy)=map(ix,iy)+T(k)*m(k)*n_e(k)*T(k)**(0.5)
                    !weights(ix,iy)=weights(ix,iy)+m(k)*n_e(k)*T(k)**(0.5)
                    map(ix,iy)=map(ix,iy)+T(k)*n_e(k)*n_e(k)*T(k)**(0.5)
                    weights(ix,iy)=weights(ix,iy)+n_e(k)*n_e(k)*T(k)**(0.5)
                 endif

                else if(unit=="n") then
                 map(ix,iy)=map(ix,iy)+n_e(k)*dist*pc*1e3
                else if(unit=="P") then
                 !write(*,*) "in P loop"
                 map(ix,iy)=map(ix,iy)+P(k)*m(k)
                 !if (weights(ix,iy)>real(2**31)) then
                     !write(*,*) "map(ix,iy) superior to 2**31"
                     !stop
                 !end if
                 !if (k==42373123) then
                     !write(*,*) "ix",ix,"iy",iy,"m(k)",m(k)
                     !write(*,*) "weights(ix,iy)",weights(ix,iy)
                 !end if
                 weights(ix,iy)=weights(ix,iy)+m(k)
                else if(unit=="y") then
                  map(ix,iy)=map(ix,iy)+P(k)*dist*pc*1e3*1.602e-10 ! multiplication by pc*1e3 to convert the distance from kpc to m, multiplication by 1.602e-10 to convert the electron pressure from kev.cm-3 to Pa (still in SI units), the compton parameter is dimensionless
                else if(unit=="v") then
                 map(ix,iy)=map(ix,iy)+vz(k)
                 weights(ix,iy)=weights(ix,iy)+1
                else if(unit=="vd") then
                 map(ix,iy)=map(ix,iy)+sqrt(vx(k)**2+vy(k)**2)
                 weights(ix,iy)=weights(ix,iy)+1
                else if(unit=="m") then ! emission measure
                  !write(*,*) "in EM loop"
                  map(ix,iy)=map(ix,iy)+((n_e(k)**2)/0.864)*(dist*1e-3) ! multiplication by 1e-3 to convert in cm-6*Mpc like Eckert et al. 2012 (before: multiplication by pc*1e3 to convert the distance from kpc to m and *1e2 to convert in cm)
                  !map(ix,iy)=map(ix,iy)+((n_e(k)**2)/0.864)*(dist*pc*1e3*1e2)**3 ! multiplication by 1e-3 to convert in cm-6*Mpc like Eckert et al. 2012 (before: multiplication by pc*1e3 to convert the distance from kpc to m and *1e2 to convert in cm)
                  !write(*,*) "ne**2", (n_e(k)**2)/0.864
                  !write(*,*) "vol", (dist*pc*1e3*1e2)**3
                  !write(*,*) "map",map(ix,iy)
                  !stop
                else if(unit=="s") then !surface density
                  map(ix,iy)=map(ix,iy)+ m(k)/(2**((19-lvl(k))*2))




                endif

           enddo
        enddo

     endif

      !call cpu_time(finish)
      !print '("Time = ",f6.3," seconds.")',finish-start
      !stop
    enddo

  !call MPI_FINALIZE(ierror)



    if(unit=="P" .or. unit=="T" .or. unit=="v") then
      !write(*,*) "max map", maxval(map(:,:))
      !write(*,*) "max weights", maxval(weights(:,:))
      map(:,:)=map(:,:)/weights(:,:)
    endif

    if(unit=="y") then
      map(:,:)=map(:,:)*(sigma_t/(me*c**2))
    end if

    if (unit=="s") then
        map(:,:)=map(:,:)/(size_box/2**(pxlvl))
    end if

  !write(*,*) "bl",bl,"br",br,"tl",tl,"tr",tr
  !stop
    !write(*,*) "rangement fini"
    !write(*,*) "max map", maxval(map)
  !stop

  !do i=1,nx
    open(2,file=fileout,form='unformatted')

    write(2) map
  !enddo

  !write(*,*) "map", map(:,:)

  !write(*,*) map

  !write(2,*) nx,ny

    close(2)

    write(*,*) "file saved"

    !if (datatype=="hy") then
    !    deallocate(n_e,T,P,x,y,z,vx,vy,vz,m,lvl,map)
    !else if (datatype=="dm") then
    !    deallocate(x,y,z,vx,vy,vz,m,lvl,map)
    !end if


    end subroutine create_map_from_random



    function random_rotations(n_i,n_f,datatype) result(phi)
      real(kind=xp), dimension(:),allocatable:: n_e,T,P,m,lvl,random_phi,random_theta,n_e_rand,T_rand,P_rand,m_rand,lvl_rand
      real(kind=xp), dimension(:,:),allocatable::pos,v,pos_rand,v_rand
      real(kind=xp) :: mw_x,mw_y,mw_z,phi,theta,s,c,zero,un,fil_x,fil_y,fil_z,vx_virgo_halo,vy_virgo_halo,vz_virgo_halo
      real(kind=xp), dimension(3,3):: rx,ry,rz
      integer:: ncell,i,file_len,n_i,n_f
      real(kind=xp), dimension(3):: mw,fil,rotaxis,virgo
      character(len=:), allocatable ::fileout
      character(len=2) :: i_str,datatype

      write(*,*) "open datafile"

      if (datatype=="hy") then
        open(1,file='virgo_xyz_hydro_l19_gal_clean_m1e8.5.dat',form='unformatted')
        !open(1,file='virgo_xyz_hydro_l19.dat',form='unformatted')

        read(1) ncell
        allocate(n_e(ncell),T(ncell),P(ncell),m(ncell),lvl(ncell),pos(3,ncell),v(3,ncell),pos_rand(3,ncell),v_rand(3,ncell),n_e_rand(ncell),T_rand(ncell),P_rand(ncell),m_rand(ncell),lvl_rand(ncell))

        read(1) n_e
        !ne=ne*(1E6/0.864)
        read(1) T
        read(1) P
        !P=P/10
        read(1) pos(1,:)
        read(1) pos(2,:)
        read(1) pos(3,:)
        read(1) v(1,:)
        read(1) v(2,:)
        read(1) v(3,:)
        read(1) m
        read(1) lvl

        close(1)

      else if (datatype=="dm") then
        open(1,file='virgo_xyz_dm_high_res.dat',form='unformatted')

        read(1) ncell

        allocate(m(ncell),pos(3,ncell),v(3,ncell),pos_rand(3,ncell),v_rand(3,ncell),n_e_rand(ncell),T_rand(ncell),P_rand(ncell),m_rand(ncell),lvl_rand(ncell))

        !allocate(m(ncell),pos(3,ncell),v(3,ncell),pos_rand(3,ncell),v_rand(3,ncell),m_rand(ncell))

        read(1) pos(1,:)
        read(1) pos(2,:)
        read(1) pos(3,:)
        read(1) v(1,:)
        read(1) v(2,:)
        read(1) v(3,:)
        read(1) m

        close(1)
      end if


      !stop




      virgo=(/0.48461068,0.50809848,0.49687076/)
      fil=(/0.497,0.4984,0.4996/)
      mw=(/0.5,0.5,0.5/)

      !pos(:,1)=(/0,0,0/)
      pos(:,1)=(fil(:)-mw(:))*(unit_l/3.08567758128E21)

      !write(*,*) "pos 1",pos(:,1)

      rotaxis=mw


      !transformation in Virgo frame of reference (i.e. Virgo at the origin)

      vx_virgo_halo = -507.8579
      vy_virgo_halo = 229.4530
      vz_virgo_halo = -136.9451

      v(1,:) = v(1,:) - vx_virgo_halo
      v(2,:) = v(2,:) - vy_virgo_halo
      v(3,:) = v(3,:) - vz_virgo_halo


      pos(1,:)=pos(1,:)-(virgo(1)-mw(1))*(unit_l/3.08567758128E21)
      pos(2,:)=pos(2,:)-(virgo(2)-mw(2))*(unit_l/3.08567758128E21)
      pos(3,:)=pos(3,:)-(virgo(3)-mw(3))*(unit_l/3.08567758128E21)


      allocate(random_phi(100),random_theta(100))


      !do i = 1,100
      !  random_phi(i) = 2*pi*ran()
      !  random_theta(i) = pi*ran()
      !end do

      open(2,file='./maps/high_res/random_proj/random_theta.dat',form='unformatted')
      read(2) random_theta
      close(2)

      open(3,file='./maps/high_res/random_proj/random_phi.dat',form='unformatted')
      read(3) random_phi
      close(3)

      !write(*,*) "random angles saved"

      !write(*,*) "random_phi",random_phi
      !write(*,*) "random_theta",random_theta

      !stop

      !write(*,*) "i",i

      !write(i_str,'(I1)') i

      !write(*,*) "i_str",i_str

      !file_len = len('./maps/high_res/random_proj/map_high_19_rand_') + len('_y_los.bin') + 1

      !allocate(character(file_len) :: fileout)

      !fileout = trim('./maps/high_res/random_proj/map_high_19_rand_')//trim(i_str)//trim('_y_los.bin')

      !write(*,*) "fileout: ",fileout

      !stop

      zero=0.0
      un=1.0

      do i = n_i,n_f
        write(*,*) "i",i

        if (datatype=="hy") then

          n_e_rand=n_e
          T_rand=T
          P_rand=P
          m_rand=m
          lvl_rand=lvl

        else if (datatype=="dm") then

            m_rand=m
            lvl_rand=0
            n_e_rand=0
            T_rand=0
            P_rand=0


        end if


          if (i==0) then
              c=cos(random_phi(100))
              s=sin(random_phi(100))

          else
              c=cos(random_phi(i))
              s=sin(random_phi(i))

          end if

          !c=cos(random_phi(i))
          !s=sin(random_phi(i))

          rz(:,1)=(/c,-s,zero/)
          rz(:,2)=(/s,c,zero/)
          rz(:,3)=(/zero,zero,un/)

          pos_rand=transpose(matmul(transpose(pos),rz))
          v_rand=transpose(matmul(transpose(v),rz))

          if (i==0) then
              c=cos(random_theta(100))
              s=sin(random_theta(100))

          else
              c=cos(random_theta(i))
              s=sin(random_theta(i))

          end if

          rx(:,1)=(/un,zero,zero/)
          rx(:,2)=(/zero,c,-s/)
          rx(:,3)=(/zero,s,c/)

          pos_rand=transpose(matmul(transpose(pos_rand),rx))
          v_rand=transpose(matmul(transpose(v_rand),rx))

          pos_rand(1,:)=pos_rand(1,:)+(virgo(1)-mw(1))*(unit_l/3.08567758128E21)
          pos_rand(2,:)=pos_rand(2,:)+(virgo(2)-mw(2))*(unit_l/3.08567758128E21)
          pos_rand(3,:)=pos_rand(3,:)+(virgo(3)-mw(3))*(unit_l/3.08567758128E21)



          file_len = len('./maps/high_res/random_proj/em/map_high_19_rand_') + len('_em_los.bin') + 2
          allocate(character(file_len) :: fileout)

          if (i<10) then
              write(i_str,'(I1)') i
              fileout = trim('./maps/high_res/random_proj/em/map_high_19_rand_0')//trim(i_str)//trim('_em_los.bin')
          else
              write(i_str,'(I2)') i
              fileout = trim('./maps/high_res/random_proj/em/map_high_19_rand_')//trim(i_str)//trim('_em_los.bin')
          end if

          write(*,*) "fileout: ",fileout

          !stop

          call create_map_from_random("z",19,"m","mw",fileout,ncell,pos_rand,v_rand,n_e_rand,T_rand,P_rand,m_rand,lvl_rand,'hy')

          deallocate(fileout)


          !a=create_map_ter("z",19,"y","mw",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l19_MW_los.dat','./maps/high_res/map_high_19_cen_y_los.bin','hy')



      end do



      !pos(:,1)=(/0,0,0/)


      !mw_x=0.5-0.48461068
      !mw_y=0.5-0.50809848
      !mw_z=0.5-0.49687076

      !fil_x=0.5-0.497
      !fil_y=0.5-0.4984
      !fil_z=0.5-0.4996





      !phi=atan(mw_x/mw_z)
      !phi=atan(mw_y/mw_z)
      !theta=atan(mw_x/mw_y)
      !theta=atan(fil_x/fil_y)

      !theta=atan((rotaxis(1)-virgo(1))/(rotaxis(2)-virgo(2)))

      !write(*,*) "phi",phi
      !write(*,*) "theta",theta

      !theta=pi/2
      !phi=0

      !write(*,*) "mw",pos(:,1)

      !c=cos(theta)
      !s=sin(theta)

      !write(*,*) "cos",c,"sin",s

      !rz(:,1)=(/c,-s,zero/)
      !rz(:,2)=(/s,c,zero/)
      !rz(:,3)=(/zero,zero,un/)
      !write(*,*) "rz",rz

      !pos=transpose(matmul(transpose(pos),rz))
      !v=transpose(matmul(transpose(v),rz))

      !write(*,*) "mw",pos(:,1)

      !phi=atan(pos(2,1)/pos(3,1))
      !phi=

      !write(*,*) "phi",phi

      !c=cos(phi)
      !s=sin(phi)



      !ry(:,1)=(/c,zero,s/)
      !ry(:,2)=(/zero,un,zero/)
      !ry(:,3)=(/-s,zero,c/)

      !rx(:,1)=(/un,zero,zero/)
      !rx(:,2)=(/zero,c,-s/)
      !rx(:,3)=(/zero,s,c/)

      !pos=transpose(matmul(transpose(pos),rx))
      !v=transpose(matmul(transpose(v),rx))

      !phi=pi/2

      !write(*,*) "phi",phi

      !c=cos(phi)
      !s=sin(phi)

      !write(*,*) "cos",c,"sin",s



      !ry(:,1)=(/c,zero,s/)
      !ry(:,2)=(/zero,un,zero/)
      !ry(:,3)=(/-s,zero,c/)

      !rx(:,1)=(/un,zero,zero/)
      !rx(:,2)=(/zero,c,-s/)
      !rx(:,3)=(/zero,s,c/)

      !pos=transpose(matmul(transpose(pos),ry))

      !write(*,*) "pos 1",pos(:,1)


      !saving data in new file

      !pos(1,:)=pos(1,:)+(virgo(1)-mw(1))*(unit_l/3.08567758128E21)
      !pos(2,:)=pos(2,:)+(virgo(2)-mw(2))*(unit_l/3.08567758128E21)
      !pos(3,:)=pos(3,:)+(virgo(3)-mw(3))*(unit_l/3.08567758128E21)

      !write(*,*) "pos 1",pos(:,1)

      !stop

      !write(*,*) "ecriture"

      !open(2,file='virgo_xyz_hydro_l19_all_bar_MW_los.dat',form='unformatted')
      !write(2) ncell
      !write(2) ne
      !write(2) T
      !write(2) P
      !write(2) pos(1,:)
      !write(2) pos(2,:)
      !write(2) pos(3,:)
      !write(2) v(1,:)
      !write(2) v(2,:)
      !write(2) v(3,:)
      !write(2) m
      !write(2) lvl
      !close(2)

      !write(*,*) "fin ecriture"

    end function random_rotations


    function stack_maps() result(size_box)
        real(kind=xp), dimension(:,:),allocatable :: map,map_sum
        integer :: nx,i,file_len
        real(kind=xp) :: lenx,dx,size_box,pxlvl
        character(len=:), allocatable ::fileout,filein
        character(len=2) :: i_str

        pxlvl = 19

        size_box=737441 !737433.333
        dx=size_box/2**pxlvl
        lenx=size_box*0.03
        nx=int(lenx/dx)

        allocate(map(nx,nx),map_sum(nx,nx))

        do i=0,99

            write(*,*) "i",i

            !write(i_str,'(I2)') i

            !file_len = len('./maps/high_res/random_proj/map_high_19_rand_') + len('_y_los.bin') + 2

            !allocate(character(file_len) :: filein)

            !fileout = trim('./maps/high_res/random_proj/map_high_19_rand_')//trim(i_str)//trim('_y_los.bin')

            if (i<10) then
                write(*,*) "i<10"

                write(i_str,'(I1)') i
                file_len = len('./maps/high_res/random_proj/sd_bar/map_high_19_rand_') + len('_sd_bar_los.bin') + 2
                allocate(character(file_len) :: filein)
                filein = trim('./maps/high_res/random_proj/sd_bar/map_high_19_rand_0')//trim(i_str)//trim('_sd_bar_los.bin')

            else
                write(*,*) "i>=10"

                write(i_str,'(I2)') i
                file_len = len('./maps/high_res/random_proj/sd_bar/map_high_19_rand_') + len('_sd_bar_los.bin') + 2
                allocate(character(file_len) :: filein)
                filein = trim('./maps/high_res/random_proj/sd_bar/map_high_19_rand_')//trim(i_str)//trim('_sd_bar_los.bin')

            end if

            write(*,*) "filein : ",filein


            open(2,file=filein,form='unformatted')

            read(2) map

            close(2)

            deallocate(filein)

            map_sum=map_sum+map

        end do

        !stop

        map_sum=map_sum/100

        fileout = './maps/high_res/random_proj/sd_bar/map_high_19_rand_sum100_sd_bar_los.bin'

        open(3,file=fileout,form='unformatted')
        write(3) map_sum
        close(3)

        write(*,*) "file saved"

    end function stack_maps


    !!!!!! project on VSF in Virgo core w S.ettori

    function VSF(filein,order,test) result(nx)

        character(len=*)::filein
        real(kind=xp), dimension(:,:),allocatable :: map
        integer :: nx,ny,nz,i,j,k,test,order



        open(1,file=filein,form='unformatted')

        read(1) nx,ny,nz
        write(*,*) "nx,ny,nz",nx,ny,nz
        read(1) map

        close(1)

        write(*,*) "reading done"


    end function VSF

    subroutine create_map_from_random_bis(proj,pxlvl,unit,weighting,fileout,ncell,pos,v,n_e,T,P,m,lvl,datatype,map_width_in_mpc,tcut)


        real(kind=xp), dimension(11)::a
        real(kind=xp), dimension(:,:),allocatable :: map,weights,map2
        real(kind=xp) :: cen_x,cen_y,cen_z,minx,maxx,miny,maxy,minz,maxz,lenx,leny,lenz,dx,dy,dz,dist,volume,ratio,i,j,l,start,finish,size_box,x_cen,y_cen,z_cen,factor
        real(kind=xp), intent(in) :: tcut
        !real(kind=xp), dimension(ncell) :: x,y,z,P,n_e,T,m,vx,vy,vz,lvl,r
        real(kind=xp), dimension(ncell) :: x,y,z,vx,vy,vz,rho
        integer ::k,distx,disty,distz,xinf,xsup,yinf,ysup,ix,iy,test,nx,ny,nz,n,ierror,ni,nj,nl,zinf,zsup
        integer :: xtest,ytest,ztest,bl,br,tl,tr,in_map,nbr_cells_in_map
        real(kind=xp), dimension(21) :: m_rad,n_rad

        character(*), intent(in)::fileout
        character(len=2), intent(in)::datatype
        character(len=4), intent(in)::weighting
        character(len=1), intent(in)::proj,unit
        integer, intent(in)::pxlvl,ncell,map_width_in_mpc
        real(kind=xp), dimension(ncell), intent(inout) :: P,n_e,T,m,lvl
        real(kind=xp), dimension(3,ncell), intent(in) :: pos,v


        x(:) = pos(1,:)
        y(:) = pos(2,:)
        z(:) = pos(3,:)

        vx(:) = v(1,:)
        vy(:) = v(2,:)
        vz(:) = v(3,:)

        !write(*,*) 'ncell in',ncell
        nbr_cells_in_map=0
        !stop

        if (datatype=='dm') then
            lvl=pxlvl
        end if



        !write(*,*) "projection: ",proj
        !write(*,*)  "unit: ",unit


        !open(1,file=filein,form='unformatted')
        !read(1) ncell
        !write(*,*) "ncell",ncell

        !if (datatype=="hy") then

        !write(*,*) "datatype:",datatype

        !allocate(n_e(ncell),T(ncell),P(ncell),x(ncell),y(ncell),z(ncell),vx(ncell),vy(ncell),vz(ncell),m(ncell),lvl(ncell),r(ncell))
        !read(1) n_e
        n_e=n_e/0.864 !conversion to electron density
        !write(*,*) "n_e",n_e(1:10)
        !read(1) T
        T=T*(kb/1.602e-16) ! conversion to kev instead of kelvin
        !read(1) P
        P=P/10
        P=P/(1.602e-10)
        P=P*(0.76/0.864) ! conversion from Pa in cgs system to kev.cm-3 in SI units + conversion to electron pressure (this line)
        !read(1) x
        !read(1) y
        !read(1) z
        !read(1) vx
        !read(1) vy
        !read(1) vz
        !read(1) m
        !read(1) lvl

        !write(*,*) "min lvl", minval(lvl(:))
        !write(*,*) "sum(m)",sum(m)

        !write(*,*) "x",x(1:10)

        !close(1)

        !else if (datatype=="dm") then

        !write(*,*) "datatype:",datatype

        !allocate(x(ncell),y(ncell),z(ncell),vx(ncell),vy(ncell),vz(ncell),m(ncell),lvl(ncell))
        !read(1) x
        !read(1) y
        !read(1) z
        !read(1) vx
        !read(1) vy
        !read(1) vz
        !read(1) m
        !lvl = pxlvl

        !write(*,*) "m",m(1:20)

        !endif

        !stop


        !x_cen = 0.48461068
        !y_cen = 0.50809848
        !z_cen = 0.49687076

        !x_cen = (x_cen - 0.5) * (unit_l / 3.08567758128E21)
        !y_cen = (y_cen - 0.5) * (unit_l / 3.08567758128E21)
        !z_cen = (z_cen - 0.5) * (unit_l / 3.08567758128E21)


        !r = log10(sqrt((x - x_cen) ** 2 + (y - y_cen) ** 2 + (z - z_cen) ** 2))
        !m_rad=0
        !n_rad=0

        !do i=1,21
        !   write(*,*) i
  !   do j=1,size(m)
  !      if(r(j)>(i+17.5)*0.1 .and. r(j)<(i+18.5)*0.1) then
  !         m_rad(i)=m_rad(i)+m(j)
  !         n_rad(i)=n_rad(i)+1
  !      endif
  !   enddo
  !enddo

  !open(3,file='/data/cluster/tlebeau/virgo/rad_log_baryons_21.dat',form='unformatted')
  !write(3) m_rad
  !write(3) n_rad

  !write(*,*) m_rad
  !write(*,*) n_rad

  !  close(3)
  !  stop


        cen_x=0.48461068    !Virgo
        cen_y=0.50809848
        cen_z=0.49687076

        cen_x= (cen_x-0.5)*(unit_l/3.08567758128E21)
        cen_y= (cen_y-0.5)*(unit_l/3.08567758128E21)
        cen_z= (cen_z-0.5)*(unit_l/3.08567758128E21)

        !write(*,*) "cen_x",cen_x

    !dx=737441.0/2**pxlvl
    !lenx=737441.0*0.03
        size_box=737441 !737433.333
        dx=size_box/2**pxlvl
        lenx=size_box*0.03
        nx=int(lenx/dx)
        dy=dx
        dz=dx
        leny=lenx
        lenz=lenx
        ny=nx
        nz=nx
        !write(*,*) "dx",dx,"lenx",lenx,"nx",nx
    !stop



        if (proj=="x") then
            !lenx = lenx/221.23
            leny = leny/(22.123/map_width_in_mpc)
            lenz = lenz/(22.123/map_width_in_mpc)
        else if (proj=="y") then
            lenx = lenx/(22.123/map_width_in_mpc)
            !leny = leny/221.23
            lenz = lenz/(22.123/map_width_in_mpc)
        else if (proj=="z") then
            lenx = lenx/(22.123/map_width_in_mpc)
            leny = leny/(22.123/map_width_in_mpc)
            !lenz = lenz/2212.3
        end if

        minx=cen_x-lenx/2
        miny=cen_y-leny/2
        minz=cen_z-lenz/2

        nx=int(lenx/dx)
        ny=int(leny/dy)
        nz=int(lenz/dz)

        if (proj=="y") then
                allocate(map(nz,nx))
                allocate(weights(nz,nx))
                if (unit=="d") then
                    allocate(map2(nz,nx))
                end if
        else if (proj=="z") then
                allocate(map(nx,ny))
                allocate(weights(nx,ny))
                if (unit=="d") then
                    allocate(map2(nx,ny))
                end if
        else
                allocate(map(nz,ny))
                allocate(weights(nz,ny))
                if (unit=="d") then
                    allocate(map2(nz,ny))
                end if
        end if

        map=0
        weights=0
        map2=0

        !write(*,*) "map", map(1:10,1:10)
        !write(*,*) "weights", weights(1:10,1:10)

        !write(*,*) "lenx",lenx,"leny",leny,"lenz",lenz

        !write(*,*) "cen_x",cen_x,"cen_y",cen_y,"cen_z",cen_z

        !write(*,*) "nx",nx,"ny",ny,"nz",nz

        !map=0
        !weights=0



        !write(*,*) "ncell",ncell
        !write(*,*) "len m",size(m)

        bl=0
        br=0
        tl=0
        tr=0

      !call MPI_INIT(ierror)

    do k=1,size(m)
        factor = 1
      in_map = 0
      xtest=0
      ytest=0
      ztest=0
      !call cpu_time(start)
     if (proj=="x") then
        i=(y(k)-miny)/dy
        j=(z(k)-minz)/dz
        l=(x(k)-minx)/dx

        ni=nz
        nj=ny
        nl=nx
     endif

     if (proj=="y") then
        i=(x(k)-minx)/dx
        j=(z(k)-minz)/dz
        l=(y(k)-miny)/dy

        ni=nz
        nj=nx
        nl=ny

     endif

     if (proj=="z") then
        i=(x(k)-minx)/dx
        j=(y(k)-miny)/dy
        l=(z(k)-minz)/dz
        !write(*,*) "i real",(x(k)-minx)/dx
        !write(*,*) "j real",(y(k)-miny)/dy

        ni=nx
        nj=ny
        nl=nz
     endif

     distx=2**(pxlvl-lvl(k)-1)
     disty=distx
     distz=distx
     dist=size_box/2**(lvl(k))
     !write(*,*) "test"
     if (int(lvl(k))>=pxlvl) then
             !write(*,*) "yes"
             xinf=int(i)
             xsup=xinf
             yinf=int(j)
             ysup=yinf
             zinf=int(l)
             zsup=zinf

     else
             !if((i-int(i))>0.5) then
             !    xinf=int(i)-distx
             !    xsup=int(i)+(distx-1)
             !    xtest=1
             !else
             !    xinf=int(i)-(distx-1)
             !    xsup=int(i)+distx
             !end if
             !if((j-int(j))>0.5) then
             !    yinf=int(j)-disty
             !    ysup=int(j)+(disty-1)
             !    ztest=xtest+2
             !else
             !    yinf=int(j)-(disty-1)
             !    ysup=int(j)+disty
             !    ztest=xtest
             !end if
             xinf=int(i)-(distx-1)
             xsup=int(i)+distx
             yinf=int(j)-(disty-1)
             ysup=int(j)+disty
             zinf=int(l)-(distz-1)
             zsup=int(l)+distz

     end if





         if(i>0 .and. j>0 .and. i<ni+1 .and. j<nj+1 .and. l>0 .and. l<nl+1) then !.and. l>int(nx/2)-71 .and. l<int(nx/2)+71

         !if (ztest==0) then
         !    bl=bl+1
         !else if (ztest==1) then
         !    br=br+1
         !else if (ztest==2) then
         !    tl=tl+1
         !else
         !    tr=tr+1
         !end if
         xinf=max(xinf,1)
         xsup=min(xsup,ni)
         yinf=max(yinf,1)
         ysup=min(ysup,nj)

         if (zinf<0) then
            zinf = 1
            factor = (real(zsup)-real(zinf))/(2*real(distz))
         end if

         if (zsup>nl) then
            zsup = nl
            factor = (real(zsup)-real(zinf))/(2*real(distz))
         end if

         !zinf=max(zinf,1)
         !zsup=min(zsup,nl)
         !factor = (real(zsup)-real(zinf))/(2*real(distz))

         !write(*,*) map(xinf:xsup,yinf:ysup)
         !map(xinf:xsup,yinf:ysup)=map(xinf:xsup,yinf:ysup)+P(k)*m(k)
         !weights(xinf:xsup,yinf:ysup)=weights(xinf:xsup,yinf:ysup)+m(k)

         !write(*,*) "test"
         !weights(ix,iy)=weights(ix,iy)+m(k)
         !malwrite(*,*) map(xinf:xsup,yinf:ysup)
         !weights(ix,iy)=weights(ix,iy)+m(k)
         !stop
         in_map = 1
     else
        if (xsup>0 .and. xinf<1) then
           xinf = 1
           !in_map = 1
        end if
        if (xinf<ni+1 .and. xsup>ni) then
           xsup = ni
           !in_map = 1
        end if
        if (ysup>0 .and. yinf<1) then
           yinf = 1
           !in_map = 1
        end if
        if (yinf<nj+1 .and. ysup>nj) then
           ysup = nj
           !in_map = 1
        end if
        if (zsup>0 .and. zinf<1) then
           zinf = 1
           !in_map = 1
           factor = (real(zsup)-real(zinf))/(2*real(distz))
           !write(*,*) "zsup",zsup,"zinf",zinf,"distx",distx,"factor",factor
           !stop
           !write(*,*) "factor",factor
        end if
        if (zinf<nl+1 .and. zsup>nl) then
           zsup = nl
           factor = (real(zsup)-real(zinf))/(2*real(distz))
           !in_map = 1
           !write(*,*) "factor",factor
           !write(*,*) "zsup",zsup,"zinf",zinf,"distx",distx,"factor",factor
           !stop
        end if

        if(xinf>0 .and. yinf>0 .and. xsup<ni+1 .and. ysup<nj+1 .and. zinf>0 .and. zsup<nl+1) then
           in_map = 1
        !   write(*,*) "in map"
           !write(*,*) "zsup",zsup,"zinf",zinf,"distz",distz,"factor",factor,"lvl",lvl(k),"i",i,"j",j,"l",l
           !stop
        !else
        !   in_map = 0
        end if

     end if

      if (in_map==1 .and. T(k)>tcut) then !.and. T(k)>1.7 .and. T(k)<12) then !
        !write(*,*) "T",T(k)
        !stop

        nbr_cells_in_map = nbr_cells_in_map + 1

        !write(*,*) "xinf",xinf,"xsup",xsup,"yinf",yinf,"ysup",ysup,"zinf",zinf,"zsup",zsup,"factor",factor

        ratio=real(k)/real(size(m))
         if(modulo(k,1000000)==0) then
            write(*,*) "ratio", ratio, "k", k, "lvl",lvl(k)
            !write(*,*) "max map", maxval(map), "min map", minval(map)
            !stop
         end if

        do ix=xinf,xsup
             do iy=yinf,ysup
                 if(unit=="T") then
                    if (weighting=="mw") then
                    map(ix,iy)=map(ix,iy)+T(k)*m(k)*factor !n_e(k) !P(k)*m(k)
                    weights(ix,iy)=weights(ix,iy)+m(k)*factor
                    !if (ix>27710 .and. ix<27720 .and. iy>29385 .and. iy<29395) then
                    !    write(*,*) "i",i,"j",j,"l",l,"xinf",xinf,"xsup",xsup,"yinf",yinf,"ysup",ysup,"factor",factor,'lvl',lvl(k)
                    !end if

                    else if (weighting=="sl") then
                    !map(ix,iy)=map(ix,iy)+T(k)*m(k)*n_e(k)*T(k)**(-0.75)
                    !weights(ix,iy)=weights(ix,iy)+m(k)*n_e(k)*T(k)**(-0.75)
                    map(ix,iy)=map(ix,iy)+T(k)*n_e(k)*n_e(k)*T(k)**(-0.75)*factor
                    weights(ix,iy)=weights(ix,iy)+n_e(k)*n_e(k)*T(k)**(-0.75)*factor

                    else if (weighting=="ew") then
                    !map(ix,iy)=map(ix,iy)+T(k)*m(k)*n_e(k)*T(k)**(0.5)
                    !weights(ix,iy)=weights(ix,iy)+m(k)*n_e(k)*T(k)**(0.5)
                    map(ix,iy)=map(ix,iy)+T(k)*n_e(k)*n_e(k)*T(k)**(0.5)
                    weights(ix,iy)=weights(ix,iy)+n_e(k)*n_e(k)*T(k)**(0.5)
                    endif

                 else if(unit=="n") then
                 map(ix,iy)=map(ix,iy)+n_e(k)*dist*pc*1e3*factor*1e2
                 !map(ix,iy)=map(ix,iy)+n_e(k)*m(k)*factor
                 !test for filament project maps
                 !map(ix,iy)=map(ix,iy)+n_e(k)*factor !n_e(k) !P(k)*m(k)
                 !weights(ix,iy)=weights(ix,iy)+factor*m(k)
                 weights(ix,iy)=weights(ix,iy)+factor

                else if(unit=="P") then
                 !write(*,*) "in P loop"
                 map(ix,iy)=map(ix,iy)+P(k)*m(k)*factor
                 !if (weights(ix,iy)>real(2**31)) then
                     !write(*,*) "map(ix,iy) superior to 2**31"
                     !stop
                 !end if
                 !if (k==42373123) then
                     !write(*,*) "ix",ix,"iy",iy,"m(k)",m(k)
                     !write(*,*) "weights(ix,iy)",weights(ix,iy)
                 !end if
                 weights(ix,iy)=weights(ix,iy)+m(k)*factor
                else if(unit=="y") then
                  map(ix,iy)=map(ix,iy)+P(k)*dist*pc*1e3*1.602e-10*factor ! multiplication by pc*1e3 to convert the distance from kpc to m, multiplication by 1.602e-10 to convert the electron pressure from kev.cm-3 to Pa (still in SI units), the compton parameter is dimensionless
                else if(unit=="v") then
                     if (weighting=="vxew") then
                         map(ix,iy)=map(ix,iy)+vx(k)*factor*n_e(k)**2
                         weights(ix,iy)=weights(ix,iy)+factor*n_e(k)**2
                         !map(ix,iy)=map(ix,iy)+vx(k)*factor*m(k)
                     else if (weighting=="vyew") then
                         map(ix,iy)=map(ix,iy)+vy(k)*factor*n_e(k)**2
                         weights(ix,iy)=weights(ix,iy)+factor*n_e(k)**2
                     else if (weighting=="vzew") then
                         map(ix,iy)=map(ix,iy)+vz(k)*factor*n_e(k)**2
                         weights(ix,iy)=weights(ix,iy)+factor*n_e(k)**2
                         !write(*,*) "test" , vz(k)*factor*n_e(k)**2
                         !stop
                     else if (weighting=="vew") then
                         map(ix,iy)=map(ix,iy)+sqrt(vx(k)**2+vy(k)**2+vz(k)**2)*factor*n_e(k)**2
                         weights(ix,iy)=weights(ix,iy)+factor*n_e(k)**2

                     else if (weighting=="vxmw") then
                         map(ix,iy)=map(ix,iy)+vx(k)*factor*m(k)
                         weights(ix,iy)=weights(ix,iy)+factor*m(k)
                     else if (weighting=="vymw") then
                         map(ix,iy)=map(ix,iy)+vy(k)*factor*m(k)
                         weights(ix,iy)=weights(ix,iy)+factor*m(k)
                     else if (weighting=="vzmw") then
                         map(ix,iy)=map(ix,iy)+vz(k)*factor*m(k)
                         weights(ix,iy)=weights(ix,iy)+factor*m(k)
                     else if (weighting=="vmw") then
                         map(ix,iy)=map(ix,iy)+sqrt(vx(k)**2+vy(k)**2+vz(k)**2)*factor*m(k)
                         weights(ix,iy)=weights(ix,iy)+factor*m(k)
                     end if
                 !map(ix,iy)=map(ix,iy)+vx(k)*factor
                     !if (weighting=="vz") then
                     !    weights(ix,iy)=weights(ix,iy)+factor*n_e(k)**2
                     !else
                     !    weights(ix,iy)=weights(ix,iy)+factor
                     !end if

                 !weights(ix,iy)=weights(ix,iy)+factor*m(k)
                else if(unit=="d") then !velocity dispersion
                     !write(*,*) "test"
                     !stop
                    if (weighting=="vxew") then
                        map(ix,iy)=map(ix,iy)+vx(k)**2*factor*n_e(k)**2
                        map2(ix,iy)=map2(ix,iy)+vx(k)*factor*n_e(k)**2
                        weights(ix,iy)=weights(ix,iy)+factor*n_e(k)**2
                    else if (weighting=="vyew") then
                        map(ix,iy)=map(ix,iy)+vy(k)**2*factor*n_e(k)**2
                        map2(ix,iy)=map2(ix,iy)+vy(k)*factor*n_e(k)**2
                        weights(ix,iy)=weights(ix,iy)+factor*n_e(k)**2
                    else if (weighting=="vzew") then
                        map(ix,iy)=map(ix,iy)+vz(k)**2*factor*n_e(k)**2
                        map2(ix,iy)=map2(ix,iy)+vz(k)*factor*n_e(k)**2
                        weights(ix,iy)=weights(ix,iy)+factor*n_e(k)**2
                    else if (weighting=="vew") then
                        map(ix,iy)=map(ix,iy)+sqrt(vx(k)**2+vy(k)**2+vz(k)**2)**2*factor*n_e(k)**2
                        map2(ix,iy)=map2(ix,iy)+sqrt(vx(k)**2+vy(k)**2+vz(k)**2)*factor*n_e(k)**2
                        weights(ix,iy)=weights(ix,iy)+factor*n_e(k)**2

                    else if (weighting=="vxmw") then
                        map(ix,iy)=map(ix,iy)+vx(k)**2*factor*m(k)
                        map2(ix,iy)=map2(ix,iy)+vx(k)*factor*m(k)
                        weights(ix,iy)=weights(ix,iy)+factor*m(k)
                    else if (weighting=="vymw") then
                        !write(*,*) "test"
                        !stop
                        map(ix,iy)=map(ix,iy)+vy(k)**2*factor*m(k)
                        map2(ix,iy)=map2(ix,iy)+vy(k)*factor*m(k)
                        weights(ix,iy)=weights(ix,iy)+factor*m(k)
                    else if (weighting=="vzmw") then
                        map(ix,iy)=map(ix,iy)+vz(k)**2*factor*m(k)
                        map2(ix,iy)=map2(ix,iy)+vz(k)*factor*m(k)
                        weights(ix,iy)=weights(ix,iy)+factor*m(k)
                    else if (weighting=="vmw") then
                        map(ix,iy)=map(ix,iy)+sqrt(vx(k)**2+vy(k)**2+vz(k)**2)**2*factor*m(k)
                        map2(ix,iy)=map2(ix,iy)+sqrt(vx(k)**2+vy(k)**2+vz(k)**2)*factor*m(k)
                        weights(ix,iy)=weights(ix,iy)+factor*m(k)
                    end if



                else if(unit=="m") then ! emission measure
                  !write(*,*) "in EM loop"
                  map(ix,iy)=map(ix,iy)+((n_e(k)**2)/0.864)*(dist*1e-3)*factor ! multiplication by 1e-3 to convert in cm-6*Mpc like Eckert et al. 2012 (before: multiplication by pc*1e3 to convert the distance from kpc to m and *1e2 to convert in cm)
                  !map(ix,iy)=map(ix,iy)+((n_e(k)**2)/0.864)*(dist*pc*1e3*1e2)**3 ! multiplication by 1e-3 to convert in cm-6*Mpc like Eckert et al. 2012 (before: multiplication by pc*1e3 to convert the distance from kpc to m and *1e2 to convert in cm)
                  !write(*,*) "ne**2", (n_e(k)**2)/0.864
                  !write(*,*) "vol", (dist*pc*1e3*1e2)**3
                  !write(*,*) "map",map(ix,iy)
                  !stop
                else if(unit=="s") then !surface density
                  map(ix,iy)=map(ix,iy)+ (m(k)*factor)/(2**((19-lvl(k))*2))

                else if(unit=="k") then
                     map(ix,iy)=map(ix,iy)+(T(k)/n_e(k)**(2/3))*m(k)*factor !n_e(k) !P(k)*m(k)
                     weights(ix,iy)=weights(ix,iy)+m(k)*factor

                else if(unit=="f") then
                     if (weighting=="vx") then
                         map(ix,iy)=map(ix,iy)+vx(k)*rho(k)*factor
                     else if (weighting=="vy") then
                         map(ix,iy)=map(ix,iy)+vy(k)*rho(k)*factor
                     else if (weighting=="vz") then
                         map(ix,iy)=map(ix,iy)+vz(k)*rho(k)*factor

                         end if
                     weights(ix,iy)=weights(ix,iy)+factor

                else if(unit=="a") then  ! compute mach number map
                     map(ix,iy)=map(ix,iy)+(sqrt(vx(k)**2+vy(k)**2+vz(k)**2)/sqrt(((5 / 3) * kb * T(k)) / (mu * mp)))*factor
                     weights(ix,iy)=weights(ix,iy)+factor



                endif

           enddo
        enddo
      end if


    end do

      !write(*,*) 'number of nan pixels',nnan

      !stop



      !call cpu_time(finish)
      !print '("Time = ",f6.3," seconds.")',finish-start
      !stop


  !call MPI_FINALIZE(ierror)

    !write(*,*) "map",map(27710:27720,29391)
    !write(*,*) "map", map(1:10,1:10)
    !write(*,*) "weights", weights(1:10,1:10)


    if(unit=="P" .or. unit=="T" .or. unit=="v" .or. unit=="k" .or. unit=="f" .or. unit=="a") then
      !write(*,*) "max map", maxval(map(:,:))
      !write(*,*) "max weights", maxval(weights(:,:))
      !write(*,*) "max map", maxval(map),'min map',minval(map)
      map(:,:)=map(:,:)/weights(:,:)
      !write(*,*) "max map", maxval(map),'min map',minval(map)

    !if(unit=="n") then !test for filament project maps
    !  map(:,:)=map(:,:)/weights(:,:)
    !end if
    endif

    if(unit=="y") then
      map(:,:)=map(:,:)*(sigma_t/(me*c**2))
    end if

    if (unit=="s") then
        map(:,:)=map(:,:)/(size_box/2**(pxlvl))
    end if

    if (unit=="d") then
        map2(:,:)=map2(:,:)/weights(:,:)
        map(:,:)=map(:,:)/weights(:,:) - map2(:,:)**2
        map(:,:)=sqrt(map(:,:))
    end if

    if (unit=="n") then
        map(:,:)=map(:,:)/weights(:,:)
    end if

  !write(*,*) "bl",bl,"br",br,"tl",tl,"tr",tr
  !stop
    !write(*,*) "rangement fini"
    !write(*,*) "max map", maxval(map)
  !stop

  !do i=1,nx
    write(*,*) "rangement fini"
    write(*,*) "max map", maxval(map), "min map", minval(map)
    write(*,*) "nbr_cells_in_map",nbr_cells_in_map
  !stop
    open(2,file=fileout,form='unformatted')
    write(2) nx,ny,nz
    write(2) cen_x,cen_y,cen_z
    write(2) map
  !enddo

  !write(*,*) "map", map(:,:)

  !write(*,*) map

  !write(2,*) nx,ny

    close(2)

    write(*,*) "file saved"

    !if (datatype=="hy") then
    !    deallocate(n_e,T,P,x,y,z,vx,vy,vz,m,lvl,map)
    !else if (datatype=="dm") then
    !    deallocate(x,y,z,vx,vy,vz,m,lvl,map)
    !end if

    deallocate(map,weights)


    end subroutine create_map_from_random_bis

    function random_rotations_bis(n_i,n_f,datatype) result(phi)
      real(kind=xp), dimension(:),allocatable:: n_e,T,P,m,lvl,random_phi,random_theta,n_e_rand,T_rand,P_rand,m_rand,lvl_rand
      real(kind=xp), dimension(:,:),allocatable::pos,v,pos_rand,v_rand
      real(kind=xp) :: mw_x,mw_y,mw_z,phi,theta,s,c,zero,un,fil_x,fil_y,fil_z,vx_virgo_halo,vy_virgo_halo,vz_virgo_halo,tcut
      real(kind=xp), dimension(3,3):: rx,ry,rz
      integer:: ncell,i,file_len,n_i,n_f
      real(kind=xp), dimension(3):: mw,fil,rotaxis,virgo
      character(len=:), allocatable ::fileout
      character(len=2) :: i_str,datatype

      write(*,*) "open datafile"

      if (datatype=="hy") then
        open(1,file='./virgo_xyz_files/virgo_xyz_hydro_l15_high.dat',form='unformatted')
        !open(1,file='virgo_xyz_hydro_l19.dat',form='unformatted')

        read(1) ncell
        write(*,*) "ncell",ncell
        allocate(n_e(ncell),T(ncell),P(ncell),m(ncell),lvl(ncell),pos(3,ncell),v(3,ncell),pos_rand(3,ncell),v_rand(3,ncell),n_e_rand(ncell),T_rand(ncell),P_rand(ncell),m_rand(ncell),lvl_rand(ncell))

        read(1) n_e
        !ne=ne*(1E6/0.864)
        read(1) T
        read(1) P
        !P=P/10
        read(1) pos(1,:)
        read(1) pos(2,:)
        read(1) pos(3,:)
        read(1) v(1,:)
        read(1) v(2,:)
        read(1) v(3,:)
        read(1) m
        read(1) lvl

        close(1)

      else if (datatype=="dm") then
        open(1,file='virgo_xyz_dm_high_res.dat',form='unformatted')

        read(1) ncell

        allocate(m(ncell),pos(3,ncell),v(3,ncell),pos_rand(3,ncell),v_rand(3,ncell),n_e_rand(ncell),T_rand(ncell),P_rand(ncell),m_rand(ncell),lvl_rand(ncell))

        !allocate(m(ncell),pos(3,ncell),v(3,ncell),pos_rand(3,ncell),v_rand(3,ncell),m_rand(ncell))

        read(1) pos(1,:)
        read(1) pos(2,:)
        read(1) pos(3,:)
        read(1) v(1,:)
        read(1) v(2,:)
        read(1) v(3,:)
        read(1) m

        close(1)
      end if


      !stop




      virgo=(/0.48461068,0.50809848,0.49687076/)
      fil=(/0.497,0.4984,0.4996/)
      mw=(/0.5,0.5,0.5/)

      !pos(:,1)=(/0,0,0/)
      pos(:,1)=(fil(:)-mw(:))*(unit_l/3.08567758128E21)

      !write(*,*) "pos 1",pos(:,1)

      rotaxis=mw


      !transformation in Virgo frame of reference (i.e. Virgo at the origin)

      vx_virgo_halo = -507.8579
      vy_virgo_halo = 229.4530
      vz_virgo_halo = -136.9451

      v(1,:) = v(1,:) - vx_virgo_halo
      v(2,:) = v(2,:) - vy_virgo_halo
      v(3,:) = v(3,:) - vz_virgo_halo


      pos(1,:)=pos(1,:)-(virgo(1)-mw(1))*(unit_l/3.08567758128E21)
      pos(2,:)=pos(2,:)-(virgo(2)-mw(2))*(unit_l/3.08567758128E21)
      pos(3,:)=pos(3,:)-(virgo(3)-mw(3))*(unit_l/3.08567758128E21)


      allocate(random_phi(100),random_theta(100))


      !do i = 1,100
      !  random_phi(i) = 2*pi*ran()
      !  random_theta(i) = pi*ran()
      !end do

      !write(*,*) "random_phi",random_phi
      !write(*,*) "random_theta",random_theta

      !stop

      !open(2,file='./maps/high_res/random_proj/random_theta.dat',form='unformatted')
      open(2,file='./maps/high_res/velocity/15f15_analysis/random_proj/random_theta.dat',form='unformatted')
      read(2) random_theta
      close(2)

      open(3,file='./maps/high_res/velocity/15f15_analysis/random_proj/random_phi.dat',form='unformatted')
      read(3) random_phi
      close(3)

      !write(*,*) "random angles saved"

      !write(*,*) "random_phi",random_phi
      !write(*,*) "random_theta",random_theta

      !write(*,*) "angles saved"

      !stop

      !write(*,*) "i",i

      !write(i_str,'(I1)') i

      !write(*,*) "i_str",i_str

      !file_len = len('./maps/high_res/random_proj/map_high_19_rand_') + len('_y_los.bin') + 1

      !allocate(character(file_len) :: fileout)

      !fileout = trim('./maps/high_res/random_proj/map_high_19_rand_')//trim(i_str)//trim('_y_los.bin')

      !write(*,*) "fileout: ",fileout

      !stop

      zero=0.0
      un=1.0

      do i = n_i,n_f
        write(*,*) "i",i

        if (datatype=="hy") then

          n_e_rand=n_e
          T_rand=T
          P_rand=P
          m_rand=m
          lvl_rand=lvl
          pos_rand=pos
          v_rand=v

        else if (datatype=="dm") then

            m_rand=m
            lvl_rand=0
            n_e_rand=0
            T_rand=0
            P_rand=0


        end if


          if (i==0) then
              c=cos(random_phi(100))
              s=sin(random_phi(100))

          else
              c=cos(random_phi(i))
              s=sin(random_phi(i))

          end if

          !c=cos(random_phi(i))
          !s=sin(random_phi(i))

          rz(:,1)=(/c,-s,zero/)
          rz(:,2)=(/s,c,zero/)
          rz(:,3)=(/zero,zero,un/)

          pos_rand=transpose(matmul(transpose(pos_rand),rz))
          v_rand=transpose(matmul(transpose(v_rand),rz))

          if (i==0) then
              c=cos(random_theta(100))
              s=sin(random_theta(100))

          else
              c=cos(random_theta(i))
              s=sin(random_theta(i))

          end if

          rx(:,1)=(/un,zero,zero/)
          rx(:,2)=(/zero,c,-s/)
          rx(:,3)=(/zero,s,c/)

          pos_rand=transpose(matmul(transpose(pos_rand),rx))
          v_rand=transpose(matmul(transpose(v_rand),rx))

          pos_rand(1,:)=pos_rand(1,:)+(virgo(1)-mw(1))*(unit_l/3.08567758128E21)
          pos_rand(2,:)=pos_rand(2,:)+(virgo(2)-mw(2))*(unit_l/3.08567758128E21)
          pos_rand(3,:)=pos_rand(3,:)+(virgo(3)-mw(3))*(unit_l/3.08567758128E21)



          !file_len = len('./maps/high_res/random_proj/em/map_high_19_rand_') + len('_em_los.bin') + 2
          file_len = len('./maps/high_res/velocity/15f15_analysis/random_proj/map_high_15f15_rand_') + 2 + len('_map_vd_ew_Tsup7_5Mpc2.bin')

          !!ew!!
          !allocate(character(file_len) :: fileout)

          !if (i<10) then
          !    write(i_str,'(I1)') i
          !    fileout = trim('./maps/high_res/velocity/15f15_analysis/random_proj/map_high_15f15_rand_0')//trim(i_str)//trim('_map_v_ew_Tsup7_5Mpc2.bin')

          !else
          !    write(i_str,'(I2)') i
              !fileout = trim('./maps/high_res/random_proj/em/map_high_19_rand_')//trim(i_str)//trim('_em_los.bin')
          !    fileout = trim('./maps/high_res/velocity/15f15_analysis/random_proj/map_high_15f15_rand_')//trim(i_str)//trim('_map_v_ew_Tsup7_5Mpc2.bin')
          !end if

          !write(*,*) "fileout: ",fileout

          !stop

          !tcut = 8.6142e-1

          !call create_map_from_random_bis("z",15,"v","vzew",fileout,ncell,pos_rand,v_rand,n_e_rand,T_rand,P_rand,m_rand,lvl_rand,'hy',5,tcut)

          !!mw!!
          !deallocate(fileout)
          allocate(character(file_len) :: fileout)

          if (i<10) then
              write(i_str,'(I1)') i
              fileout = trim('./maps/high_res/velocity/15f15_analysis/random_proj/map_high_15f15_rand_0')//trim(i_str)//trim('_map_vd_ew_Tsup7_5Mpc2.bin')

          else
              write(i_str,'(I2)') i
              !fileout = trim('./maps/high_res/random_proj/em/map_high_19_rand_')//trim(i_str)//trim('_em_los.bin')
              fileout = trim('./maps/high_res/velocity/15f15_analysis/random_proj/map_high_15f15_rand_')//trim(i_str)//trim('_map_vd_ew_Tsup7_5Mpc2.bin')
          end if

          write(*,*) "fileout: ",fileout

          !stop

          tcut = 8.6142e-1

          !write(*,*) "ncell out 2",ncell


          call create_map_from_random_bis("z",15,"d","vzew",fileout,ncell,pos_rand,v_rand,n_e_rand,T_rand,P_rand,m_rand,lvl_rand,'hy',5,tcut)


          deallocate(fileout)


          !a=create_map_ter("z",19,"y","mw",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l19_MW_los.dat','./maps/high_res/map_high_19_cen_y_los.bin','hy')



      end do



      !pos(:,1)=(/0,0,0/)


      !mw_x=0.5-0.48461068
      !mw_y=0.5-0.50809848
      !mw_z=0.5-0.49687076

      !fil_x=0.5-0.497
      !fil_y=0.5-0.4984
      !fil_z=0.5-0.4996





      !phi=atan(mw_x/mw_z)
      !phi=atan(mw_y/mw_z)
      !theta=atan(mw_x/mw_y)
      !theta=atan(fil_x/fil_y)

      !theta=atan((rotaxis(1)-virgo(1))/(rotaxis(2)-virgo(2)))

      !write(*,*) "phi",phi
      !write(*,*) "theta",theta

      !theta=pi/2
      !phi=0

      !write(*,*) "mw",pos(:,1)

      !c=cos(theta)
      !s=sin(theta)

      !write(*,*) "cos",c,"sin",s

      !rz(:,1)=(/c,-s,zero/)
      !rz(:,2)=(/s,c,zero/)
      !rz(:,3)=(/zero,zero,un/)
      !write(*,*) "rz",rz

      !pos=transpose(matmul(transpose(pos),rz))
      !v=transpose(matmul(transpose(v),rz))

      !write(*,*) "mw",pos(:,1)

      !phi=atan(pos(2,1)/pos(3,1))
      !phi=

      !write(*,*) "phi",phi

      !c=cos(phi)
      !s=sin(phi)



      !ry(:,1)=(/c,zero,s/)
      !ry(:,2)=(/zero,un,zero/)
      !ry(:,3)=(/-s,zero,c/)

      !rx(:,1)=(/un,zero,zero/)
      !rx(:,2)=(/zero,c,-s/)
      !rx(:,3)=(/zero,s,c/)

      !pos=transpose(matmul(transpose(pos),rx))
      !v=transpose(matmul(transpose(v),rx))

      !phi=pi/2

      !write(*,*) "phi",phi

      !c=cos(phi)
      !s=sin(phi)

      !write(*,*) "cos",c,"sin",s



      !ry(:,1)=(/c,zero,s/)
      !ry(:,2)=(/zero,un,zero/)
      !ry(:,3)=(/-s,zero,c/)

      !rx(:,1)=(/un,zero,zero/)
      !rx(:,2)=(/zero,c,-s/)
      !rx(:,3)=(/zero,s,c/)

      !pos=transpose(matmul(transpose(pos),ry))

      !write(*,*) "pos 1",pos(:,1)


      !saving data in new file

      !pos(1,:)=pos(1,:)+(virgo(1)-mw(1))*(unit_l/3.08567758128E21)
      !pos(2,:)=pos(2,:)+(virgo(2)-mw(2))*(unit_l/3.08567758128E21)
      !pos(3,:)=pos(3,:)+(virgo(3)-mw(3))*(unit_l/3.08567758128E21)

      !write(*,*) "pos 1",pos(:,1)

      !stop

      !write(*,*) "ecriture"

      !open(2,file='virgo_xyz_hydro_l19_all_bar_MW_los.dat',form='unformatted')
      !write(2) ncell
      !write(2) ne
      !write(2) T
      !write(2) P
      !write(2) pos(1,:)
      !write(2) pos(2,:)
      !write(2) pos(3,:)
      !write(2) v(1,:)
      !write(2) v(2,:)
      !write(2) v(3,:)
      !write(2) m
      !write(2) lvl
      !close(2)

      !write(*,*) "fin ecriture"

    end function random_rotations_bis



    

      
  end module func
