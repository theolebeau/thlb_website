program mainf90

  use const
  use func
  implicit none

  integer:: i,ncell,j,namas,file_len
  real(kind=xp) a,x_cen,y_cen,z_cen,start,finish,tcut
  real(kind=xp), dimension(:), allocatable :: x,y,z,P,n_e,T,m,vx,vy,vz,lvl,r
  real(kind=xp), dimension(21) :: m_rad,n_rad
  real(kind=xp), dimension(3,3) :: mass_tensor
  character(len=128)::arg
  character(len=:), allocatable:: titre
  character(len=2):: i_str
  !real(kind=xp), dimension(2) :: atest,btest

  !atest(:)=(/1,4/)
  !btest(:)=(/2,3/)

  !mass_tensor=0

  !write(*,*) 'a',atest,'b',btest,'a+b',atest+btest,'a*b', atest(:)*btest(:)

  !write(*,*) "test", 1e2, 1d2

  !stop


  !write(*,*) "mass tensor" , mass_tensor

  !stop

  !call getarg(1,arg)
  !read(arg,*) namas

  !write(*,*) "namas",namas

  !stop


  i=0
  !call cpu_time(start)
  !i=test_mpi()
  !stop
  !i=gas_in_gal()
  !stop
  !i=cut_lines_in_file('./maps/high_res/map_high_21_fil_P_los.bin','./maps/high_res/map_high_21_fil_P_los_nline.bin')
  !i=cut_lines_in_file('./maps/high_res/P_maps/map_high_19_fil_P_los.bin','./maps/high_res/P_maps/map_high_19_fil_P_los_nline.bin')
  !a=rotation()
  !stop
  !i=ascii_to_bin()
  !i=remove_gal()

  !stop

  !i=rescale_data("virgo_xyz_hydro_l15_high.dat","virgo_xyz_hydro_l15_high_rescale.dat",15,0)

  !i=rescale_data("./virgo_xyz_files/virgo_xyz_hydro_l15_high.dat","virgo_xyz_hydro_l15_high_rescale.dat",15,"all")

  !a=stack_maps()

  !a=random_rotations_bis(1,100,'hy')
  !stop


  !a=create_map_ter("z",19,"s","mw",'/data/cluster/tlebeau/virgo/virgo_xyz_dm_high_res_MW_los.dat','./maps/high_res/map_high_19_cen_SD_DM_los.bin','dm')
  !a=create_map_ter("z",19,"s","mw",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l19_all_bar_MW_los.dat','./maps/high_res/map_high_19_cen_SD_bar_los.bin','hy')
  !a=create_map_ter("z",19,"P","mw",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l19.dat','./maps/high_res/map_high_19_z_P_filzoom.bin','hy')


  !a=create_map_ter("y",21,"T","mw",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21_all_bar_fil_los.dat','./maps/high_res/filament/map_high_21_fil_map_T_y_1.38d.bin','hy',1)

  !a=create_map_3D_bis("z",18,"v","mw",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21.dat','./maps/high_res/filament/3D_fields/map_3D_high_18_xyz_core_map_v_z.bin','hy',1)

  !a=create_map_ter("x",19,"n","mw",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21.dat','./maps/high_res/filament/map_high_19_xyz_left_minus_8Mpc_map_ne_x_d0.02.bin','hy',1,0)
  !a=create_map_ter("x",19,"f","vz",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21.dat','./maps/high_res/filament/map_high_19_xyz_left_minus_8Mpc_map_fz_x_d0.02.bin','hy',1,0)
  a=0
  a=create_map_ter("x",19,"l","vx",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l21.dat','./maps/high_res/filament/map_high_19_xyz_left_minus_2Mpc_map_lvl_x_d0.02.bin','hy',1,22,a)
  !a=create_map_ter("z",19,"a","vz",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21.dat','./maps/high_res/filament/map_high_19_xyz_full_map_mach_z_d0.02.bin','hy',1,0)


  !a=create_map_ter("x",19,"T","mw",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21.dat','./maps/high_res/filament/map_high_19_xyz_left_minus_0Mpc_map_T_x_d0.02.bin','hy',1,0)
  !a=create_map_ter("x",19,"v","vz",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21.dat','./maps/high_res/filament/map_high_19_xyz_left_plus_2Mpc_map_vz_x_d0.02.bin','hy',1,0)
  !a=create_map_ter("x",19,"v","vy",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21.dat','./maps/high_res/filament/map_high_19_xyz_left_plus_2Mpc_map_vy_x_d0.02.bin','hy',1,0)

  stop

  !a=create_map_ter("x",19,"n","mw",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21.dat','./maps/high_res/velocity/map_high_19_x_core_map_ne.bin','hy',0,0)
  !a=create_map_ter("y",19,"n","mw",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21.dat','./maps/high_res/velocity/map_high_19_y_core_map_ne.bin','hy',0,0)
  !a=create_map_ter("z",17,"v","vz",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21.dat','./maps/high_res/velocity/map_high_17_z_map_vz.bin','hy',0,0)
  !a=create_map_ter("z",17,"v","vz",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l17_high.dat','./maps/high_res/velocity/map_high_17f17_z_map_vz_ew.bin','hy',0,0)
  !a=create_map_ter("z",21,"d","vz",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21.dat','./maps/high_res/velocity/map_high_21_z_map_vdz_ew.bin','hy',0,0)

  !a=create_map_ter("z",16,"v","vz",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l16_high.dat','./maps/high_res/velocity/map_high_16f16_z_map_vz_mw_5Mpc2.bin','hy',0,0)
  !a=create_map_ter("y",16,"v","vy",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l16_high.dat','./maps/high_res/velocity/map_high_16f16_y_map_vy_mw_5Mpc2.bin','hy',0,0)
  !a=create_map_ter("x",16,"v","vx",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l16_high.dat','./maps/high_res/velocity/map_high_16f16_x_map_vx_mw_5Mpc2.bin','hy',0,0)
  !a=create_map_ter("x",16,"d","vx",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l16_high.dat','./maps/high_res/velocity/map_high_16f16_x_map_vdx_ew.bin','hy',0,0)
  !a=create_map_ter("y",16,"v","vy",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l16_high.dat','./maps/high_res/velocity/map_high_16f16_analysis/map_high_16f16_y_map_vy_ew_Tsup7.bin','hy',0,0)

  !a=create_map_ter("z",16,"d","vzmw",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l16_high_MW_los.dat','./maps/high_res/velocity/map_high_16f16_analysis/map_high_16f16_cen_map_vdz_mw_Tsup7_5Mpc2.bin','hy',0,0)
  !a=create_map_ter("z",16,"d","vzmw",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l16_high.dat','./maps/high_res/velocity/map_high_16f16_analysis/map_high_16f16_z_map_vdz_mw_Tsup7_5Mpc2.bin','hy',0,0)
  !a=create_map_ter("z",21,"v","vz",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21.dat','./maps/high_res/velocity/map_high_21_z_pm_1.3_Mpc_map_vz_d0.02.bin','hy',1,0)

  !a=create_map_ter("z",21,"v","vz",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21_all_bar_MW_los.dat','./maps/high_res/velocity/map_high_21_z_core170kpc_4kev_cut_map_vz.bin','hy',0,0)

  !a=create_map_ter("z",15,"v","v",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l15_high_rescale_vel.dat','./maps/high_res/filament/map_high_15_xyz_full_map_v_z_d0.02.bin','ve',1,0)

  !virgo_xyz_hydro_l21_all_bar_MW_los.dat
  !'./maps/high_res/filament/T_slices/map_high_19_xyz_full_map_T_x_d0.08.bin'

  !i=VSF("./maps/high_res/velocity/map_high_19_z_core_map_vz.bin",1,1)

  !a=create_map_ter("z",15,"v","vzmw",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l15_high.dat','./maps/high_res/velocity/15f15_analysis/map_high_15f15_x_map_vx_mw_Tsup7_5Mpc2.bin','hy',0,0)

  !a=create_map_ter("z",15,"v","vzew",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l15_MW_los.dat','./maps/high_res/velocity/15f15_analysis/map_high_15f15_cen_map_vz_ew_Tsup7_5Mpc2.bin','hy',0,0)

  !a=create_map_ter("x",15,"d","vxew",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l15_high.dat','./maps/high_res/velocity/15f15_analysis/map_high_15f15_x_map_vdx_ew_Tsup7_5Mpc2.bin','hy',0,0)
  !a=create_map_ter("y",15,"d","vyew",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l15_high.dat','./maps/high_res/velocity/15f15_analysis/map_high_15f15_y_map_vdy_ew_Tsup7_5Mpc2.bin','hy',0,0)
  !a=create_map_ter("z",15,"d","vzew",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l15_high.dat','./maps/high_res/velocity/15f15_analysis/map_high_15f15_z_map_vdz_ew_Tsup7_5Mpc2.bin','hy',0,0)
  !a=create_map_ter("z",15,"d","vzew",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l15_MW_los.dat','./maps/high_res/velocity/15f15_analysis/map_high_15f15_cen_map_vdz_ew_Tsup7_5Mpc2.bin','hy',0,0)

  !a=create_map_ter("x",15,"P","vxmw",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l15_high_gal_clean_m1e8.5.dat','./maps/high_res/velocity/15f15_analysis/map_high_15f15_x_map_P_mw_Tsup5_gal8.5_10Mpc2.bin','hy',0,10,tcut)
  !a=create_map_ter("y",15,"P","vymw",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l15_high_gal_clean_m1e8.5.dat','./maps/high_res/velocity/15f15_analysis/map_high_15f15_y_map_P_mw_Tsup5_gal8.5_10Mpc2.bin','hy',0,10,tcut)
  !a=create_map_ter("z",15,"P","vzmw",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l15_high_gal_clean_m1e8.5.dat','./maps/high_res/velocity/15f15_analysis/map_high_15f15_z_map_P_mw_Tsup5_gal8.5_10Mpc2.bin','hy',0,10,tcut)
  !a=create_map_ter("z",15,"P","vzmw",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l15_high_gal_clean_m1e8.5_MW_los.dat','./maps/high_res/velocity/15f15_analysis/map_high_15f15_cen_map_P_mw_Tsup5_gal8.5_10Mpc2.bin','hy',0,10,tcut)
  tcut = 8.6142e-1
  !tcut = 0.2
  !a=create_map_ter("x",15,"T","mw",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l15_high.dat','./maps/high_res/velocity/15f15_analysis/map_high_15f15_x_map_T_mw_Tsup7_5Mpc2.bin','hy',0,5,tcut)
  !a=create_map_ter("y",15,"n","mw",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l15_high.dat','./maps/high_res/velocity/15f15_analysis/map_high_15f15_y_map_ne_mw_Tsup_0.2kev_15Mpc2.bin','hy',0,15,tcut)
  !a=create_map_ter("z",15,"n","mw",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l15_high.dat','./maps/high_res/velocity/15f15_analysis/map_high_15f15_z_map_ne_mw_Tsup_0.2kev_15Mpc2.bin','hy',0,15,tcut)
  !a=create_map_ter("z",15,"n","mw",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l15_high_MW_los.dat','./maps/high_res/velocity/15f15_analysis/map_high_15f15_cen_map_ne_mw_Tsup_0.2kev_15Mpc2.bin','hy',0,15,tcut)
  !a=create_map_ter("x",15,"n","mw",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l15_high_gal_clean_m1e8.5.dat','./maps/high_res/velocity/15f15_analysis/map_high_15f15_x_map_ne_mw_Tsup5_gal8.5_10Mpc2.bin','hy',0,10,tcut)

  a=create_map_ter("x",15,"T","sl",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l15_high.dat','./maps/high_res/velocity/15f15_analysis/map_high_15f15_x_map_T_sl_Tsup7_5Mpc2.bin','hy',0,5,tcut)
  a=create_map_ter("y",15,"T","sl",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l15_high.dat','./maps/high_res/velocity/15f15_analysis/map_high_15f15_y_map_T_sl_Tsup7_5Mpc2.bin','hy',0,5,tcut)
  a=create_map_ter("z",15,"T","sl",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l15_high.dat','./maps/high_res/velocity/15f15_analysis/map_high_15f15_z_map_T_sl_Tsup7_5Mpc2.bin','hy',0,5,tcut)
  a=create_map_ter("z",15,"T","sl",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l15_high_MW_los.dat','./maps/high_res/velocity/15f15_analysis/map_high_15f15_cen_map_T_sl_Tsup7_5Mpc2.bin','hy',0,5,tcut)

  stop

  do i=1,10
    write(*,*) "i",i
    if (i<10) then
      write(i_str,'(I1)') i
      i_str='0'//trim(i_str)
    else
      write(i_str,'(I2)') i

    end if
    !file_len = len('./maps/high_res/filament/T_slices/x/map_high_19_xyz_full_map_T_x-') + len('_d0.08.bin') + 2
    file_len = len('./maps/high_res/velocity/15f15_analysis/random_proj/map_high_15f15_map_v_ew_Tsup7_5Mpc2.bin') + 2
    end do
    allocate(character(file_len) :: titre)
    titre = trim('./maps/high_res/velocity/15f15_analysis/random_proj/map_high_15f15_map_v_')//trim(i_str)//trim('ew_Tsup7_5Mpc2.bin')
    write(*,*) titre
    tcut = 8.6142e-1
    a=create_map_ter("z",15,"v","vzew",'/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l15_high.dat',titre,'hy',0,10,tcut)

  !  deallocate(titre)
  !enddo
  !i=r200_subhalos()

  !a=create_map_3D("z",15,"P","mw",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l19_fil_los.dat','./maps/3D_pyvista_lvl15_without_gals_P_fil_los.bin')

  a=create_map_3D_bis("z",15,"v","vx",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21.dat','./maps/high_res/filament/3D_fields/map_3D_high_15_xyz_right_minus_0Mpc_vx_10x10x20Mpc_cube_rdr15.bin','hy',1)
  a=create_map_3D_bis("z",15,"v","vy",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21.dat','./maps/high_res/filament/3D_fields/map_3D_high_15_xyz_right_minus_0Mpc_vy_10x10x20Mpc_cube_rdr15.bin','hy',1)
  a=create_map_3D_bis("z",15,"v","vz",'/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21.dat','./maps/high_res/filament/3D_fields/map_3D_high_15_xyz_right_minus_0Mpc_vz_10x10x20Mpc_cube_rdr15.bin','hy',1)
  stop

 
  
  !ncell=166999251
  !ncell=338012253
  

  !open(1,file='/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21_gal_clean.dat',form='unformatted')
  !open(2,file='/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21_gal_clean_ncell.dat',form='unformatted')
  open(1,file='/data/cluster/tlebeau/virgo/virgo_xyz_dm_high_res.dat',form='unformatted')
  read(1) ncell
  write(*,*) "ncell",ncell

  allocate(x(ncell),y(ncell),z(ncell),vx(ncell),vy(ncell),vz(ncell),m(ncell),r(ncell))
  
  !read(1) n_e
  !write(*,*) size(n_e)
  !read(1) T
  !read(1) P
  read(1) x
  read(1) y
  read(1) z
  !read(1) vx
  !read(1) vy
  !read(1) vz
  !read(1) m
  !read(1) lvl
 
  !write(2) ncell
  !write(2) n_e
  !write(2) T
  !write(2) P
  !write(2) x
  !write(2) y
  !write(2) z
  !write(2) vx
  !write(2) vy
  !write(2) vz
  !write(2) m
  !write(2) lvl

  close(1)

  write(*,*) "data loaded"
  !stop
  !close(2)

  x_cen = 0.48461068
  y_cen = 0.50809848
  z_cen = 0.49687076

  x_cen = (x_cen - 0.5) * (unit_l / 3.08567758128E21)
  y_cen = (y_cen - 0.5) * (unit_l / 3.08567758128E21)
  z_cen = (z_cen - 0.5) * (unit_l / 3.08567758128E21)

  
  r = sqrt((x - x_cen) ** 2 + (y - y_cen) ** 2 + (z - z_cen) ** 2)

  mass_tensor=0

  do i=1,size(x)
    if (r(i)<2187) then
      mass_tensor(1,1)=mass_tensor(1,1)+x(i)*x(i)
      mass_tensor(1,2)=mass_tensor(1,2)+x(i)*y(i)
      mass_tensor(1,3)=mass_tensor(1,3)+x(i)*z(i)
      mass_tensor(2,2)=mass_tensor(1,1)+y(i)*y(i)
      mass_tensor(3,3)=mass_tensor(1,1)+z(i)*z(i)
      mass_tensor(2,3)=mass_tensor(1,1)+y(i)*z(i)
    end if
  end do

  mass_tensor(2,1)=mass_tensor(1,2)
  mass_tensor(3,1)=mass_tensor(1,3)
  mass_tensor(3,2)=mass_tensor(2,3)

  write(*,*) "mass tensor",mass_tensor
  stop
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

  !write(*,*) m_rad
  !write(*,*) n_rad

  !open(3,file='/data/cluster/tlebeau/virgo/rad_log_dm_21.dat',form='unformatted')
  !write(3) m_rad
  !write(3) n_rad

  !close(3)
     
 
  
end program mainf90
