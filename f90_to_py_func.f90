module f90_to_py
      use const
      implicit none

contains

      subroutine dm_mass_tensor(ncell,mass_tensor) !,y,z,r)
            doubleprecision, dimension(ncell) :: x,y,z,r
            !real(kind=xp), dimension(ncell) :: databis
            doubleprecision :: x_cen,y_cen,z_cen
            integer, intent(in) :: ncell
            integer :: i,counter
            integer :: ncell_confirm
            doubleprecision, dimension(3,3), intent(out) :: mass_tensor

            open(1,file='/data/cluster/tlebeau/virgo/virgo_xyz_dm_low_res.dat',form='unformatted')
            read(1) ncell_confirm

            read(1) x
            read(1) y
            read(1) z

            close(1)

            x_cen = 0.48461068
            y_cen = 0.50809848
            z_cen = 0.49687076
            x_cen = (x_cen - 0.5) * (unit_l / 3.08567758128E21)
            y_cen = (y_cen - 0.5) * (unit_l / 3.08567758128E21)
            z_cen = (z_cen - 0.5) * (unit_l / 3.08567758128E21)

            x=x-x_cen
            y=y-y_cen
            z=z-z_cen

            x_cen=0
            y_cen=0
            z_cen=0

            r = sqrt((x - x_cen) ** 2 + (y - y_cen) ** 2 + (z - z_cen) ** 2)

            mass_tensor=0

            counter=0

            do i=1,size(x)
                  if (r(i)<2187) then
                        counter=counter+1
                        mass_tensor(1,1)=mass_tensor(1,1)+x(i)*x(i)
                        mass_tensor(1,2)=mass_tensor(1,2)+x(i)*y(i)
                        mass_tensor(1,3)=mass_tensor(1,3)+x(i)*z(i)
                        mass_tensor(2,2)=mass_tensor(2,2)+y(i)*y(i)
                        mass_tensor(3,3)=mass_tensor(3,3)+z(i)*z(i)
                        mass_tensor(2,3)=mass_tensor(2,3)+y(i)*z(i)
                  end if
            end do

            write(*,*) "counter",counter
            stop

            mass_tensor(2,1)=mass_tensor(1,2)
            mass_tensor(3,1)=mass_tensor(1,3)
            mass_tensor(3,2)=mass_tensor(2,3)

            !write(*,*) "mass tensor",mass_tensor

      end subroutine dm_mass_tensor

      subroutine kinetic_energy(ncell_dm,ncell_hydro,rlim,T)
            integer, intent(in) :: ncell_dm,ncell_hydro
            doubleprecision, dimension(ncell_dm) :: x_dm,y_dm,z_dm,vx_dm,vy_dm,vz_dm,r_dm,m_dm
            integer :: ncell_confirm,i
            doubleprecision, intent(out) :: T
            doubleprecision, intent(in) :: rlim
            doubleprecision :: v_square, x_cen,y_cen,z_cen

            open(1,file='/data/cluster/tlebeau/virgo/virgo_xyz_dm_high_res.dat',form='unformatted')
            read(1) ncell_confirm

            read(1) x_dm
            read(1) y_dm
            read(1) z_dm
            read(1) vx_dm
            read(1) vy_dm
            read(1) vz_dm
            read(1) m_dm

            close(1)

            T=0

            m_dm=m_dm*m_sun
            vx_dm=vx_dm*1e3
            vy_dm=vy_dm*1e3
            vz_dm=vz_dm*1e3

            x_cen = 0.48461068
            y_cen = 0.50809848
            z_cen = 0.49687076
            x_cen = (x_cen - 0.5) * (unit_l / 3.08567758128E21)
            y_cen = (y_cen - 0.5) * (unit_l / 3.08567758128E21)
            z_cen = (z_cen - 0.5) * (unit_l / 3.08567758128E21)

            r_dm = sqrt((x_dm - x_cen) ** 2 + (y_dm - y_cen) ** 2 + (z_dm - z_cen) ** 2)


            do i=1,ncell_dm
                  if(r_dm(i)<rlim) then
                        v_square=vx_dm(i)**2+vy_dm(i)**2+vz_dm(i)**2
                        T=T+m_dm(i)*v_square
                  end if
            end do

            T=T/2

            write(*,*) "kinetic energy",T




      end subroutine kinetic_energy

      subroutine potential_energy(ncell_dm,ncell_hydro,n_subbox,width,rlim,U)
            integer, intent(in) :: ncell_dm,ncell_hydro,n_subbox
            doubleprecision, intent(in)::width,rlim
            doubleprecision, dimension(ncell_dm) :: x_dm,y_dm,z_dm,r_dm,vx_dm,vy_dm,vz_dm,m_dm
            integer :: ncell_confirm,i,j,k,l
            doubleprecision, intent(out) :: U
            doubleprecision :: dist,x_cen,y_cen,z_cen,len,minx,maxx,miny,maxy,minz,maxz,dx
            doubleprecision, dimension(n_subbox,n_subbox,n_subbox,4) :: sub_boxes
            doubleprecision, dimension(n_subbox**3,4) :: sub_boxes_list
            doubleprecision :: sum_mdm_rvir

            open(1,file='/data/cluster/tlebeau/virgo/virgo_xyz_dm_high_res.dat',form='unformatted')
            read(1) ncell_confirm

            read(1) x_dm
            read(1) y_dm
            read(1) z_dm
            read(1) vx_dm
            read(1) vy_dm
            read(1) vz_dm
            read(1) m_dm

            close(1)

            write(*,*) "data loaded"

            !write(*,*) "dm mass sum" , sum(m_dm)
            !write(*,*) "dm mass mean", sum(m_dm)/size(m_dm)

            U=0

            x_cen = 0.48461068
            y_cen = 0.50809848
            z_cen = 0.49687076
            x_cen = (x_cen - 0.5) * (unit_l / 3.08567758128E21)
            y_cen = (y_cen - 0.5) * (unit_l / 3.08567758128E21)
            z_cen = (z_cen - 0.5) * (unit_l / 3.08567758128E21)

            r_dm = sqrt((x_dm - x_cen) ** 2 + (y_dm - y_cen) ** 2 + (z_dm - z_cen) ** 2)

            !m_dm=m_dm*m_sun
            !x_dm=x_dm*pc*1e3
            !y_dm=y_dm*pc*1e3
            !z_dm=z_dm*pc*1e3

            len=2*width
            minx=x_cen-width
            maxx=x_cen+width
            miny=y_cen-width
            maxy=y_cen+width
            minz=z_cen-width
            maxz=z_cen+width
            dx=len/n_subbox

            sub_boxes=0
            sub_boxes_list=0
            U=0

            !do k=1,n_subbox
            !      sub_boxes(k,:,:,1)=minx+dx*(0.5+(k-1))
            !      sub_boxes(:,k,:,2)=miny+dx*(0.5+(k-1))
            !      sub_boxes(:,:,k,3)=minz+dx*(0.5+(k-1))
            !end do

            !write(*,*) "z positions", sub_boxes(:,:,:,3)

            !stop

            !write(*,*) "positions calculated"



            sum_mdm_rvir=0

            do k=1,ncell_dm
                  if(r_dm(k)<rlim) then
                        i=int(((x_dm(k)-minx)/len)*n_subbox)
                        j=int(((y_dm(k)-miny)/len)*n_subbox)
                        l=int(((z_dm(k)-minz)/len)*n_subbox)
                        sub_boxes(i,j,l,4)=sub_boxes(i,j,l,4)+m_dm(k)
                        sub_boxes(i,j,l,1)=sub_boxes(i,j,l,1)+m_dm(k)*x_dm(k)
                        sub_boxes(i,j,l,2)=sub_boxes(i,j,l,2)+m_dm(k)*y_dm(k)
                        sub_boxes(i,j,l,3)=sub_boxes(i,j,l,3)+m_dm(k)*z_dm(k)
                        sum_mdm_rvir=sum_mdm_rvir+m_dm(k)
                  end if
            end do

            !write(*,*) "sum mass m_dm ds rvir",sum_mdm_rvir
            !write(*,*) "sum mass subboxes", sum(sub_boxes(:,:,:,4))
            !write(*,*) "ratio", sum_mdm_rvir/sum(sub_boxes(:,:,:,4))

            !stop

            !write(*,*) "masses in subboxes calculated", sub_boxes(:,:,:,4)

            !stop

            !sub_boxes(:,:,:,1)=sub_boxes(:,:,:,1)/sub_boxes(:,:,:,4)
            !sub_boxes(:,:,:,2)=sub_boxes(:,:,:,2)/sub_boxes(:,:,:,4)
            !sub_boxes(:,:,:,3)=sub_boxes(:,:,:,3)/sub_boxes(:,:,:,4)


            !write(*,*) "masses in subboxes calculated", sub_boxes(:,1,1,1)

            !stop

            !write(*,*) "masses",sub_boxes(:,:,:,4)

            do i=1,n_subbox
                  do j=1,n_subbox
                        do l=1,n_subbox
                              if (sub_boxes(i,j,l,4)>0) then
                                    !write(*,*) "in"
                                    sub_boxes(i,j,l,1)=sub_boxes(i,j,l,1)/sub_boxes(i,j,l,4)
                                    sub_boxes(i,j,l,2)=sub_boxes(i,j,l,2)/sub_boxes(i,j,l,4)
                                    sub_boxes(i,j,l,3)=sub_boxes(i,j,l,3)/sub_boxes(i,j,l,4)
                              !else
                              !      sub_boxes(k,:,:,1)=minx+dx*(0.5+(k-1))
                              !      sub_boxes(:,k,:,2)=miny+dx*(0.5+(k-1))
                              !      sub_boxes(:,:,k,3)=minz+dx*(0.5+(k-1))

                              end if
                              k=i+n_subbox*(j-1)+(n_subbox**2)*(l-1)
                              sub_boxes_list(k,:)=sub_boxes(i,j,l,:)
                        end do
                  end do
                  !write(*,*) "i",i
            end do

            !write(*,*) "pos 1", sub_boxes(1,1,1,:)
            !write(*,*) "pos 2", sub_boxes(1,1,2,:)
            !write(*,*) "pos 3", sub_boxes(1,1,3,:)

            !write(*,*) "x sub boxes", sub_boxes(:,:,:,1)
            !write(*,*) "y sub boxes", sub_boxes(:,:,:,2)
            !write(*,*) "z sub boxes", sub_boxes(:,:,:,3)
            !write(*,*) "m sub boxes", sub_boxes(:,:,:,4)

            !write(*,*) "m sub_boxes_list",sub_boxes_list(:,4)


            !stop

            !write(*,*) "subboxes list created" , sub_boxes_list(:,1)


            do i=1,size(sub_boxes_list(:,1))
                  !write(*,*) "i",i
                  do j=i+1,size(sub_boxes_list(:,1))
                        if(sub_boxes_list(i,4)>0 .and. sub_boxes_list(j,4)>0) then
                        dist=sqrt((sub_boxes_list(i,1)-sub_boxes_list(j,1))**2+(sub_boxes_list(i,2)- &
                        sub_boxes_list(j,2))**2+(sub_boxes_list(i,3)-sub_boxes_list(j,3))**2)
                        !write(*,*) "dist", dist
                        U=U+(sub_boxes_list(i,4)*sub_boxes_list(j,4))/dist
                        end if
                        !write(*,*) "U",U
                        !if (j==20) then
                        !      stop
                        !end if

                  end do
            end do


            !do i=1,ncell_dm
            !      if (MOD(i,int(1e4))==0) then
            !            write(*,*) "ratio", float(i)/float(ncell_dm)
            !            !write(*,*) "U",U
            !      end if

            !      if(r_dm(i)<5031) then
            !            write(*,*) "in"
            !            do j=i+1,ncell_dm
            !                  if (MOD(j,int(1e6))==0) then
            !                  !write(*,*) "ratio", float(i)/float(ncell_dm)
            !                  !write(*,*) "U",U
            !                  end if
            !                  if(r_dm(j)<5031) then
            !                        dist=sqrt((x_dm(i)-x_dm(j))**2+(y_dm(i)-y_dm(j))**2+(z_dm(i)-z_dm(j))**2)
            !                        U=U+(m_dm(i)*m_dm(j))/dist
            !                        !write(*,*) "in 2"
            !                  end if
            !            end do
            !            !write(*,*) "out"
            !            !stop

            !      end if

            !end do

            !write(*,*) "U",U

            U=-U*G
            U=U*(m_sun**2/(pc*1e3))

            write(*,*) "potential energy",U


      end subroutine potential_energy

      subroutine find_minimum_potential(ncell_dm,x_min,y_min,z_min)
            integer, intent(in) :: ncell_dm
            doubleprecision, intent(out) :: x_min,y_min,z_min
            doubleprecision, dimension(ncell_dm) :: x_dm,y_dm,z_dm,r_dm,vx_dm,vy_dm,vz_dm,m_dm
            integer :: ncell_confirm,i,n_it,j
            doubleprecision :: x_test,y_test,z_test,dx,dy,dz,searching_radius,phi,phi_min,x_sum,y_sum,z_sum,m_sum


            open(1,file='/data/cluster/tlebeau/virgo/virgo_xyz_dm_high_res.dat',form='unformatted')
            read(1) ncell_confirm

            read(1) x_dm
            read(1) y_dm
            read(1) z_dm
            read(1) vx_dm
            read(1) vy_dm
            read(1) vz_dm
            read(1) m_dm

            close(1)

            write(*,*) "data loaded"

            x_min = 0.48461068
            y_min = 0.50809848
            z_min = 0.49687076
            x_min = (x_min - 0.5) * (unit_l / 3.08567758128E21)
            y_min = (y_min - 0.5) * (unit_l / 3.08567758128E21)
            z_min = (z_min - 0.5) * (unit_l / 3.08567758128E21)

            r_dm = sqrt((x_dm - x_min) ** 2 + (y_dm - y_min) ** 2 + (z_dm - z_min) ** 2)

            x_sum=0
            y_sum=0
            z_sum=0
            m_sum=0

            !do i=1,ncell_dm
            !      if (r_dm<5031) then
            !      x_sum = x_sum + x_dm(i)*m_dm(i)
            !      y_sum = y_sum + y_dm(i)*m_dm(i)
            !      z_sum = z_sum + z_dm(i)*m_dm(i)
            !      m_sum = m_sum + m_dm(i)
            !      end if
            !end do

            !x_sum=x_sum/m_sum
            !y_sum=y_sum/m_sum
            !z_sum=z_sum/m_sum

            !write(*,*) "x_com = ",x_sum
            !write(*,*) "y_com = ",y_sum
            !write(*,*) "z_com = ",z_sum

            !stop

            !x_min = 0.48461068
            !y_min = 0.50809848
            !z_min = 0.49687076
            !x_min = (x_min - 0.5) * (unit_l / 3.08567758128E21)
            !y_min = (y_min - 0.5) * (unit_l / 3.08567758128E21)
            !z_min = (z_min - 0.5) * (unit_l / 3.08567758128E21)

            !Après 10000 itérations :

            x_min = -11334.407242699199
            y_min = 5964.429627777705
            z_min = -2303.1811641473237

            !Après avoir recalculé le centre de masse :

            !x_min = -11653.577688278487
            !y_min = 6274.5095544073174
            !z_min = -2197.2118388121430


            r_dm = sqrt((x_dm - x_min) ** 2 + (y_dm - y_min) ** 2 + (z_dm - z_min) ** 2)

            do j=1,ncell_dm
                  phi_min=phi_min-G*m_dm(j)/r_dm(j)
            end do

            write(*,*) "initial potential",phi_min


            searching_radius=10

            n_it=100


            do i=1,n_it
                  write(*,*) "i",i
                  phi=0
                  r_dm=0
                  x_test=0
                  y_test=0
                  z_test=0
                  !if(mod(i,int(n_it/5))==0) then
                  !      searching_radius=searching_radius*(1-0.2)
                  !      write(*,*) "i",i,"searching radius",searching_radius
                  !end if
                  call random_number(dx)
                  call random_number(dy)
                  call random_number(dz)
                  dx=(dx-0.5)*searching_radius
                  dy=(dy-0.5)*searching_radius
                  dz=(dz-0.5)*searching_radius

                  x_test=x_min+dx
                  y_test=y_min+dy
                  z_test=z_min+dz

                  r_dm = sqrt((x_dm - x_test) ** 2 + (y_dm - y_test) ** 2 + (z_dm - z_test) ** 2)

                  do j=1,ncell_dm
                        phi=phi-G*m_dm(j)/r_dm(j)
                  end do

                  write(*,*) "phi",phi

                  if(phi<phi_min) then
                        write(*,*) "Found lower potential: ",phi," at x,y,z:",x_test,y_test,z_test
                        write(*,*) "Distance to previous minimum: ",sqrt((x_min-x_test)**2+(y_min-y_test)**2+(z_min-z_test)**2)
                        x_min=x_test
                        y_min=y_test
                        z_min=z_test
                        phi_min=phi
                  end if




            end do

      end subroutine find_minimum_potential

      subroutine read_dm(ncell,file_name,x,y,z,m,vx,vy,vz)
            character(len=*), intent(in)::file_name
            doubleprecision, dimension(ncell) ,intent(out):: x,y,z,m,vx,vy,vz
            integer, intent(in) :: ncell
            integer :: ncell_confirm
            doubleprecision :: x_cen,y_cen,z_cen

            open(1,file=file_name,form='unformatted')
            read(1) ncell_confirm

            read(1) x
            read(1) y
            read(1) z
            read(1) vx
            read(1) vy
            read(1) vz
            read(1) m


            close(1)

            x_cen = 0.48461068
            y_cen = 0.50809848
            z_cen = 0.49687076
            x_cen = (x_cen - 0.5) * (unit_l / 3.08567758128E21)
            y_cen = (y_cen - 0.5) * (unit_l / 3.08567758128E21)
            z_cen = (z_cen - 0.5) * (unit_l / 3.08567758128E21)

            x=x-x_cen
            y=y-y_cen
            z=z-z_cen

      end subroutine read_dm

       subroutine read_map_file(ncell,file_name,map,three_D)
            character(len=*), intent(in)::file_name
            doubleprecision, dimension(ncell),intent(out):: map !7598896696
            integer :: nx,ny,nz,three_D
            doubleprecision :: cen_x,cen_y,cen_z
            !doubleprecision, dimension(ncell) ::vx,vy,vz
            !real(kind=xp), dimension(ncell) :: databis
            !doubleprecision :: x_cen,y_cen,z_cen
            integer, intent(in) :: ncell
            !unsigned :: ncell
            !integer :: i
            !integer :: ncell_confirm
            !doubleprecision :: x_cen,y_cen,z_cen
            !doubleprecision, dimension(3,3), intent(out) :: mass_tensor

            !write(*,*) "ncell",ncell
            !allocate(map(ncell))
            !write(*,*) "len(map)",size(map)

            nx=0
            ny=0
            nz=0

            open(1,file=file_name,form='unformatted')
            read(1) nx,ny,nz
            !write(*,*) "nx,ny,nz",nx,ny,nz
            if (three_D==0) then
                  read(1) cen_x,cen_y,cen_z
            end if
            read(1) map

            !read(1) x
            !read(1) y
            !read(1) z
            !read(1) vx
            !read(1) vy
            !read(1) vz
            !read(1) m


            close(1)

            !write(*,*) "reading done"

            !x_cen = 0.48461068
            !y_cen = 0.50809848
            !z_cen = 0.49687076
            !x_cen = (x_cen - 0.5) * (unit_l / 3.08567758128E21)
            !y_cen = (y_cen - 0.5) * (unit_l / 3.08567758128E21)
            !z_cen = (z_cen - 0.5) * (unit_l / 3.08567758128E21)

            !x=x-x_cen
            !y=y-y_cen
            !z=z-z_cen

      end subroutine read_map_file

      subroutine read_vel_file(ncell,file_name,x,y,z,vx,vy,vz)
            character(len=*), intent(in)::file_name
            doubleprecision, dimension(ncell) ,intent(out):: x,y,z,vx,vy,vz
            integer, intent(in) :: ncell
            integer :: ncell_confirm
            doubleprecision :: x_cen,y_cen,z_cen

            open(1,file=file_name,form='unformatted')
            read(1) ncell_confirm

            read(1) x
            read(1) y
            read(1) z
            read(1) vx
            read(1) vy
            read(1) vz
            !read(1) m


            close(1)

            x_cen = 0.48461068
            y_cen = 0.50809848
            z_cen = 0.49687076
            x_cen = (x_cen - 0.5) * (unit_l / 3.08567758128E21)
            y_cen = (y_cen - 0.5) * (unit_l / 3.08567758128E21)
            z_cen = (z_cen - 0.5) * (unit_l / 3.08567758128E21)

            x=x-x_cen
            y=y-y_cen
            z=z-z_cen

      end subroutine read_vel_file

end module f90_to_py