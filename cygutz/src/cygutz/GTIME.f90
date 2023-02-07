module gtime
    use gprec
    implicit none
    private

    type,public::time_ob
        real::ta(2,7) ! Depth of timing tags.
        real(q)::tb(2,7)

        contains
        procedure,public::start=>start_time
        procedure,public::finish=>finish_time
        procedure,public::print_usage=>print_usage
    end type time_ob

    contains


    subroutine start_time(this, i)
    class(time_ob) :: this
    integer,intent(in)::i

    call set_time(this,1,i)
    return

    end subroutine start_time


    subroutine finish_time(this, i)
    class(time_ob) :: this
    integer,intent(in)::i

    call set_time(this,2,i)
    return

    end subroutine finish_time


    subroutine print_usage(this,task,i,io)
    class(time_ob) :: this
    character(*),intent(in)::task
    integer,intent(in)::i,io

    integer nh,nm
    real(q) time
    real(q),parameter::s_per_m=60._4,s_per_h=3600._4
    character(100) fmt

    if(io<0)return
    time=this%ta(2,i)-this%ta(1,i)
    nh=int(time/s_per_h); time=mod(time,s_per_h)
    nm=int(time/s_per_m); time=mod(time,s_per_m)
    write(fmt,'("(1x,a",i0,",a6,i5,a3,i2,a3,f5.2,a8)")')len(task)
    write(io,fmt)task," takes",nh," h ",nm," m ",time," s(cpu)."
    time=this%tb(2,i)-this%tb(1,i)
    nh=int(time/s_per_h); time=mod(time,s_per_h)
    nm=int(time/s_per_m); time=mod(time,s_per_m)
    write(fmt,'("(1x,",i0,"x,6x,i5,a3,i2,a3,f5.2,a9)")')len(task)
    write(io,fmt)nh," h ",nm," m ",time," s(wall)."
    return

    end subroutine print_usage


    subroutine set_time(this,it,i)
    class(time_ob) :: this
    integer,intent(in)::i,it

    integer::tib,tirate

    call cpu_time(this%ta(it,i))
    call system_clock(tib,tirate)
    this%tb(it,i)=dble(tib)/dble(tirate)
    return

    end subroutine set_time


end module gtime
