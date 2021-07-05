### Purpose and working flow with codes under this path

This path contains the user (mod_usr.t) and parameter (amrvac.par) files for the PFSS model.

The files have to be used in combination with MPI-AMRVAC version 2.0+ (commit b7e967e), which can be found via [this link](https://github.com/amrvac/amrvac/tree/b7e967ecfbaa027a683fd54525f3a83cd0ad9251). The applicability to other versions needs further tests.

The boundary condition (synoptic map for the radial magnetic field) is provided in ./potential_boundary. The format should conform with the reading subroutine under MPI_AMRVAC/src/physics/mod_pfss.t
>        call MPI_FILE_OPEN(MPI_COMM_SELF,mapname,MPI_MODE_RDONLY,MPI_INFO_NULL,&
>                           file_handle,ierrmpi)
>        call MPI_FILE_READ(file_handle,xm,1,MPI_INTEGER,statuss,ierrmpi)
>        call MPI_FILE_READ(file_handle,ym,1,MPI_INTEGER,statuss,ierrmpi)
>        allocate(b_r0(xm,ym))
>        allocate(theta(ym))
>        allocate(phi(xm))
>        call MPI_FILE_READ(file_handle,theta,ym,MPI_DOUBLE_PRECISION,&
>                           statuss,ierrmpi)
>        call MPI_FILE_READ(file_handle,phi,xm,MPI_DOUBLE_PRECISION,&
>                           statuss,ierrmpi)
>        call MPI_FILE_READ(file_handle,b_r0,xm*ym,MPI_DOUBLE_PRECISION,&
>                           statuss,ierrmpi)
>        call MPI_FILE_CLOSE(file_handle,ierrmpi)

Note that xm is the size in the Phi direction, ym is that in the Theta direction. And one should set the following values in amrvac.par
>     domain_nx1 = ym
>     domain_nx2 = ym
>     domain_nx3 = xm
And note that block_nx1, block_nx2, block_nx3 should be factors of domain_nx1, domain_nx2, domain_nx3, respectively. And block_nx1, block_nx2, block_nx3 should be even numbers.
