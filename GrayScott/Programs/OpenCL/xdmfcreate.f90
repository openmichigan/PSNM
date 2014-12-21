	PROGRAM XdmfCreate
        !--------------------------------------------------------------------------------
        ! .. Purpose ..
        ! XdmfCreate is a postprocessing program which creates header files for
        ! Paraview and Visit
        ! It uses the INPUTFILE and assumes that the filenames in the program
        ! are
        ! consistent with those in the current file.
        !
        ! .. PARAMETERS .. INITIALIZED IN INPUTFILE
        ! Nx = power of two, number of modes in the x direction
        ! Ny = power of two, number of modes in the y direction
        ! Nz = power of two, number of modes in the z direction
        ! Nt = the number of timesteps
        ! plotgap = the number of timesteps to take before plotting
        ! Lx = definition of the periodic domain of computation in x direction
        ! Ly = definition of the periodic domain of computation in y direction
        ! Lz = definition of the periodic domain of computation in z direction
        ! Dt = the time step
        !
        ! REFERENCES
        !
        ! www.xdmf.org
        ! Paraview users mailing list
        ! VisIt users mailing list
        !
        ! ACCURACY
        !
        ! ERROR INDICATORS AND WARNINGS
        !
        ! FURTHER COMMENTS
        !------------------------------------------------------------------------------------
        ! EXTERNAL ROUTINES REQUIRED
        IMPLICIT NONE
        ! .. Scalar Arguments ..
        INTEGER(KIND=4) :: Nx, Ny, Nz, Nt, plotgap
        REAL(KIND=8)    :: Lx, Ly, Lz, DT
        ! .. Local scalars ..
        INTEGER(KIND=4) :: stat,plotnum,ind,n,numplots
        REAL(KIND=8), PARAMETER :: pi=3.14159265358979323846264338327950288419716939937510d0
        ! .. Local Arrays ..
        CHARACTER*50    :: InputFileName, OutputFileName, OutputFileName2
        CHARACTER*10    :: number_file
        InputFileName='INPUTFILE'
        OPEN(unit=11,FILE=trim(InputFileName),status="OLD")
        REWIND(11)
        READ(11,*) Nx, Ny, Nz, Nt, plotgap, Lx, Ly, Lz, DT
        CLOSE(11)       
        plotnum=0
        numplots=1+Nt/plotgap
        OutputFileName='udata.xmf'
        OPEN(unit=11,FILE=trim(OutputFileName),status="UNKNOWN")
        REWIND(11)
        WRITE(11,'(a)') "<?xml version=""1.0"" ?>"
        WRITE(11,'(a)') "<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>"
        WRITE(11,'(a)') "<Xdmf xmlns:xi=""http://www.w3.org/2001/XInclude"" Version=""2.0"">"
        WRITE(11,'(a)') "<Domain>"
        WRITE(11,'(a)') " <Topology name=""topo"" TopologyType=""3DCoRectMesh"""
        WRITE(11,*) "  Dimensions=""",Nx,Ny,Nz,""">"
        WRITE(11,'(a)') " </Topology>"
        WRITE(11,'(a)') " <Geometry name=""geo"" Type=""ORIGIN_DXDYDZ"">"
        WRITE(11,'(a)') "  <!-- Origin -->"
        WRITE(11,'(a)') "  <DataItem Format=""XML"" Dimensions=""3"">"
        WRITE(11,*) -pi*Lx, -pi*Ly, -pi*Lz 
        WRITE(11,'(a)') "  </DataItem>"
        WRITE(11,'(a)') "  <!-- DxDyDz -->"
        WRITE(11,'(a)') "  <DataItem Format=""XML"" Dimensions=""3"">"
        WRITE(11,*) REAL(2.0*pi*Lx/Nx,kind(0d0)), REAL(2.0*pi*Ly/Ny,kind(0d0)), REAL(2.0*pi*Lz/Nz,kind(0d0))  
        WRITE(11,'(a)') "   </DataItem>"
        WRITE(11,'(a)') "  </Geometry>"
        WRITE(11,'(a)')
        WRITE(11,'(a)') "       <Grid Name=""TimeSeries"" GridType=""Collection"" CollectionType=""Temporal"">"
        WRITE(11,'(a)') "               <Time TimeType=""HyperSlab"">"
        WRITE(11,'(a)') "                       <DataItem Format=""XML"" NumberType=""Float"" Dimensions=""3"">"
        WRITE(11,*) "                   0.0", dt," 2" 
        WRITE(11,'(a)') "                       </DataItem>"
        WRITE(11,'(a)') "               </Time>"

        DO n=1,numplots

                OutputFileName = 'data/u'
                ind = index(OutputFileName,' ') - 1
                WRITE(number_file,'(i0)') 10000000+plotnum
                OutputFileName = OutputFileName(1:ind)//number_file
                ind = index(OutputFileName,' ') - 1
                OutputFileName = OutputFileName(1:ind)//'.datbin'
                
                OutputFileName2 = "<Grid Name=""T"
                WRITE(number_file,'(i0)') n 
                OutputFileName2 = OutputFileName2(1:13)//number_file
                ind =  index(number_file,' ') - 1
                OutputFileName2 = OutputFileName2(1:13+ind)//'" GridType="Uniform">'
                plotnum=plotnum+1
                
                WRITE(11,'(a)')
                WRITE(11,'(3X,a)') "   ",OutputFileName2
                WRITE(11,'(a)') "       <Topology Reference=""/Xdmf/Domain/Topology[1]""/>"
                WRITE(11,'(a)') "       <Geometry Reference=""/Xdmf/Domain/Geometry[1]""/>"
                WRITE(11,'(a)') "        <Attribute Name=""Displacement"" Center=""Node"">"
                WRITE(11,'(a)') "               <DataItem Format=""Binary"" "
                WRITE(11,'(a)') "                DataType=""Float"" Precision=""8"" Endian=""Little"" "
                WRITE(11,*) "            Dimensions=""",Nx, Ny, Nz,""">"
                WRITE(11,'(16X,a)') "           ",OutputFileName
                WRITE(11,'(a)') "               </DataItem>"
                WRITE(11,'(a)') "       </Attribute>"
                WRITE(11,'(3X,a)') "</Grid>"
        END DO
        WRITE(11,'(a)')"        </Grid>"
        WRITE(11,'(a)') "</Domain>"
        WRITE(11,'(a)') "</Xdmf>"
        CLOSE(11)
		!same for v
        plotnum=0
        numplots=1+Nt/plotgap
        OutputFileName='vdata.xmf'
        OPEN(unit=11,FILE=trim(OutputFileName),status="UNKNOWN")
        REWIND(11)
        WRITE(11,'(a)') "<?xml version=""1.0"" ?>"
        WRITE(11,'(a)') "<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>"
        WRITE(11,'(a)') "<Xdmf xmlns:xi=""http://www.w3.org/2001/XInclude"" Version=""2.0"">"
        WRITE(11,'(a)') "<Domain>"
        WRITE(11,'(a)') " <Topology name=""topo"" TopologyType=""3DCoRectMesh"""
        WRITE(11,*) "  Dimensions=""",Nx,Ny,Nz,""">"
        WRITE(11,'(a)') " </Topology>"
        WRITE(11,'(a)') " <Geometry name=""geo"" Type=""ORIGIN_DXDYDZ"">"
        WRITE(11,'(a)') "  <!-- Origin -->"
        WRITE(11,'(a)') "  <DataItem Format=""XML"" Dimensions=""3"">"
        WRITE(11,*) -pi*Lx, -pi*Ly, -pi*Lz 
        WRITE(11,'(a)') "  </DataItem>"
        WRITE(11,'(a)') "  <!-- DxDyDz -->"
        WRITE(11,'(a)') "  <DataItem Format=""XML"" Dimensions=""3"">"
        WRITE(11,*) REAL(2.0*pi*Lx/Nx,kind(0d0)), REAL(2.0*pi*Ly/Ny,kind(0d0)), REAL(2.0*pi*Lz/Nz,kind(0d0))  
        WRITE(11,'(a)') "   </DataItem>"
        WRITE(11,'(a)') "  </Geometry>"
        WRITE(11,'(a)')
        WRITE(11,'(a)') "       <Grid Name=""TimeSeries"" GridType=""Collection"" CollectionType=""Temporal"">"
        WRITE(11,'(a)') "               <Time TimeType=""HyperSlab"">"
        WRITE(11,'(a)') "                       <DataItem Format=""XML"" NumberType=""Float"" Dimensions=""3"">"
        WRITE(11,*) "                   0.0", dt," 2" 
        WRITE(11,'(a)') "                       </DataItem>"
        WRITE(11,'(a)') "               </Time>"

        DO n=1,numplots

                OutputFileName = 'data/v'
                ind = index(OutputFileName,' ') - 1
                WRITE(number_file,'(i0)') 10000000+plotnum
                OutputFileName = OutputFileName(1:ind)//number_file
                ind = index(OutputFileName,' ') - 1
                OutputFileName = OutputFileName(1:ind)//'.datbin'
                
                OutputFileName2 = "<Grid Name=""T"
                WRITE(number_file,'(i0)') n 
                OutputFileName2 = OutputFileName2(1:13)//number_file
                ind =  index(number_file,' ') - 1
                OutputFileName2 = OutputFileName2(1:13+ind)//'" GridType="Uniform">'
                plotnum=plotnum+1
                
                WRITE(11,'(a)')
                WRITE(11,'(3X,a)') "   ",OutputFileName2
                WRITE(11,'(a)') "       <Topology Reference=""/Xdmf/Domain/Topology[1]""/>"
                WRITE(11,'(a)') "       <Geometry Reference=""/Xdmf/Domain/Geometry[1]""/>"
                WRITE(11,'(a)') "        <Attribute Name=""Displacement"" Center=""Node"">"
                WRITE(11,'(a)') "               <DataItem Format=""Binary"" "
                WRITE(11,'(a)') "                DataType=""Float"" Precision=""8"" Endian=""Little"" "
                WRITE(11,*) "            Dimensions=""",Nx, Ny, Nz,""">"
                WRITE(11,'(16X,a)') "           ",OutputFileName
                WRITE(11,'(a)') "               </DataItem>"
                WRITE(11,'(a)') "       </Attribute>"
                WRITE(11,'(3X,a)') "</Grid>"
        END DO
        WRITE(11,'(a)')"        </Grid>"
        WRITE(11,'(a)') "</Domain>"
        WRITE(11,'(a)') "</Xdmf>"
        CLOSE(11)
        END PROGRAM XdmfCreate
