       program writepdb

!     program to find the cations and their coordinated atoms in a PDB file
!     by Pablo and modified by Silvia

       implicit none

       integer :: i,nresiold,iresi,j,nclosecat,id,ipd,dresi,pa,pb,pg,po,poa,pob,pac,pbc,pgc
       integer :: ifile,ncat,nphos,jcat,jphos,img,ip,itot,icoord,jclosecat,img2
       integer :: nresiphosold,numphosres, numphosresmax
       character*8 :: chainphosold,dchain
       character*9 :: distot,mdistot,pacoord,pbcoord,pgcoord

       real    :: rfix,rdis,rdisp
       real, dimension       (1:2000000) :: rx,ry,rz,rb
       integer, dimension    (1:2000000) :: iatom,nresi
       integer, dimension    (1:2000) ::    icat,iphos     
       character*5, dimension(1:2000000) :: typename
       character*1, dimension(1:2000000) :: disorder
       character*8, dimension(1:2000000) :: namechain,letchain
       character*4, dimension(1:2000000) :: nameres 

       real, dimension (1:3) ::  p0,p1,p2,p3,p4,p5,p6,p7,p8,p9
       real    :: angdeg,dih1deg,dih2deg,dih3deg,dih4deg,dih5deg,distsgc

       character*4 :: pdbname,chardum2
       character*14:: inputfile,checkfile
       character*15:: outputfile,coordfile
       character*8 :: namechainold
       character*6 :: atom
       character*64::charbig,order,chartot,chartot2,charbig2,chartotini
       character*64::chartot3,charbig3
       character*80::fileinp,fileout
       character*1 :: space


       space=' '

       print*,'YOU SHOULD HAVE RUN atom.sh script before'

       open(unit=5,file='pdblist.inp',status='old')
       open(unit=28,file='summarycat.out')
       open(unit=29,file='summarycat_P.out')
       open(unit=27,file='summarycat_P3_OK.out')
       open(unit=41,file='summarycat_P3_OK-disord.out')
       open(unit=26,file='summarycat_P3_disord.out')
       open(unit=24,file='summarycat_P3_coord-le-4.out')
       open(unit=23,file='angle_PA-PB-PG.txt')
       open(unit=22,file='angle_check.txt')
       open(unit=31,file='dihedral_O5-PA-O3A-PB.txt')
       open(unit=33,file='dihedral_PA-O3A-PB-O3B.txt')
       open(unit=34,file='dihedral_O3A-PB-O3B-PG.txt')
       open(unit=32,file='dihedral_check.txt')
       open(unit=36,file='dihedral_OA-PA-PB-OB.txt')
       open(unit=37,file='dihedral_OB-PB-PG-OG.txt')
       open(unit=38,file='dihedral_check_coordinated.txt')
       open(unit=39,file='distance_Mg-O3B.txt')

!      Start the loop to read all the pdb files (SGC)
       do 20 ifile=1,100000
!        Read pdbname from pdblist.inp file (SGC)
         read(5,*,end=21) pdbname
!        Define input, output and temporary files (SGC)
         inputfile='input/'//pdbname//'.dat'
         outputfile='output/'//pdbname//'.out'
         checkfile='check/'//pdbname//'.tmp'
         coordfile='coord/'//pdbname//'.dat'
!        Print on the screen the input file that is opening (SGC)
         print*, 'opening  ', inputfile
!        Open files (SGC)
         open(unit=15,file=inputfile,status='old')
         open(unit= 7,file=outputfile)               
         open(unit=8,file= checkfile   )
         open(unit=10,file= coordfile   )
         open(unit=18,file= 'checkorder.tmp'   )  
         open(unit=9,file='preli.tmp' )
!        Define initial values for the program (SGC)
         nresiold=-1
         namechainold='1234567'
         iresi=0
         rx=0.0d0
         ry=0.0d0
         rz=0.0d0
         iatom=-1
         nresi=-1
         typename='xxxxx'
         disorder='x'
         namechain='xxxxxxxx'
         nameres= 'xxxx'
         ncat=0
         nphos=0
         icat=-1      
         iphos=-1      
         nresiphosold=0
         chainphosold='zzzzzzzz'
         numphosres=0
         numphosresmax=0
!        Start the loop that reads the input file, writes the check file and finds the &
!          cations and P atoms (SGC)
         do i=1,100000000
!          Read each row from input file (.dat file)  (SGC)
!            Definitions of the columns (SGC):
!              atom: word ATOM or HETATOM
!              iatom(i): number of the atom in the pdb file
!              typename(i): type of element (with conectivity, e.g.PG)
!              disorder(i): code the write when there is disorder
!              nameres(i): residue name (e.g. ALA, SER,...)
!              letchain(i): chain (sometimes the asymmetric unit can have more than one)
!              nresi(i): residue number in the protein
!              chardum2: 4 spaces
!              rx(i): cartesian coordinate x
!              ry(i): cartesian coordinate y
!              rz(i): cartesian coordinate z
!              rb(i): occupation number (for disordered residues)
!              rfix(i):???
!              namechain(i): type of element (no conectivity but sometimes charge)
           read(15,200,err=10,end=11) atom,iatom(i),typename(i),disorder(i),           &
            nameres(i),letchain(i), nresi(i),chardum2 ,rx(i),ry(i),rz(i),        &
            rb(i),rfix,namechain(i)
!
!          Every time a metal defined here is found in a row it adds 1 to ncat  (SGC)
!            ncat: number of cations (SGC)
!            icat(ncat): row number for each cation ncat (SGC)
!
!          ADD HERE ALL KIND OF CATIONS...
!
            if  ((namechain(i).eq.'    MG  ')              &
             .or.(namechain(i).eq.'    CA  ')              &
             .or.(namechain(i).eq.'    MN  ')              &
             .or.(namechain(i).eq.'    MG+2')              & 
             .or.(namechain(i).eq.'    MG2+')              &
             .or.(namechain(i).eq.'    MN2+')              &
             .or.(namechain(i).eq.'    CA2+')              &
             .or.(namechain(i).eq.'    CA+2')              &
             .or.(namechain(i).eq.'    MN+2')) then
               ncat=ncat+1             
               icat(ncat)=i
!              It writes in the file fort.555 the metal, the row number, the row number again, 
!                the pdbname and the chain in the PDB   (SGC)
               write(555,*) namechain(i),i,icat(ncat),pdbname,letchain(i)
            endif

!           Every time a P atom is found in a row it adds 1 to nphos  (SGC)
!             nphos: number of P atoms (SGC)
!             iphos(nphos): row number for each P atom (SGC)

            if((namechain(i).eq.'     P  ')                &
             .or.(namechain(i).eq.'     P+4')              &
             .or.(namechain(i).eq.'     P4+')) then
               nphos=nphos+1
               iphos(nphos)=i
!              It writes in the file fort.556 all the row for the P atom   (SGC)
               write(556,200) atom,iatom(i),typename(i),disorder(i), nameres(i),letchain(i),      &
                    nresi(i) ,chardum2,rx(i),ry(i),rz(i),rb(i),rfix,namechain(i)

!              max number of Phosph. in that residue
!              Defined before(SGC): 
!                nresiphosold=0 
!                chainphosold='zzzzzzzz' 
!                numphosres=0
!                numphosresmax=0 
!              If the residue number in the protein (nresi) is equal to nresiphosold and the     &
!                protein chain (letchain(i)) is equal to chainphosold and it adds one to the     &
!                numphosres (number of P atoms in the residue) (SGC)
!              If not it define numphosres=1 (is the first P in that residue) (SGC)
               if((nresi(i).eq.nresiphosold).and.(letchain(i).eq.chainphosold)) then     
                  numphosres=numphosres+1
                  else
                  numphosres=1
               endif
!              If numphosres is greater than numphosresmax it change numphosresmax to numphosres (SGC)
               if(numphosres.gt.numphosresmax) numphosresmax=numphosres 
!              It change the values of nresiphosold and chainphosold to the actual ones (SGC)
               nresiphosold=nresi(i)
               chainphosold=letchain(i)
            endif

!           It writes all the rows it read in the check file (SGC)
            write(8,200) atom,iatom(i),typename(i),disorder(i),nameres(i), letchain(i),      &
              nresi(i) ,chardum2,rx(i),ry(i),rz(i),rb(i),rfix,namechain(i)

!           I don't know what this does or where is the endif (SGC) 
            if(mod(i,10000).eq.0) print*,'line=',i

!        End the loop that reads the input file, writes the check file and finds the &
!          cations and P atoms (SGC)
         enddo

!        If there is an error reading the input file it goes here and print where the error is (SGC)
 10      continue
           print*,'ERROR'
           print*,'in ',inputfile
           print*, 'line=',i
           print*,'     :-(    '
           stop
!        When the input file ends it goes here and print in the screen nphos and ncat (SGC)
 11      continue 

!        It prints on the screen the number of P atoms and cations on thet PDB file (SGC)
         print*,'coord. and MG distances...'
         print*,'n phospho=',nphos
         print*,'n cations=',ncat
       
!        loop over the selected cations...
!        Start the loop over the founded cations (ncat) (SGC)
         do 35 jcat=1,ncat
!          Definition of the parameters and symbols for the output (SGC)
!            icoord: Coordination number of the cation      
           icoord=0
!           nclosecat: number of cation at a distance less than 9.9 (SGC)
           nclosecat=0
           charbig='-->'
           charbig2=' '
           charbig3=' '
           chartot2=charbig2
           chartot3=charbig3
           chartot=charbig
           chartotini=chartot
           distot='no disor'
           mdistot='no disor'
           pbcoord='NO coord'
           pgcoord='NO coord'

!          Parameters for angle calculation (SGC)
           angdeg=0.0d0
           dih1deg=0.0d0
           dih2deg=0.0d0
           dih3deg=0.0d0
           dih4deg=0.0d0
           dih5deg=0.0d0
           distsgc=0.0d0
           pa=0
           pb=0
           pg=0
           po=0
           poa=0
           pob=0
           pac=0
           pbc=0
           pgc=0
           p0(1)=0.0d0
           p0(2)=0.0d0
           p0(3)=0.0d0
           p1(1)=0.0d0
           p1(2)=0.0d0
           p1(3)=0.0d0
           p2(1)=0.0d0
           p2(2)=0.0d0
           p2(3)=0.0d0
           p3(1)=0.0d0
           p3(2)=0.0d0
           p3(3)=0.0d0
           p4(1)=0.0d0
           p4(2)=0.0d0
           p4(3)=0.0d0
           p5(1)=0.0d0
           p5(2)=0.0d0
           p5(3)=0.0d0
           p6(1)=0.0d0
           p6(2)=0.0d0
           p6(3)=0.0d0
           p7(1)=0.0d0
           p7(2)=0.0d0
           p7(3)=0.0d0
           p8(1)=0.0d0
           p8(2)=0.0d0
           p8(3)=0.0d0
           p9(1)=0.0d0
           p9(2)=0.0d0
           p9(3)=0.0d0
           dresi=0
           dchain=' '


!          Remember icat(ncat): row number for each cation ncat, so img is indicating the row number 
!            of the studied cation in the loop (SGC)
           img=icat(jcat)
!          It checks if there is disorder in the cation and if yes it change the status of mdidtot (SGC)
           if (rb(img).lt.1.0d0) then
             mdistot='M disord'
           endif
!          Loop over each cation to see if there are other cations close to them and be able to define
!            if there are 1, 2 or 3 metal ions (SGC)
!          It compares the rx coordinate and if is smaller than 9.9 it calculates the distance between 
!            them (rdis), if it is smaller than 9.9 add 1 to nclosecat (SGC)
           do jclosecat=1,ncat
             img2=icat(jclosecat)
             if(abs(rx(img2)-rx(img)).le.9.9d0) then
               rdis=(rx(img2)-rx(img))**2 + (ry(img2)-ry(img))**2 +(rz(img2)-rz(img))**2 
               rdis=sqrt(rdis)
               if(rdis.le.9.90) then
                 nclosecat=nclosecat+1
               endif
             endif
           enddo

!          Here writes in the output file (SGC)
!            Initial lines, type of element and arrow (SGC)
           write(7,*)
           write(7,*) '++++++++++++++++++++++++++++++++++++++++++++++++', '++++++++++++++++++++++++++++++++++++++++'
           write(7,*)
           write(7,*) typename(img),rb(img),chartot  
!          Writes the Mg coordinates in the coordfile (SGC)
           write(10,*) pdbname,letchain(i)
           write(10,208) namechain(img),rx(img),ry(img),rz(img)
           write(7,*) 

!          Start loop over all the atoms to check if the distance is smaller than 3.0 (SGC)
!            It compares the rx coordinate and continue if it is smaller than 3.0 (SGC)         
           do 30 itot=1,i-1
             if(abs(rx(itot)-rx(img)).le.3.0d0) then
               rdis=(rx(itot)-rx(img))**2 + (ry(itot)-ry(img))**2 +(rz(itot)-rz(img))**2 
               rdis=sqrt(rdis)
!              Check if the element is a P or C atom, if yes it will continue, if not it will check it is not 
!                the same atom and add one to the coordination number (SGC)
               if((namechain(itot).eq.'     P  ')                &
                .or.(namechain(itot).eq.'     C  ')              &
                .or.(namechain(itot).eq.'     P+4')              &
                .or.(namechain(itot).eq.'     P4+')) then
                 continue
                 else
!                Check that the measured distance is smaller than 3.0 and the atom is not itself (SGC)
                 if((rdis.le.3.0d0).and.(itot.ne.img)) then
!                  Add 1 to the coordination number (SGC) 
                   icoord=icoord+1
!                  Write in the output (SGC):
!                    rdis: distance from the metal 
!                    typename(itot) type of element (with conectivity, e.g.PG)
!                    nameres(itot): the residue name 
!                    iatom(itot): number of the atom in the pdb file 
!                    nresi(itot): residue number in the protein
!                    rb(itot): occupation
!                    letchain(itot): chain (A, B, etc)
!                    disorder(itot): disorder code
                   write(7,205) rdis, typename(itot), nameres(itot), iatom(itot),nresi(itot), &
                    rb(itot),letchain(itot),disorder(itot) 
!                  Write the coordinates of the coordinated atom in the coordfile (SGC)
                   write(10,208) namechain(itot),rx(itot),ry(itot),rz(itot)
!                  Check if the coordination number is smaller than 10 and define chartot2 and charbig2 (SGC)
                   if (icoord.le.10) then
!                    chartot2 will be the residues coordinated to the metal in the files summarycat*out (SGC)
!                    It adds here the name of the residue that is coordinated (SGC)
!                    charbig2 is used to storage the residues names and add the next one (SGC)        
                     chartot2=TRIM(charbig2)//space//nameres(itot)
                     charbig2=chartot2
                     else
!                    If more than 10 atoms are coordinated then it changes the residues names by: (SGC)
                     chartot2='XXXXXXX more than 10 XXXXXXXXXX'
                     charbig2=chartot2
                   endif
!                  If there is disorder in the coordinated atom change distot status
                   if (rb(itot).lt.1.0) then 
                     distot='disorder'
                   endif
!                  Loop to look is the atom is connected to a P atom (looking for the coordinated O atoms &
!                    conected to the P atom) (SGC)
!                  the distance between the coordinated atom and the P (rdisp) should be less than 2.0 (SGC)
                   do 25 jphos=1,nphos
                     ip=iphos(jphos)
                     if((abs(rx(ip)-rx(itot)).le.2).and.(nresi(itot).eq.nresi(ip)) &
                      .and.(disorder(itot).eq.disorder(ip))) then
                        rdisp=(rx(ip)-rx(itot))**2 + (ry(ip)-ry(itot))**2 +(rz(ip)-rz(itot))**2 
                        rdisp=sqrt(rdisp)
                        if(rdisp.le.2) then
!                       If there is an atom with a distance smaller that 2 the program defines chartot and charbig 
!                         adding the typename of the P atomat which the coordinated atom is connected (SGC)
                          chartot=TRIM(charbig)//space//space//typename(ip)
                          charbig=chartot
!                         Define and store the residue number and chain of this residue (NTP) (SGC)
                          dresi=nresi(ip)
                          dchain=letchain(ip)
! --- Selection of the coordinated Oxigens bond to P atoms of the NTP ---
!                          If the O is bonded to the PA atom then define p7(O) (the 3D points for the dihedral) (SGC)
                           if ((typename(ip).eq.' PA ')) then
                             p7(1)=rx(itot)
                             p7(2)=ry(itot)
                             p7(3)=rz(itot)
!                            Defines chartot3 and charbig3 writing the coordinated O and the P atom bonded to it (SGC)
                             chartot3=TRIM(charbig3)//space//typename(itot)//space//typename(ip)
                             charbig3=chartot3
!                            Changes the status of pbcoord to PB coord (SGC)
                             pacoord='PA coord'
!                            Define pbc to know how many coordinated O are bonded to the PB (SGC)
                             pac=pac+1
                           endif
!                          If the O is bonded to the PB atom then define p8(O) (the 3D points for the dihedral) (SGC)
                           if ((typename(ip).eq.' PB ')) then
                             p8(1)=rx(itot)
                             p8(2)=ry(itot)
                             p8(3)=rz(itot)
!                            Defines chartot3 and charbig3 writing the coordinated O and the P atom bonded to it (SGC)
                             chartot3=TRIM(charbig3)//space//typename(itot)//space//typename(ip)
                             charbig3=chartot3
!                            Changes the status of pbcoord to PB coord (SGC)
                             pbcoord='PB coord'
!                            Define pb to know how many coordinated O are bonded to the PB (SGC)
                             pbc=pbc+1
                           endif
!                          If the O is bonded to the PB atom then define p4(O) and p3(P) (the 3D points for the dihedral)
                           if ((typename(ip).eq.' PG ')) then
                             p9(1)=rx(itot)
                             p9(2)=ry(itot)
                             p9(3)=rz(itot)
!                            Defines chartot3 and charbig3 writing the coordinated O and the P atom bonded to it (SGC)
                             chartot3=TRIM(charbig3)//space//typename(ip)//space//typename(itot)
                             charbig3=chartot3
!                            Changes the status of pgcoord to PG coord (SGC)
                             pgcoord='PG coord'
!                            Define pg to know how many coordinated O are bonded to the PG (SGC)
                             pgc=pgc+1
                           endif
! --- End Selection of the coordinated Oxigens bond to P atoms of the NTP ---

                        endif
                     endif
!                  End loop to look is the atom is connected to a P atom (SGC)
  25               continue
!                End checking that the measured distance is smaller than 3.0 and the atom is not itself (SGC)
                 endif
!              end checking if the element is a P or C atom (SGC)
               endif
!            end comparing the rx coordinate  (SGC) 
             endif
!          End loop over all the atoms to check if the distance is smaller than 3.0 (SGC)
  30       continue

! ------ Angle and Dihedral selection of atoms and calculation -------

!          Loop over all the atoms to define the points we will use for the angle calculation (SGC)
!           The angle is PA-PB-PG we will use the coordinates of these 3 atoms
!           in the NTP defined before (dresi) in the same chain that the metal cation (SGC) 
           do 40 id=1,i-1
             if ((nresi(id).eq.dresi).and.(typename(id).eq." PA ") &
              .and.(letchain(id).eq.letchain(img))) then
                p1(1)=rx(id)
                p1(2)=ry(id)
                p1(3)=rz(id)
                pa=pa+1
             endif
             if ((nresi(id).eq.dresi).and.(typename(id).eq." PB ") &
              .and.(letchain(id).eq.letchain(img))) then
                p2(1)=rx(id)
                p2(2)=ry(id)
                p2(3)=rz(id)
                pb=pb+1
             endif
             if ((nresi(id).eq.dresi).and.(typename(id).eq." PG ") &
              .and.(letchain(id).eq.letchain(img))) then
                p3(1)=rx(id)
                p3(2)=ry(id)
                p3(3)=rz(id)
                pg=pg+1
             endif
             if ((nresi(id).eq.dresi).and.(typename(id).eq." O5'") &
              .and.(letchain(id).eq.letchain(img))) then
                p4(1)=rx(id)
                p4(2)=ry(id)
                p4(3)=rz(id)
                po=po+1
             endif
!            p5 and p6 are defined in 3 different ways because in some of the cases in the NTP the O3A is substituted
!             by N or C to crystalize it (SGC)
             if ((nresi(id).eq.dresi).and.(typename(id).eq." O3A") &
              .and.(letchain(id).eq.letchain(img))) then
                p5(1)=rx(id)
                p5(2)=ry(id)
                p5(3)=rz(id)
                poa=poa+1
             endif
             if ((nresi(id).eq.dresi).and.(typename(id).eq." C3A") &
              .and.(letchain(id).eq.letchain(img))) then
                p5(1)=rx(id)
                p5(2)=ry(id)
                p5(3)=rz(id)
                poa=poa+1
             endif
             if ((nresi(id).eq.dresi).and.(typename(id).eq." N3A") &
              .and.(letchain(id).eq.letchain(img))) then
                p5(1)=rx(id)
                p5(2)=ry(id)
                p5(3)=rz(id)
                poa=poa+1
             endif
             if ((nresi(id).eq.dresi).and.(typename(id).eq." O3B") &
              .and.(letchain(id).eq.letchain(img))) then
                p6(1)=rx(id)
                p6(2)=ry(id)
                p6(3)=rz(id)
                pob=pob+1
             endif
             if ((nresi(id).eq.dresi).and.(typename(id).eq." C3B") &
              .and.(letchain(id).eq.letchain(img))) then
                p6(1)=rx(id)
                p6(2)=ry(id)
                p6(3)=rz(id)
                pob=pob+1
             endif
             if ((nresi(id).eq.dresi).and.(typename(id).eq." N3B") &
              .and.(letchain(id).eq.letchain(img))) then
                p6(1)=rx(id)
                p6(2)=ry(id)
                p6(3)=rz(id)
                pob=pob+1
             endif

!          End loop over all the atoms to define the points we will use for the angle calculation (SGC)
!          pa,pb,pg,po,poa,pob will tell us if the NTP is OK or there is something wrong with the structure and we need
!           check it manually (SGC)
  40       continue

!          ANGLE calculation: PA-PB-PG angle (SGC)
!          The calculation of the angle will be done using the acos function (SGC)
!          It will do it only if pa, pb and pg are equal to 1 (SGC) 
           if ((pa.eq.1).and.(pb.eq.1).and.(pg.eq.1)) then
             call obtain_angle(p1,p2,p3,angdeg)
           endif
!          End ANGLE calculation: PA-PB-PG angle (SGC)

!          DIHEDRAL ANGLES calculation: (SGC)

!      O5-PA-O3A-PB angle (SGC)
!           See how to calculate dihedral angles in the word document Dihedral_angle_calculation_info (SGC)
!           It will do it only if po, pa, poa and pb are equal to 1 (SGC) 
           if ((po.eq.1).and.(pa.eq.1).and.(poa.eq.1).and.(pb.eq.1)) then
             call obtain_dihedral(p4,p1,p5,p2,dih1deg)
           endif
!      O5-PA-O3A-PB angle (SGC)

!      PA-O3A-PB-O3B angle (SGC)
!           See how to calculate dihedral angles in the word document Dihedral_angle_calculation_info (SGC)
!           It will do it only if pa, poa, pb and pob are equal to 1 (SGC) 
           if ((pa.eq.1).and.(poa.eq.1).and.(pb.eq.1).and.(pob.eq.1)) then
             call obtain_dihedral(p1,p5,p2,p6,dih2deg)
           endif
!      PA-O3A-PB-O3B angle (SGC)

!      O3A-PB-O3B-PG angle (SGC)
!           See how to calculate dihedral angles in the word document Dihedral_angle_calculation_info (SGC)
!           It will do it only if pa, poa, pb and pob are equal to 1 (SGC) 
           if ((poa.eq.1).and.(pb.eq.1).and.(pob.eq.1).and.(pg.eq.1)) then
             call obtain_dihedral(p5,p2,p6,p3,dih3deg)
           endif
!      O3A-PB-O3B-PG angle (SGC)

!      OA-PA-PB-OB angle (SGC)
!           See how to calculate dihedral angles in the word document Dihedral_angle_calculation_info (SGC)
!           It will do it only if PA and PB are "coordinated" (SGC) 
           if ((pacoord.eq.'PA coord').and.(pbcoord.eq.'PB coord')) then
             call obtain_dihedral(p7,p1,p2,p8,dih4deg)
           endif
!      OA-PA-PB-OB angle (SGC)

!      OB-PB-PG-OG angle (SGC)
!           See how to calculate dihedral angles in the word document Dihedral_angle_calculation_info (SGC)
!           It will do it only if PB and PG are "coordinated" (SGC) 
           if ((pbcoord.eq.'PB coord').and.(pgcoord.eq.'PG coord')) then
             call obtain_dihedral(p8,p2,p3,p9,dih5deg)
           endif
!      OB-PB-PG-OG angle (SGC)


!          End DIHEDRAL ANGLES calculation (SGC)

! ------ End Angle and Dihedral selection of atoms and calculation -------

! ------ Obtain distance between O3B atom and Mg ion
!       To easily differenciate between 2 coordinations that are very
!       similar
          p0(1)=rx(img)
          p0(2)=ry(img)
          p0(3)=rz(img)
          call obtain_distance(p6,p0,distsgc)

! ------ Write angles and dihedrals files ------

!          Write the angle in the angle.out file (SGC)
           if((nphos.gt.0).and.(numphosresmax.ge.3).and.(icoord.gt.4)   .and.(chartot.ne.chartotini) &
            .and.(mdistot.eq.'no disor').and.(distot.eq.'no disor').and. &
            (pa.eq.1).and.(pb.eq.1).and.(pg.eq.1)) &
            write(23,223) pdbname,typename(img),nresi(img),letchain(img),dresi,dchain,"angle = ",angdeg
!          Write the pbdnames of the angle to check, if pa,pb or pg are not equal to 1 (SGC)
           if((nphos.gt.0).and.(numphosresmax.ge.3).and.(icoord.gt.4)   .and.(chartot.ne.chartotini) &
            .and.(mdistot.eq.'no disor').and.(distot.eq.'no disor').and. &
            ((pa.ne.1).or.(pb.ne.1).or.(pg.ne.1))) &
            write(22,224) pdbname,typename(img),nresi(img),letchain(img),dresi,dchain,pa,pb,pg,"angle = ",angdeg

!          Write the angle in the dihedral_O5-PA-O3A-PB.out file (SGC)
           if((nphos.gt.0).and.(numphosresmax.ge.3).and.(icoord.gt.4)   .and.(chartot.ne.chartotini) &
            .and.(mdistot.eq.'no disor').and.(distot.eq.'no disor').and. &
            (po.eq.1).and.(poa.eq.1).and.(pa.eq.1).and.(pb.eq.1)) &
            write(31,223) pdbname,typename(img),nresi(img),letchain(img),dresi,dchain,"angle = ",dih1deg
!          Write the angle in the dihedral_PA-O3A-PB-O3B.out file (SGC)
           if((nphos.gt.0).and.(numphosresmax.ge.3).and.(icoord.gt.4)   .and.(chartot.ne.chartotini) &
            .and.(mdistot.eq.'no disor').and.(distot.eq.'no disor').and. &
            (pa.eq.1).and.(poa.eq.1).and.(pb.eq.1).and.(pob.eq.1)) &
            write(33,223) pdbname,typename(img),nresi(img),letchain(img),dresi,dchain,"angle = ",dih2deg
!          Write the angle in the dihedral_O3A-PB-O3B-PG.out file (SGC)
           if((nphos.gt.0).and.(numphosresmax.ge.3).and.(icoord.gt.4)   .and.(chartot.ne.chartotini) &
            .and.(mdistot.eq.'no disor').and.(distot.eq.'no disor').and. &
            (poa.eq.1).and.(pb.eq.1).and.(pob.eq.1).and.(pg.eq.1)) &
            write(34,223) pdbname,typename(img),nresi(img),letchain(img),dresi,dchain,"angle = ",dih3deg

!          Write the pbdnames of the dihedral to check, if po,pa,poa or pb are not equal to 1 (SGC)
           if((nphos.gt.0).and.(numphosresmax.ge.3).and.(icoord.gt.4)   .and.(chartot.ne.chartotini) &
            .and.(mdistot.eq.'no disor').and.(distot.eq.'no disor').and. &
            ((po.ne.1).or.(pa.ne.1).or.(poa.ne.1).or.(pb.ne.1).or.      (pob.ne.1).or.(pg.ne.1))) &
            write(32,225) pdbname,typename(img),nresi(img),letchain(img),dresi,dchain,po,pa,poa,pb,pob,pg,"angles = ", &
            dih1deg,dih2deg,dih3deg

!          Write the angle in the dihedral_OA-PA-PB-OB.out file
           if((nphos.gt.0).and.(numphosresmax.ge.3).and.(icoord.gt.4)   .and.(chartot.ne.chartotini) &
            .and.(mdistot.eq.'no disor').and.(distot.eq.'no disor').and. &
            (pacoord.eq.'PA coord').and.(pbcoord.eq.'PB coord')    &
            .and.(pac.eq.1).and.(pbc.eq.1)) &
            write(36,226) pdbname,typename(img),nresi(img),letchain(img),dresi,dchain,chartot3,"angle = ",dih4deg
!          Write the angle in the dihedral_OB-PB-PG-OG.out file
           if((nphos.gt.0).and.(numphosresmax.ge.3).and.(icoord.gt.4)   .and.(chartot.ne.chartotini) &
            .and.(mdistot.eq.'no disor').and.(distot.eq.'no disor').and. &
            (pbcoord.eq.'PB coord').and.(pgcoord.eq.'PG coord')    &
            .and.(pbc.eq.1).and.(pgc.eq.1)) &
            write(37,226) pdbname,typename(img),nresi(img),letchain(img),dresi,dchain,chartot3,"angle = ",dih5deg
!          Write the pbdnames of the dihedral to check
           if((nphos.gt.0).and.(numphosresmax.ge.3).and.(icoord.gt.4)   .and.(chartot.ne.chartotini) &
            .and.(mdistot.eq.'no disor').and.(distot.eq.'no disor').and. &
            (pbcoord.eq.'PB coord').and.((pacoord.eq.'PA coord').or.(pgcoord.eq.'PG coord'))    &
            .and.((pac.ne.1).or.(pbc.ne.1).or.(pgc.ne.1))) &
            write(38,227) pdbname,typename(img),nresi(img),letchain(img),dresi,dchain,chartot3,"angles = ",dih4deg,dih5deg

! ------ End Write angles and dihedrals files ------

!          Write distance between the Metal ion and the O3B atom
           if((nphos.gt.0).and.(numphosresmax.ge.3).and.(icoord.gt.4)   .and.(chartot.ne.chartotini) &
            .and.(mdistot.eq.'no disor').and.(distot.eq.'no disor')) &
            write(39,228) pdbname,typename(img),nresi(img),letchain(img),dresi,dchain,"distance = ",distsgc


! ------ Write summary files ------

!          Now it writes in the summarycat files the pdbname, type of cation and if it is distorted (occupation lt 1) (SGC)
           write(28,206,advance="no") pdbname,typename(img),iatom(img),letchain(img),mdistot,numphosresmax 
           if(nphos.gt.0) write(29,206,advance="no") pdbname,  typename(img),iatom(img),letchain(img),mdistot,numphosresmax
           if((nphos.gt.0).and.(numphosresmax.ge.3).and.(icoord.gt.4)   .and.(chartot.ne.chartotini) &
            .and.(mdistot.eq.'no disor').and.(distot.eq.'no disor')) &
            write(27,207,advance="no") pdbname,typename(img),letchain(img),mdistot
           if((nphos.gt.0).and.(numphosresmax.ge.3).and.(icoord.gt.4)   .and.(chartot.ne.chartotini)) &
            write(41,207,advance="no") pdbname,typename(img),letchain(img),mdistot
           if((nphos.gt.0).and.(numphosresmax.ge.3).and.(icoord.gt.4)   .and.(chartot.ne.chartotini) &
            .and.((mdistot.eq.'M disord').or.(distot.eq.'disorder'))) &
            write(26,207,advance="no") pdbname,typename(img),letchain(img),mdistot
           if((nphos.gt.0).and.(numphosresmax.ge.3).and.(icoord.lt.4)   .and.(chartot.ne.chartotini)) &
            write(24,207,advance="no") pdbname,typename(img),letchain(img),mdistot


!          Write in the summarycat files the coordinated P atoms, residues and if the coordinated atoms are distorted  (SGC)
           write(28,210,advance="no") chartot, chartot2, distot
           if(nphos.gt.0)  write(29,210,advance="no") chartot, chartot2, distot
           if((nphos.gt.0).and.(numphosresmax.ge.3).and.(icoord.gt.4)   .and.(chartot.ne.chartotini) &
            .and.(mdistot.eq.'no disor').and.(distot.eq.'no disor')) &
            write(27,210,advance="no") chartot, chartot2, distot
           if((nphos.gt.0).and.(numphosresmax.ge.3).and.(icoord.gt.4)   .and.(chartot.ne.chartotini)) &
            write(41,210,advance="no") chartot, chartot2, distot
           if((nphos.gt.0).and.(numphosresmax.ge.3).and.(icoord.gt.4)   .and.(chartot.ne.chartotini) &
            .and.((mdistot.eq.'M disord').or.(distot.eq.'disorder'))) &
            write(26,210,advance="no") chartot, chartot2, distot
           if((nphos.gt.0).and.(numphosresmax.ge.3).and.(icoord.lt.4)   .and.(chartot.ne.chartotini)) &
            write(24,210,advance="no") chartot, chartot2, distot

!          Write in the summarycat files the coordination number and number of metal ions (SGC)
           write(28,211) "  ncoord= ",icoord,nclosecat," metal ions"
           if(nphos.gt.0) write(29,211) "  ncoord= ",icoord,nclosecat,  " metal ions"    
           if((nphos.gt.0).and.(numphosresmax.ge.3).and.(icoord.gt.4)   .and.(chartot.ne.chartotini) &
              .and.(mdistot.eq.'no disor').and.(distot.eq.'no disor'))     &
              write(27,211)"  ncoord= ",icoord,nclosecat," metal ions"
           if((nphos.gt.0).and.(numphosresmax.ge.3).and.(icoord.gt.4)   .and.(chartot.ne.chartotini)) &
              write(41,211)"  ncoord= ",icoord,nclosecat," metal ions"
           if((nphos.gt.0).and.(numphosresmax.ge.3).and.(icoord.gt.4)   .and.(chartot.ne.chartotini) &
              .and.((mdistot.eq.'M disord').or.(distot.eq.'disorder')))     &
              write(26,211)"  ncoord= ",icoord,nclosecat," metal ions"
           if((nphos.gt.0).and.(numphosresmax.ge.3).and.(icoord.lt.4)   .and.(chartot.ne.chartotini)) &
              write(24,211)"  ncoord= ",icoord,nclosecat," metal ions"

!          Write in the output the coordination number and "coordinated P" (SGC)
           write(7,*)
           write(7,*)
           write(7,*) '     ncoord  =', icoord
           write(7,*)
           write(7,*) '     nphos coord', chartot     
           write(7,*)
!        End the loop over the founded cations (ncat) (SGC)
  35     continue
!        Print done for that PDB file and close the files for that PDB
         print*,'DONE FOR FILE', ifile 
         close(unit=15)
         close(unit=10)
         close(unit=7)
         close(unit=8)
         write(18,*) 'diff ',inputfile,' ', checkfile
!      End the loop to read all the pdb files (SGC)
  20   continue
!      It is done!!!
  21   print*,'done for every input file :-D'
  12   print*,'TACHAN!!!'

!      FORMAT for the different files (SGC):
!        Format 200 - check file (SGC)
 200   format(a6,i5,1x,a4,a1,a3,1x,a1,i4,a2,2x,f8.3,f8.3,f8.3,f6.2,f6.2,6x,a8)
!        Format 205 - output file (SGC)
 205   format(5x,f8.2,2x,a5,2x,a4,2x,i9,2x,i8,2x,f6.2,2x,a8,2x,a1)
!        Format 206 - summarycat* files 4 first columns (SGC)
 206   format(2x,a4,2x,a5,i5,1x,a1,1x,a10,2x,i3) 
!        Format 207 - summarycat* files 4 first columns (SGC)
 207   format(2x,a4,2x,a5,1x,a1,1x,a10) 
!        Format 208 - coord file (SGC)
 208   format(a8,2x,f8.3,f8.3,f8.3)
!        Format 210 - summarycat* files chartot section (SGC)
 210   format(2x,a30,2x,a50,2x,a10)
!        Format 211 - summarycat* files final section (SGC)
 211   format(a11,2x,i3,2x,i3,2x,a13)
!        Format 223 - dihedral-angle.out file final section (SGC)
 223   format(2x,a4,2x,a5,2x,i4,2x,a8,i4,2x,a8,a8,f8.3)
!        Format 224 - angle_check.out file final section (SGC)
 224   format(2x,a4,2x,a5,2x,i4,2x,a8,i4,2x,a8,i2,1x,i2,1x,i2,2x,a8,f8.3)
!        Format 225 - dihedral_check.out file final section (SGC)
 225   format(2x,a4,2x,a5,2x,i4,2x,a8,i4,2x,a8,i2,1x,i2,1x,i2,1x,i2,2x,i2,2x,i2,2x,a9,f8.3,f8.3,f8.3)
!        Format 226 - dihedral_coordinated.out file final section (SGC)
 226   format(2x,a4,2x,a5,2x,i4,2x,a8,i4,2x,a8,a50,2x,a8,f8.3)
!        Format 227 - dihedral_check_coordinated.out file final section (SGC)
 227   format(2x,a4,2x,a5,2x,i4,2x,a8,i4,2x,a8,a50,2x,a9,f8.3,f8.3)
!        Format 228 - distance_M-O3B.out file final section (SGC)
 228   format(2x,a4,2x,a5,2x,i4,2x,a8,i4,2x,a8,2x,a11,f8.3)

       end program writepdb

       subroutine obtain_vector(a,b,v)
!        It calculates the vector v from the sustraction of the xyz coordinates of points a and b 
         implicit none
         real, dimension(3), intent(in) :: a,b
         real, dimension(3), intent(out) :: v
         integer :: icol
           do icol=1,3
             v(icol)=b(icol)-a(icol)
           enddo
       end subroutine obtain_vector

       subroutine obtain_distance(a,b,d)
!        It calculates the distance between two points
         implicit none
         real, dimension(3), intent(in) :: a,b
         real, intent(out) :: d
         d=sqrt((b(1)-a(1))**2+(b(2)-a(2))**2+(b(3)-a(3))**2)
       end subroutine obtain_distance

       subroutine normalize_vector(a,n)
!        It normalize the vector a and give it you the vector n
         implicit none
         real, dimension(3), intent(in) :: a
         real, dimension(3), intent(out) :: n
         real :: md
         md=sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
         n(1)=a(1)/md
         n(2)=a(2)/md
         n(3)=a(3)/md
       end subroutine normalize_vector

       subroutine dot_prod(a,b,x)
!        It does the dot product of vectors a and b and gives you the value x
         implicit none
         real, dimension(3), intent(in) :: a,b
         real, intent(out) :: x
         x=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
       end subroutine dot_prod


       subroutine cross_prod(a,b,v)
!        It does the cross product of vectors a and b and give you the new vector v
!          which is perpendicular to a and b
         implicit none
         real, dimension(3), intent(in) :: a,b
         real, dimension(3), intent(out) :: v
         v(1)=(a(2)*b(3)-b(2)*a(3))
         v(2)=(-a(1)*b(3)+b(1)*a(3))
         v(3)=(a(1)*b(2)-b(1)*a(2))
       end subroutine cross_prod

       subroutine obtain_angle(a,b,c,x)
         implicit none
         real, dimension(3), intent(in) :: a,b,c
         real, intent(out) :: x
         real, dimension (1:3) ::  v1,v2,nv1,nv2
         real    :: z,ang
!        Calculate vectors between the 3 atoms (SGC)
!        PB will be used as center
         call obtain_vector(b,a,v1)
         call obtain_vector(b,c,v2)
!        Normalize the calculated vectors (SGC)
         call normalize_vector(v1,nv1)
         call normalize_vector(v2,nv2)
!        Calculate the dot product of the 2 normalized vectors
         call dot_prod(nv1,nv2,z)
!        Calculate angle using acos function (SGC)
         ang=acos(z)
!        Change the angle from radians to degrees (SGC)
         x=ang*(180/3.1415927)
       end subroutine obtain_angle

       subroutine obtain_dihedral(a,b,c,d,z)
         implicit none
         real, dimension(3), intent(in) :: a,b,c,d
         real, intent(out) :: z
         real, dimension (1:3) ::  v1,v2,v3,nv1,nv2,nv3,npl1,npl2,m1
         real    :: x,y,dih
!        Calculate vectors between the 4 atoms (SGC)
         call obtain_vector(a,b,v1)
         call obtain_vector(b,c,v2)
         call obtain_vector(c,d,v3)
!        Normalize the calculated vectors (SGC)
         call normalize_vector(v1,nv1)
         call normalize_vector(v2,nv2)
         call normalize_vector(v3,nv3)
!        Calculate the perpendicular vector to the plane composed by 2 vectors (SGC)
         call cross_prod(nv1,nv2,npl1)
         call cross_prod(nv2,nv3,npl2)
!        Calculate the extra reference vector m1 (SGC)
         call cross_prod(npl1,nv2,m1)
!        Calculate x and y (SGC)
         call dot_prod(npl1,npl2,x)
         call dot_prod(m1,npl2,y)
!        Calculate angle using atan2 function (SGC)
         dih=atan2(y,x)
!        Change the angle from radians to degrees and change the sign because the atan2 function
!         define the angle with the opposite sign than the IUPAC (SGC)
         z=dih*(-180/3.1415927)
       end subroutine obtain_dihedral



