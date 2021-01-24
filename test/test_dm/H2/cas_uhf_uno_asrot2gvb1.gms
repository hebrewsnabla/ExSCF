----- GAMESS execution script 'rungms' -----
This job is running on host xn01
under operating system Linux at Wed Jan 20 23:21:34 CST 2021
Available scratch disk space (Kbyte units) at beginning of the job is
Filesystem     1K-blocks  Used Available Use% Mounted on
/dev/sda5      150914052 54728 150859324   1% /scratch
GAMESS temporary binary files will be written to /scratch/scr/srwang/gamess
GAMESS supplementary output files will be written to /scratch/scr/srwang/gamess
Copying input file ./cas_uhf_uno_asrot2gvb1.inp to your run's scratch directory...
cp ./cas_uhf_uno_asrot2gvb1.inp /scratch/scr/srwang/gamess/cas_uhf_uno_asrot2gvb1.F05
unset echo
/share/home/srwang/gamess/ddikick.x /share/home/srwang/gamess/gamess.01.x cas_uhf_uno_asrot2gvb1 -ddi 1 1 xn01 -scr /scratch/scr/srwang/gamess

 Distributed Data Interface kickoff program.
 Initiating 1 compute processes on 1 nodes to run the following command:
 /share/home/srwang/gamess/gamess.01.x cas_uhf_uno_asrot2gvb1 

          *******************************************************
          *                                                     *
          *          GAMESS VERSION = 30 JUN 2020 (R1)          *
          *                                                     *
          *              doi.org/10.1063/5.0005188              *
          *                                                     *
          **************** 64 BIT INTEL VERSION *****************

  GAMESS HAS BEEN MADE POSSIBLE WITH IMPORTANT CONTRIBUTIONS FROM
  THE FOLLOWING INDIVIDUALS (IN ALPHABETICAL ORDER):

     IVANA ADAMOVIC, CHRISTINE AIKENS, TOMOKO AKAMA,
     YURI ALEXEEV, YURIKO AOKI, POOJA ARORA, TOSHIO ASADA,
     ANDREY ASADCHEV, KIM K. BALDRIDGE, PRADIPTA BANDYOPADHYAY,
     MARIA BARYSZ, ROB BELL, JONATHAN BENTZ, COLLEEN BERTONI,
     JERRY A. BOATZ, BRETT BODE, KURT BRORSEN, KURT R. BRORSEN,
     LAIMUTIS BYTAUTAS, CALEB CARLIN, LAURA CARRINGTON,
     GALINA CHABAN, BENOIT CHAMPAGNE, WEI CHEN, MAHITO CHIBA,
     DAN CHIPMAN, CHEOL HO CHOI, TANNER CULPITT, PAUL DAY,
     ALBERT DEFUSCO, NUWAN DESILVA, J. EMILIANO DEUSTUA,
     TIM DUDLEY, MICHEL DUPUIS, STEVEN T. ELBERT,
     DMITRI FEDOROV, ALEX FINDLATER, GRAHAM FLETCHER,
     MARK FREITAG, CHRISTIAN FRIEDL, DAVID GARMER,
     IGOR S. GERASIMOV, KURT GLAESEMANN, MARK S. GORDON,
     JEFFREY GOUR, FENG LONG GU, EMILIE B. GUIDEZ,
     ANASTASIA GUNINA, SHARON HAMMES-SCHIFFER, MASATAKE HARUTA,
     KIMIHIKO HIRAO, YASUHIRO IKABATA, TZVETELIN IORDANOV,
     STEPHAN IRLE, KAZUYA ISHIMURA, JOE IVANIC, FRANK JENSEN,
     JAN H. JENSEN, VISVALDAS KAIRYS, MUNEAKI KAMIYA,
     MICHIO KATOUDA, NAOAKI KAWAKAMI, DAN KEMP, BERNARD KIRTMAN,
     KAZUO KITAURA, MARIUSZ KLOBUKOWSKI, MASATO KOBAYASHI,
     PRAKASHAN KORAMBATH, JACEK KORCHOWIEC, SHIRO KOSEKI,
     KAROL KOWALSKI, JIMMY KROMANN, STANISLAW KUCHARSKI,
     HENRY KURTZ, SAROM SOK LEANG, HUI LI, SHUHUA LI, WEI LI,
     JESSE LUTZ, ALEKSANDR O. LYKHIN, MARCIN MAKOWSKI,
     JOANI MATO, NIKITA MATSUNAGA, BENEDETTA MENNUCCI,
     GRANT MERRILL, NORIYUKI MINEZAWA, VLADIMIR MIRONOV,
     EISAKU MIYOSHI, JOHN A. MONTGOMERY JR., HIROTOSHI MORI,
     JONATHAN MULLIN, MONIKA MUSIAL, SHIGERU NAGASE,
     TAKESHI NAGATA, HIROMI NAKAI, TAKAHITO NAKAJIMA,
     YUYA NAKAJIMA, HARUYUKI NAKANO, HIROYA NAKATA, SEAN NEDD,
     HEATHER NETZLOFF, KIET A. NGUYEN, YOSHIO NISHIMOTO,
     BOSILJKA NJEGIC, RYAN OLSON, MIKE PAK, ROBERTO PEVERATI,
     BUU PHAM, PIOTR PIECUCH, ANNA POMOGAEVA, DAVID POOLE,
     SPENCER PRUITT, OLIVIER QUINET, LUKE ROSKOP,
     KLAUS RUEDENBERG, ANDREW SAND, TOSAPORN SATTASATHUCHANA,
     NOZOMI SAWADA, MICHAEL W. SCHMIDT, PATRICK E. SCHNEIDER,
     JUNJI SEINO, PRACHI SHARMA, JUN SHEN, JIM SHOEMAKER,
     YINAN SHU, DEJUN SI, JONATHAN SKONE, LYUDMILA SLIPCHENKO,
     TONY SMITH, JIE SONG, MARK SPACKMAN, CASPER STEINMANN,
     WALT STEVENS, PEIFENG SU, SHUJUN SU, CHET SWALINA,
     TETSUYA TAKETSUGU, ZHEN TAO, NANDUN THELLAMUREGE
     SEIKEN TOKURA, JACOPO TOMASI, TSUGUKI TOUMA, TAKAO TSUNEDA,
     HIROAKI UMEDA, JORGE LUIS GALVEZ VALLEJO, YALI WANG,
     SIMON WEBB, AARON WEST, THERESA L. WINDUS, MARTA WLOCH,
     PENG XU, KIYOSHI YAGI, SUSUMU YANAGISAWA, YANG YANG,
     SOOHAENG YOO, TAKESHI YOSHIKAWA, FEDERICO ZAHARIEV,
     TOBY ZENG

  WHO ARE SUPPORTED BY THEIR INSTITUTION/UNIVERSITY/COMPANY/GROUP
  (IN ALPHABETICAL ORDER):

     EP ANALYTICS, FACULTES UNIVERSITAIRES NOTRE-DAME DE LA PAIX,
     INSTITUTE FOR MOLECULAR SCIENCE, IOWA STATE UNIVERSITY,
     JACKSON STATE UNIVERSITY, JOHANNES KEPLER UNIVERSITY LINZ,
     KYUSHU UNIVERSITY, MICHIGAN STATE UNIVERSITY,
     MIE UNIVERSITY, MOSCOW STATE UNIVERSITY,
     N. COPERNICUS UNIVERSITY, NAGOYA UNIVERSITY,
     NANJING UNIVERSITY, NAT. INST. OF ADVANCED INDUSTRIAL
     SCIENCE AND TECHNOLOGY, NATIONAL INST. OF STANDARDS AND
     TECHNOLOGY, NESMEYANOV INSTITUTE OF ORGANOELEMENT COMPOUNDS
     OF RUSSIAN ACADEMY OF SCIENCES, OSAKA PREFECTURE UNIVERSITY,
     PENNSYLVANIA STATE UNIVERSITY, TOKYO INSTITUTE OF TECHNOLOGY,
     UNIVERSITY OF AARHUS, UNIVERSITY OF ALBERTA,
     UNIVERSITY OF CALIFORNIA AT SANTA BARBARA,
     UNIVERSITY OF COPENHAGEN, UNIVERSITY OF IOWA,
     UNIVERSITY OF MEMPHIS, UNIVERSITY OF MINNESOTA,
     UNIVERSITY OF NEBRASKA, UNIVERSITY OF NEW ENGLAND,
     UNIVERSITY OF NOTRE DAME, UNIVERSITY OF PISA,
     UNIVERSITY OF SILESIA, UNIVERSITY OF TOKYO,
     UNIVERSITY OF ZURICH, WASEDA UNIVERSITY, YALE UNIVERSITY

  GAMESS SOFTWARE MANAGEMENT TEAM FOR THIS RELEASE:

     MARK S. GORDON (TEAM LEAD),
     BRETT BODE (SENIOR ADVISOR),
     GIUSEPPE BARCA,
     COLLEEN BERTONI,
     KRISTOPHER KEIPERT (ADVISOR),
     SAROM S. LEANG (DEVELOPMENT LEAD),
     BUU PHAM,
     JORGE LUIS GALVEZ VALLEJO (WEBSITE ADMINISTRATOR),
     PENG XU

 EXECUTION OF GAMESS BEGUN Wed Jan 20 23:21:34 2021

            ECHO OF THE FIRST FEW INPUT CARDS -
 INPUT CARD> $CONTRL SCFTYP=GVB RUNTYP=ENERGY ICHARG=0 MULT=1 NOSYM=1 ICUT=11               
 INPUT CARD>  MAXIT=500 RELWFN=DK $END                                                      
 INPUT CARD> $SYSTEM MWORDS=500 $END                                                        
 INPUT CARD> $SCF NCO=0 NPAIR=1 DIRSCF=.TRUE. $END                                          
 INPUT CARD> $GUESS GUESS=MOREAD NORB=4 $END                                                
 INPUT CARD> $DATA                                                                          
 INPUT CARD>GAMESS inp format file produced by MOKIT, nbf=4                                 
 INPUT CARD>C1   1                                                                          
 INPUT CARD>H     1.         0.00000000         0.00000000         0.00000000               
 INPUT CARD>   S  2                                                                         
 INPUT CARD>   1   5.44717800E+00   1.56284979E-01   0.00000000E+00                         
 INPUT CARD>   2   8.24547240E-01   9.04690877E-01   0.00000000E+00                         
 INPUT CARD>   S  1                                                                         
 INPUT CARD>   1   1.83191580E-01   1.00000000E+00   0.00000000E+00                         
 INPUT CARD>                                                                                
 INPUT CARD>H     1.         0.00000000         0.00000000         2.00000001               
 INPUT CARD>   S  2                                                                         
 INPUT CARD>   1   5.44717800E+00   1.56284979E-01   0.00000000E+00                         
 INPUT CARD>   2   8.24547240E-01   9.04690877E-01   0.00000000E+00                         
 INPUT CARD>   S  1                                                                         
 INPUT CARD>   1   1.83191580E-01   1.00000000E+00   0.00000000E+00                         
 INPUT CARD>                                                                                
 INPUT CARD> $END                                                                           
 INPUT CARD> $VEC                                                                           
 INPUT CARD> 1  1 2.40886664E-01 4.68755560E-01 2.40886664E-01 4.68755560E-01               
 INPUT CARD> 2  1 2.92126569E-01 5.59225642E-01-2.92126569E-01-5.59225642E-01               
 INPUT CARD> 3  1 9.07098966E-01-9.63816400E-01-9.07098966E-01 9.63816400E-01               
 INPUT CARD> 4  1-8.85652529E-01 6.67047984E-01-8.85652529E-01 6.67047984E-01               
 INPUT CARD> $END                                                                           
  500000000 WORDS OF MEMORY AVAILABLE


     RUN TITLE
     ---------
 GAMESS inp format file produced by MOKIT, nbf=4                                 

 THE POINT GROUP OF THE MOLECULE IS C1      
 THE ORDER OF THE PRINCIPAL AXIS IS     1

 ATOM      ATOMIC                      COORDINATES (BOHR)
           CHARGE         X                   Y                   Z
 H           1.0     0.0000000000        0.0000000000        0.0000000000
 H           1.0     0.0000000000        0.0000000000        3.7794519943

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                1 H          2 H     

   1 H       0.0000000    2.0000000 *
   2 H       2.0000000 *  0.0000000  

  * ... LESS THAN  3.000


     ATOMIC BASIS SET
     ----------------
 THE CONTRACTED PRIMITIVE FUNCTIONS HAVE BEEN UNNORMALIZED
 THE CONTRACTED BASIS FUNCTIONS ARE NOW NORMALIZED TO UNITY

  SHELL TYPE  PRIMITIVE        EXPONENT          CONTRACTION COEFFICIENT(S)

 H         

      1   S       1             5.4471780    0.156284978917
      1   S       2             0.8245472    0.904690876522

      2   S       3             0.1831916    1.000000000000

 H         

      3   S       4             5.4471780    0.156284978917
      3   S       5             0.8245472    0.904690876522

      4   S       6             0.1831916    1.000000000000

 TOTAL NUMBER OF BASIS SET SHELLS             =    4
 NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS =    4
 NUMBER OF ELECTRONS                          =    2
 CHARGE OF MOLECULE                           =    0
 SPIN MULTIPLICITY                            =    1
 NUMBER OF OCCUPIED ORBITALS (ALPHA)          =    1
 NUMBER OF OCCUPIED ORBITALS (BETA )          =    1
 TOTAL NUMBER OF ATOMS                        =    2
 THE NUCLEAR REPULSION ENERGY IS        0.2645886233

     ---------------------------------
     RELATIVISTIC WAVEFUNCTION OPTIONS
     ---------------------------------
     METHOD=DK       NORDER= 2       NESOC=  0      MODEQR=      1
     NRATOM=  0      QMTTOL= 1.000E-10      QRTOL= 1.000E-02
     CLIGHT=137.03598950     TAU= 3.500E+00

 INTERNAL BASIS SET UNCONTRACTION FOR RESC/DK/IOTC/LUT-IOTC:
     6 SHELLS AND      6 AOS ARE CONSTRUCTED, COMPARED TO
     4 SHELLS AND      4 AOS IN THE CONTRACTED BASIS SET


 THIS MOLECULE IS RECOGNIZED AS BEING LINEAR,
 ORBITAL LZ DEGENERACY TOLERANCE ETOLLZ= 1.00E-06

     $CONTRL OPTIONS
     ---------------
 SCFTYP=GVB          RUNTYP=ENERGY       EXETYP=RUN     
 MPLEVL=       0     CITYP =NONE         CCTYP =NONE         VBTYP =NONE    
 DFTTYP=NONE         TDDFT =NONE    
 MULT  =       1     ICHARG=       0     NZVAR =       0     COORD =UNIQUE  
 PP    =NONE         RELWFN=DK           LOCAL =NONE         NUMGRD=       F
 ISPHER=      -1     NOSYM =       1     MAXIT =     500     UNITS =ANGS    
 PLTORB=       F     MOLPLT=       F     AIMPAC=       F     FRIEND=        
 NPRINT=       7     IREST =       0     GEOM  =INPUT   
 NORMF =       0     NORMP =       0     ITOL  =      20     ICUT  =      11
 INTTYP=BEST         GRDTYP=BEST         QMTTOL= 1.0E-06

     $SYSTEM OPTIONS
     ---------------
  REPLICATED MEMORY=   500000000 WORDS (ON EVERY NODE).
 DISTRIBUTED MEMDDI=           0 MILLION WORDS IN AGGREGATE,
 MEMDDI DISTRIBUTED OVER   1 PROCESSORS IS           0 WORDS/PROCESSOR.
 TOTAL MEMORY REQUESTED ON EACH PROCESSOR=   500000000 WORDS.
 TIMLIM=      525600.00 MINUTES, OR     365.0 DAYS.
 PARALL= F  BALTYP=  DLB     KDIAG=    0  COREFL= F
 MXSEQ2=     300 MXSEQ3=     150  mem10=         0  mem22=         0

 Using Gaussians for PCM contribution to MEP, params: 0.10000E+01 0.10000E-07

          ----------------
          PROPERTIES INPUT
          ----------------

     MOMENTS            FIELD           POTENTIAL          DENSITY
 IEMOM =       1   IEFLD =       0   IEPOT =       0   IEDEN =       0
 WHERE =COMASS     WHERE =NUCLEI     WHERE =NUCLEI     WHERE =NUCLEI  
 OUTPUT=BOTH       OUTPUT=BOTH       OUTPUT=BOTH       OUTPUT=BOTH    
 IEMINT=       0   IEFINT=       0                     IEDINT=       0
                                                       MORB  =       0
          EXTRAPOLATION IN EFFECT
 ORBITAL PRINTING OPTION: NPREO=     1     4     2     1

          *************************
          ROHF-GVB INPUT PARAMETERS
          *************************

          NORB   =    2          NCO    =    0
          NPAIR  =    1          NSETO  =    0
          PAIR ORBITALS
          PAIR    1 HAS ORBS    1    2

          ----------------------------
          ROHF-GVB COUPLING PARAMETERS
          ----------------------------

          F VECTOR (OCCUPANCIES)
    1   0.9529411765
    2   0.0470588235
          ALPHA COUPLING COEFFICEINTS

             1           2

    1    0.9529412
    2    0.0000000   0.0470588
           BETA COUPLING COEFFICIENTS

             1           2

    1    0.0000000
    2   -0.2117647   0.0000000

          NATURAL ORBITAL COEFFICIENTS
          N.O.                PAIR CICOEF-S
            1        0.9761870602    -0.2169304578

     -------------------------------
     INTEGRAL TRANSFORMATION OPTIONS
     -------------------------------
     NWORD  =            0
     CUTOFF = 1.0E-09     MPTRAN =       0
     DIRTRF =       T     AOINTS =DUP     

          ----------------------
          INTEGRAL INPUT OPTIONS
          ----------------------
 NOPK  =       1 NORDER=       0 SCHWRZ=       T

     ------------------------------------------
     THE POINT GROUP IS C1 , NAXIS= 1, ORDER= 1
     ------------------------------------------

     DIMENSIONS OF THE SYMMETRY SUBSPACES ARE
 A   =    4

 ..... DONE SETTING UP THE RUN .....
 STEP CPU TIME =     0.01 TOTAL CPU TIME =          0.0 (      0.0 MIN)
 TOTAL WALL CLOCK TIME=          0.0 SECONDS, CPU UTILIZATION IS    50.00%

          ********************
          1 ELECTRON INTEGRALS
          ********************
 TIME TO DO ORDINARY INTEGRALS=      0.00

 ------------------------------------------------------------------
 2ND ORDER DOUGLAS-KROLL-HESS RELATIVISTIC ONE-ELECTRON HAMILTONIAN
            CODED BY TAKAHITO NAKAJIMA AND DMITRI FEDOROV
 ------------------------------------------------------------------

 TIME TO DO  RELATIVISTIC INTS=      0.00
 ...... END OF ONE-ELECTRON INTEGRALS ......
 STEP CPU TIME =     0.00 TOTAL CPU TIME =          0.0 (      0.0 MIN)
 TOTAL WALL CLOCK TIME=          0.0 SECONDS, CPU UTILIZATION IS    33.33%

          -------------
          GUESS OPTIONS
          -------------
          GUESS =MOREAD            NORB  =       4          NORDER=       0
          MIX   =       F          PRTMO =       F          PUNMO =       F
          TOLZ  = 1.0E-08          TOLE  = 1.0E-05
          SYMDEN=       F          PURIFY=       F

 INITIAL GUESS ORBITALS GENERATED BY MOREAD   ROUTINE.

 SYMMETRIES FOR INITIAL GUESS ORBITALS FOLLOW.   BOTH SET(S).
     2 ORBITALS ARE OCCUPIED (    0 CORE ORBITALS).
     1=A        2=A        3=A        4=A   
 ...... END OF INITIAL ORBITAL SELECTION ......
 STEP CPU TIME =     0.00 TOTAL CPU TIME =          0.0 (      0.0 MIN)
 TOTAL WALL CLOCK TIME=          0.0 SECONDS, CPU UTILIZATION IS    33.33%

                    ----------------------
                    AO INTEGRAL TECHNOLOGY
                    ----------------------
     S,P,L SHELL ROTATED AXIS INTEGRALS, REPROGRAMMED BY
        KAZUYA ISHIMURA (IMS) AND JOSE SIERRA (SYNSTAR).
     S,P,D,L SHELL ROTATED AXIS INTEGRALS PROGRAMMED BY
        KAZUYA ISHIMURA (INSTITUTE FOR MOLECULAR SCIENCE).
     S,P,D,F,G SHELL TO TOTAL QUARTET ANGULAR MOMENTUM SUM 5,
        ERIC PROGRAM BY GRAHAM FLETCHER (ELORET AND NASA ADVANCED
        SUPERCOMPUTING DIVISION, AMES RESEARCH CENTER).
     S,P,D,F,G,L SHELL GENERAL RYS QUADRATURE PROGRAMMED BY
        MICHEL DUPUIS (PACIFIC NORTHWEST NATIONAL LABORATORY).

          --------------------
          2 ELECTRON INTEGRALS
          --------------------

 DIRECT SCF METHOD SKIPS INTEGRAL STORAGE ON DISK.
 DIRECT TRANSFORMATION SKIPS AO INTEGRAL STORAGE ON DISK.
  ...... END OF TWO-ELECTRON INTEGRALS .....
 STEP CPU TIME =     0.01 TOTAL CPU TIME =          0.0 (      0.0 MIN)
 TOTAL WALL CLOCK TIME=          0.1 SECONDS, CPU UTILIZATION IS    33.33%

          ------------------------
          ROHF-GVB SCF CALCULATION
          ------------------------
 GVB STEP WILL USE     44720 WORDS OF MEMORY.

     MAXIT= 500   NPUNCH= 2   SQCDF TOL=2.0000E-05
     NUCLEAR ENERGY=        0.2645886233
     EXTRAP=T   DAMP=F   SHIFT=F   RSTRCT=F   DIIS=F  SOSCF=F

 DIRECT SCF CALCULATION, SCHWRZ=T   FDIFF=F
 SCHWARZ INEQUALITY OVERHEAD:        10 INTEGRALS, T=        0.00

                                                                           NONZERO    BLOCKS
 ITER EX     TOTAL ENERGY       E CHANGE        SQCDF       DIIS ERROR      INTEGRALS   SKIPPED
   0  0       -1.009649214    -1.009649214   0.019645598   0.000000000             55         0
   1  1       -1.009954864    -0.000305651   0.001079952   0.000000000             55         0
   2  2       -1.009955793    -0.000000929   0.000060810   0.000000000             55         0
   3  3       -1.009955795    -0.000000003   0.000002666   0.000000000             55         0
   4  4       -1.009955795    -0.000000000   0.000000697   0.000000000             55         0

          -----------------
          DENSITY CONVERGED
          -----------------

 FINAL GVB ENERGY IS       -1.0099557954 AFTER   4 ITERATIONS

          ----------------------------
          SCF STATISTICS PER ITERATION
          ----------------------------
          NUMBER OF INTEGRAL PASSES       1
          FOCK FORMATION TIME         0.000
          GEMINAL OPT TIME            0.000
          MIXORB  OPT TIME            0.000
          OCBSE   OPT TIME            0.000

          ----------------
          PAIR INFORMATION
          ----------------
      ORBITAL   CI COEFFICIENTS    OCCUPATION NUMBERS   GVB      ENERGY
 PAIR 1   2     ORB 1     ORB 2      ORB 1     ORB 2    OVERLAP  LOWERING
  1   1   2    0.874860 -0.484376   1.53076   0.46924   0.28728  -0.10838

     THE MAXIMUM LAGRANGIAN ASYMMETRY IS 0.0000000E+00

 LZ VALUE ANALYSIS FOR THE MOS
 ----------------------------------------
 MO     1 (    1) HAS LZ(WEIGHT)= 0.00(100.0%) 
 MO     2 (    2) HAS LZ(WEIGHT)= 0.00(100.0%) 
 MO     3 (    3) HAS LZ(WEIGHT)= 0.00(100.0%) 
 MO     4 (    4) HAS LZ(WEIGHT)= 0.00(100.0%) 

          ------------
          EIGENVECTORS
          ------------

                      1          2          3          4
                   -0.2744    -0.0319     1.1530     1.1943
                     A          A          A          A   
    1  H  1  S    0.235981   0.308212   0.886972  -0.901761
    2  H  1  S    0.472440   0.541997  -0.664443   0.973609
    3  H  2  S    0.235981  -0.308212   0.886972   0.901761
    4  H  2  S    0.472440  -0.541997  -0.664443  -0.973609

          -----------
          GI ORBITALS
          -----------


                    PAIR   1

                      1          2

    1  H  1  S    0.373310  -0.005331
    2  H  1  S    0.702576  -0.055476
    3  H  2  S    0.005331  -0.373310
    4  H  2  S    0.055476  -0.702576

 ... END OF ROHF-GVB SCF CALCULATION ...
 STEP CPU TIME =     0.00 TOTAL CPU TIME =          0.0 (      0.0 MIN)
 TOTAL WALL CLOCK TIME=          0.1 SECONDS, CPU UTILIZATION IS    28.57%

     ----------------------------------------------------------------
     PROPERTY VALUES FOR THE GVB   SELF-CONSISTENT FIELD WAVEFUNCTION
     ----------------------------------------------------------------

          -----------------
          ENERGY COMPONENTS
          -----------------

         WAVEFUNCTION NORMALIZATION =       1.0000000000

                ONE ELECTRON ENERGY =      -1.5601229474
                TWO ELECTRON ENERGY =       0.2855785287
           NUCLEAR REPULSION ENERGY =       0.2645886233
                                      ------------------
                       TOTAL ENERGY =      -1.0099557954

 ELECTRON-ELECTRON POTENTIAL ENERGY =       0.2855785287
  NUCLEUS-ELECTRON POTENTIAL ENERGY =      -2.5156481225
   NUCLEUS-NUCLEUS POTENTIAL ENERGY =       0.2645886233
                                      ------------------
             TOTAL POTENTIAL ENERGY =      -1.9654809706
               TOTAL KINETIC ENERGY =       0.9555251752
                 VIRIAL RATIO (V/T) =       2.0569640881

          ---------------------------------------
          MULLIKEN AND LOWDIN POPULATION ANALYSES
          ---------------------------------------

     ATOMIC MULLIKEN POPULATION IN EACH MOLECULAR ORBITAL

                      1          2

                  1.530759   0.469241

    1             0.765379   0.234621
    2             0.765379   0.234621

               ----- POPULATIONS IN EACH AO -----
                             MULLIKEN      LOWDIN
              1  H  1  S      0.29766     0.33681
              2  H  1  S      0.70234     0.66319
              3  H  2  S      0.29766     0.33681
              4  H  2  S      0.70234     0.66319

          ----- MULLIKEN ATOMIC OVERLAP POPULATIONS -----
          (OFF-DIAGONAL ELEMENTS NEED TO BE MULTIPLIED BY 2)

             1           2

    1    0.9310459
    2    0.0689541   0.9310459

          TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS
       ATOM         MULL.POP.    CHARGE          LOW.POP.     CHARGE
    1 H             1.000000   -0.000000         1.000000   -0.000000
    2 H             1.000000   -0.000000         1.000000    0.000000

          ---------------------
          ELECTROSTATIC MOMENTS
          ---------------------

 POINT   1           X           Y           Z (BOHR)    CHARGE
                 0.000000    0.000000    1.889726        0.00 (A.U.)
         DX          DY          DZ         /D/  (DEBYE)
     0.000000    0.000000    0.000000    0.000000
 ...... END OF PROPERTY EVALUATION ......
 STEP CPU TIME =     0.01 TOTAL CPU TIME =          0.0 (      0.0 MIN)
 TOTAL WALL CLOCK TIME=          0.1 SECONDS, CPU UTILIZATION IS    37.50%
                580000  WORDS OF DYNAMIC MEMORY USED
 EXECUTION OF GAMESS TERMINATED NORMALLY Wed Jan 20 23:21:34 2021
 DDI: 263640 bytes (0.3 MB / 0 MWords) used by master data server.

 ----------------------------------------
 CPU timing information for all processes
 ========================================
 0: 0.13 + 0.22 = 0.36
 ----------------------------------------
 ddikick.x: exited gracefully.
unset echo
----- accounting info -----
Files used on the master node xn01 were:
-rw-rw-r-- 1 srwang srwang    1360 Jan 20 23:21 /scratch/scr/srwang/gamess/cas_uhf_uno_asrot2gvb1.dat
-rw-rw-r-- 1 srwang srwang    1036 Jan 20 23:21 /scratch/scr/srwang/gamess/cas_uhf_uno_asrot2gvb1.F05
-rw-rw-r-- 1 srwang srwang 4482640 Jan 20 23:21 /scratch/scr/srwang/gamess/cas_uhf_uno_asrot2gvb1.F10
-rw-rw-r-- 1 srwang srwang     616 Jan 20 23:21 /scratch/scr/srwang/gamess/cas_uhf_uno_asrot2gvb1.F23
ls: No match.
ls: No match.
ls: No match.
Wed Jan 20 23:21:37 CST 2021
0.211u 0.130s 0:03.43 9.9%	0+0k 0+8io 0pf+0w
